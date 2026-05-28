import numpy as np
from ase import Atoms as ASE_Atoms
from ase.data import covalent_radii, atomic_numbers

_ATOMIC_MASSES: dict[str, float] = {}

def _get_mass(symbol: str) -> float:
    m = _ATOMIC_MASSES.get(symbol)
    if m is None:
        m = ASE_Atoms(symbol).get_masses()[0]
        _ATOMIC_MASSES[symbol] = m
    return m


class NormalModeAnalyzer:
    def __init__(
        self,
        normal_modes: np.ndarray,
        geom: np.ndarray,
        atoms: tuple[str, ...],
    ):
        modes = np.asarray(normal_modes, dtype=float)
        # Filtro per scartare traslazioni/rotazioni
        keep = np.array([np.sum(m ** 2) > 1e-14 for m in modes])
        # Manteniamo i vettori pesati per la massa per coerenza termodinamica
        self.normal_modes = modes[keep]
        
        self.masses = np.array([_get_mass(a) for a in atoms], dtype=float)
        self.geom = np.asarray(geom, dtype=float)
        self.atoms = atoms
        self.n_atoms = len(atoms)
        self.n_modes = self.normal_modes.shape[0]

    def localize_mode(self, mode: int) -> np.ndarray:
        """Calcola la percentuale di spostamento geometrico cartesiano per atomo."""
        real_displacements = self.normal_modes[mode] / np.sqrt(self.masses[:, None])
        distances = np.linalg.norm(real_displacements, axis=1)
        total_distance = np.sum(distances)
        if total_distance == 0:
            return np.zeros_like(distances)
        return (distances / total_distance) * 100.0

    @staticmethod
    def derive_internals_from_connectivity(
        symbols: tuple[str, ...],
        positions: np.ndarray,
        scale: float = 1.3,
    ) -> list[list[int]]:
        """Genera legami, angoli e diedri dalla topologia spaziale usando i raggi covalenti."""
        n = len(symbols)
        if n < 2:
            return []

        bonds: set[tuple[int, int]] = set()
        for i in range(n):
            for j in range(i + 1, n):
                # Bug ASE risolto: indice corretto per covalent_radii
                r_cov = covalent_radii[atomic_numbers[symbols[i]]] + covalent_radii[atomic_numbers[symbols[j]]]
                d = np.linalg.norm(positions[i] - positions[j])
                if d < scale * r_cov:
                    bonds.add((i, j))

        adj: list[set[int]] = [set() for _ in range(n)]
        for i, j in bonds:
            adj[i].add(j)
            adj[j].add(i)

        internals: list[list[int]] = []

        # Legami (2 atomi)
        for i, j in bonds:
            internals.append([i, j])

        # Angoli (3 atomi)
        for j in range(n):
            neighbors = sorted(adj[j])
            for p in range(len(neighbors)):
                for q in range(p + 1, len(neighbors)):
                    i, k = neighbors[p], neighbors[q]
                    internals.append([i, j, k])

        # Diedri (4 atomi)
        for j in range(n):
            for k in adj[j]:
                if k <= j:
                    continue
                for i in adj[j]:
                    if i == k:
                        continue
                    for l in adj[k]:
                        if l == j or l == i:
                            continue
                        internals.append([i, j, k, l])

        return internals

    def _format_internal_label(self, indices: list[int]) -> str:
        prefixes = {2: "B", 3: "A", 4: "D"}
        prefix = prefixes.get(len(indices), "?")
        symbols = ",".join(f"{self.atoms[i]}{i}" for i in indices)
        return f"{prefix}({symbols})"

    def _mean_bond_length(self, atoms: ASE_Atoms, indices: list[int]) -> float:
        bonds = []
        for i in range(len(indices) - 1):
            bonds.append(atoms.get_distance(indices[i], indices[i + 1]))
        prod = 1.0
        for b in bonds:
            prod *= b
        return prod ** (1.0 / len(bonds)) if bonds else 1.0

    def _compute_raw_components(
        self,
        mode: int,
        internals: list[list[int]],
        step_size: float,
    ) -> list[tuple[str, float, list[int]]]:
        atoms = ASE_Atoms(
            symbols="".join(self.atoms),
            positions=self.geom.copy(),
        )
        
        # Bug Fisico risolto: Estrae gli spostamenti geometrici puri (rimuove la massa)
        real_displacements = self.normal_modes[mode] / np.sqrt(self.masses[:, None])
        displacement = real_displacements * step_size

        atoms_plus = atoms.copy()
        atoms_minus = atoms.copy()
        atoms_plus.positions += displacement
        atoms_minus.positions -= displacement

        raw: list[tuple[str, float, list[int]]] = []
        for idx in internals:
            n = len(idx)
            if n == 2:
                q_plus = atoms_plus.get_distance(*idx)
                q_minus = atoms_minus.get_distance(*idx)
                comp = (q_plus - q_minus) / (2 * step_size)
            elif n == 3:
                q_plus = atoms_plus.get_angle(*idx)
                q_minus = atoms_minus.get_angle(*idx)
                comp = (q_plus - q_minus) / (2 * step_size)
            elif n == 4:
                q_plus = atoms_plus.get_dihedral(*idx)
                q_minus = atoms_minus.get_dihedral(*idx)
                diff = q_plus - q_minus
                if diff > 180: diff -= 360
                elif diff < -180: diff += 360
                comp = diff / (2 * step_size)
            else:
                continue

            label = self._format_internal_label(idx)
            raw.append((label, comp, idx))

        return raw

    def _arc_weight(self, comp: float, idx: list[int]) -> float:
        """Converte le variazioni angolari in lunghezza d'arco equivalente (Ångström)."""
        if len(idx) == 2:
            return abs(comp)
        deg2rad = np.pi / 180.0
        atoms = ASE_Atoms(
            symbols="".join(self.atoms),
            positions=self.geom.copy(),
        )
        mean_r = self._mean_bond_length(atoms, idx)
        return abs(comp * deg2rad) * mean_r

    def get_dominant_internals(
        self, 
        mode: int, 
        threshold_pct: float = 1.0,
        step_size: float = 1e-3
    ) -> dict[str, list[tuple[str, float, float]]]:
        """
        Deriva automaticamente tutte le coordinate, calcola le proiezioni, 
        e restituisce SOLO quelle che superano la soglia percentuale (es. > 1%).
        Ideale per scartare il rumore a 0.0% e isolare i moti chiave.
        
        Ritorna un dizionario diviso per categoria: {'Bonds': [...], 'Angles': [...], 'Dihedrals': [...]}
        """
        # 1. Genera tutta la topologia corretta
        all_internals = self.derive_internals_from_connectivity(self.atoms, self.geom)
        
        # 2. Calcola componenti raw e pesi ad arco
        raw_norm = self._compute_raw_components(mode, all_internals, step_size)
        if not raw_norm:
            return {"Bonds": [], "Angles": [], "Dihedrals": []}

        arc_weights = [self._arc_weight(comp, idx) for _, comp, idx in raw_norm]
        total = sum(arc_weights)
        
        results = {"Bonds": [], "Angles": [], "Dihedrals": []}
        
        if total < 1e-30:
            return results

        # 3. Filtra e smista i risultati
        for (label, comp, idx), arc in zip(raw_norm, arc_weights):
            pct = float((arc / total) * 100.0)
            
            if pct >= threshold_pct:
                entry = (label, comp, pct)
                if len(idx) == 2:
                    results["Bonds"].append(entry)
                elif len(idx) == 3:
                    results["Angles"].append(entry)
                elif len(idx) == 4:
                    results["Dihedrals"].append(entry)
                    
        # Ordina ogni lista per percentuale decrescente
        for key in results:
            results[key].sort(key=lambda x: x[2], reverse=True)
            
        return results

    def print_mode_decomposition(self, mode: int, threshold_pct: float = 1.0):
        """Metodo di utilità per stampare la decomposizione in modo leggibile (come una tabella)."""
        dominant = self.get_dominant_internals(mode, threshold_pct)
        
        print(f"=== Decomposizione Modo {mode} (Top Contributi > {threshold_pct}%) ===")
        found_any = False
        
        for category, items in dominant.items():
            if items:
                found_any = True
                print(f"\n-- {category} --")
                for label, comp, pct in items:
                    print(f"{label:20s} | Variazione: {comp:>10.4f} | Peso: {pct:>6.2f} %")
                    
        if not found_any:
            print("\nNessuna coordinata supera la soglia di taglio. Il moto potrebbe essere una traslazione/rotazione rigida globale.")
        print("==========================================================")