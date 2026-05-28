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
        keep = np.array([np.sum(m ** 2) > 1e-14 for m in modes])
        self.normal_modes = modes[keep]
        
        self.masses = np.array([_get_mass(a) for a in atoms], dtype=float)
        self.geom = np.asarray(geom, dtype=float)
        self.atoms = atoms
        self.n_atoms = len(atoms)
        self.n_modes = self.normal_modes.shape[0]

    def localize_mode(self, mode: int) -> np.ndarray:
        real_displacements = self.normal_modes[mode] / np.sqrt(self.masses[:, None])
        distances = np.linalg.norm(real_displacements, axis=1)
        total_distance = np.sum(distances)
        if total_distance == 0:
            return np.zeros_like(distances)
            
        return (distances / total_distance) * 100.0

    def localize_mode_fragment(
        self,
        mode: int,
        fragments: list[list[int]],
    ) -> dict[str, float]:
        atomic = self.localize_mode(mode)
        return {
            f"Frag{i+1}": float(np.sum(atomic[indices]))
            for i, indices in enumerate(fragments)
        }

    def displace_geometry(
        self, mode: int, scale: float = 0.3
    ) -> np.ndarray:
        real_displacements = self.normal_modes[mode] / np.sqrt(self.masses[:, None])
        return self.geom + scale * real_displacements

    @staticmethod
    def classify_negative_freqs(
        freqs: np.ndarray,
        threshold: float = 20.0,
    ) -> tuple[list[int], list[int]]:
        neg = np.where(freqs < 0)[0]
        significant = [int(i) for i in neg if abs(freqs[i]) > threshold]
        noise = [int(i) for i in neg if abs(freqs[i]) <= threshold]
        return significant, noise

    @staticmethod
    def derive_internals_from_fragments(
        fragments: list[list[int]],
    ) -> list[list[int]]:
        internals: list[list[int]] = []
        for frag in fragments:
            n = len(frag)
            if n >= 2:
                for i in range(n - 1):
                    internals.append([frag[i], frag[i + 1]])
            if n >= 3:
                for i in range(n - 2):
                    internals.append([frag[i], frag[i + 1], frag[i + 2]])
            if n >= 4:
                for i in range(n - 3):
                    internals.append([frag[i], frag[i + 1], frag[i + 2], frag[i + 3]])
        return internals

    @staticmethod
    def derive_internals_from_connectivity(
        symbols: tuple[str, ...],
        positions: np.ndarray,
        scale: float = 1.3,
    ) -> list[list[int]]:
        n = len(symbols)
        if n < 2:
            return []

        bonds: set[tuple[int, int]] = set()
        for i in range(n):
            for j in range(i + 1, n):
                r_cov = covalent_radii[atomic_numbers[symbols[i]]] + covalent_radii[atomic_numbers[symbols[j]]]
                d = np.linalg.norm(positions[i] - positions[j])
                if d < scale * r_cov:
                    bonds.add((i, j))

        adj: list[set[int]] = [set() for _ in range(n)]
        for i, j in bonds:
            adj[i].add(j)
            adj[j].add(i)

        internals: list[list[int]] = []

        for i, j in bonds:
            internals.append([i, j])

        for j in range(n):
            neighbors = sorted(adj[j])
            for p in range(len(neighbors)):
                for q in range(p + 1, len(neighbors)):
                    i, k = neighbors[p], neighbors[q]
                    internals.append([i, j, k])

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
                if diff > 180:
                    diff -= 360
                elif diff < -180:
                    diff += 360
                comp = diff / (2 * step_size)
            else:
                continue

            label = self._format_internal_label(idx)
            raw.append((label, comp, idx))

        return raw

    def _arc_weight(self, comp: float, idx: list[int]) -> float:
        if len(idx) == 2:
            return abs(comp)
        deg2rad = np.pi / 180.0
        atoms = ASE_Atoms(
            symbols="".join(self.atoms),
            positions=self.geom.copy(),
        )
        mean_r = self._mean_bond_length(atoms, idx)
        return abs(comp * deg2rad) * mean_r

    def project_on_internals(
        self,
        mode: int,
        internals: list[list[int]],
        normalization_internals: list[list[int]] | None = None,
        step_size: float = 1e-3,
    ) -> list[tuple[str, float, float]]:
        norm_set = normalization_internals if normalization_internals is not None else internals
        raw_norm = self._compute_raw_components(mode, norm_set, step_size)
        if not raw_norm:
            return []

        arc_weights = [self._arc_weight(comp, idx) for _, comp, idx in raw_norm]
        total = sum(arc_weights)
        if total < 1e-30:
            return [(self._format_internal_label(idx), 0.0, 0.0) for idx in internals]

        pct_map: dict[str, float] = {}
        for (label, _, _), arc in zip(raw_norm, arc_weights):
            pct_map[label] = float(arc / total * 100.0)

        result: list[tuple[str, float, float]] = []
        for idx in internals:
            label = self._format_internal_label(idx)
            comp = 0.0
            for lbl, c, i in raw_norm:
                if lbl == label:
                    comp = c
                    break
            result.append((label, comp, pct_map.get(label, 0.0)))

        return result

    def imag_mode_summary(
        self,
        mode: int,
        fragments: list[list[int]] | None = None,
    ) -> dict:
        atomic = self.localize_mode(mode)
        top_idx = np.argsort(atomic)[-5:][::-1]
        top_atoms = [
            (int(idx), self.atoms[idx], float(atomic[idx]))
            for idx in top_idx
            if atomic[idx] > 1.0
        ]
        result = {
            "mode": mode,
            "top_atoms": top_atoms,
        }
        if fragments:
            result["fragments"] = self.localize_mode_fragment(mode, fragments)
        return result

    def get_dominant_internals(
        self, 
        mode: int, 
        threshold_pct: float = 1.0,
        step_size: float = 1e-3
    ) -> dict[str, list[tuple[str, float, float]]]:
        all_internals = self.derive_internals_from_connectivity(self.atoms, self.geom)
        raw_norm = self._compute_raw_components(mode, all_internals, step_size)
        
        results = {"Bonds": [], "Angles": [], "Dihedrals": []}
        if not raw_norm:
            return results

        arc_weights = [self._arc_weight(comp, idx) for _, comp, idx in raw_norm]
        total = sum(arc_weights)
        
        if total < 1e-30:
            return results

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
                    
        for key in results:
            results[key].sort(key=lambda x: x[2], reverse=True)
            
        return results

    def print_mode_decomposition(self, mode: int, threshold_pct: float = 1.0):
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