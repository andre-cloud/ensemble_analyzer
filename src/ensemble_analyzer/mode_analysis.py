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
        # Filtra i modi traslazionali e rotazionali
        keep = np.array([np.sum(m ** 2) > 1e-14 for m in modes])
        # Conserva la matrice originaria (pesata per la massa)
        self.normal_modes = modes[keep]
        
        self.masses = np.array([_get_mass(a) for a in atoms], dtype=float)
        self.geom = np.asarray(geom, dtype=float)
        self.atoms = atoms
        self.n_atoms = len(atoms)
        self.n_modes = self.normal_modes.shape[0]

    def localize_mode(self, mode: int) -> np.ndarray:
        # Ricava gli spostamenti cartesiani puri (geometria)
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
        # Nota: se per la visualizzazione vuoi usare gli spostamenti cartesiani reali, 
        # dovresti dividere per la massa anche qui. L'ho lasciato come l'originale
        # per non alterare le tue pipeline esterne.
        return self.geom + scale * self.normal_modes[mode]

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
        """Derive bond/angle/dihedral definitions from loc_freq fragments.

        For each fragment, generates internal coordinates between consecutive
        atoms: bonds (pairs), angles (triples), dihedrals (quadruples).

        Args:
            fragments: List of atom groups, e.g. [[0,1,2,3], [4,5,6,7]]

        Returns:
            List of internal definitions, each a list of 2, 3, or 4 atom indices.
        """
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
        """Generate all bond/angle/dihedral definitions from molecular connectivity.

        Uses covalent radii to determine bonds, then generates angles and
        dihedrals from consecutive bonds in the connectivity graph.

        Args:
            symbols: Tuple of element symbols.
            positions: Geometry array of shape (n_atoms, 3).
            scale: Multiplier for covalent radii cutoff.

        Returns:
            List of internal definitions, each a list of 2, 3, or 4 atom indices.
        """
        n = len(symbols)
        if n < 2:
            return []

        # Determine bonds via covalent radii
        bonds: set[tuple[int, int]] = set()
        for i in range(n):
            for j in range(i + 1, n):
                # CORREZIONE BUG ASE: rimosso il "- 1" dagli indici dei raggi covalenti
                r_cov = covalent_radii[atomic_numbers[symbols[i]]] + covalent_radii[atomic_numbers[symbols[j]]]
                d = np.linalg.norm(positions[i] - positions[j])
                if d < scale * r_cov:
                    bonds.add((i, j))

        # Build adjacency list
        adj: list[set[int]] = [set() for _ in range(n)]
        for i, j in bonds:
            adj[i].add(j)
            adj[j].add(i)

        internals: list[list[int]] = []

        # Bonds
        for i, j in bonds:
            internals.append([i, j])

        # Angles: i-j-k where i!=k and both bonded to j
        for j in range(n):
            neighbors = sorted(adj[j])
            for p in range(len(neighbors)):
                for q in range(p + 1, len(neighbors)):
                    i, k = neighbors[p], neighbors[q]
                    internals.append([i, j, k])

        # Dihedrals: i-j-k-l where i bonded to j, j to k, k to l
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
        """Format label for an internal coordinate.

        Returns e.g. "B(C0,C1)", "A(C0,C1,C2)", "D(C0,C1,C2,C3)".
        """
        prefixes = {2: "B", 3: "A", 4: "D"}
        prefix = prefixes.get(len(indices), "?")
        symbols = ",".join(f"{self.atoms[i]}{i}" for i in indices)
        return f"{prefix}({symbols})"

    def _mean_bond_length(self, atoms: ASE_Atoms, indices: list[int]) -> float:
        """Geometric mean of consecutive bond lengths for a set of atom indices."""
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
        """Compute raw finite-difference components for a list of internals.

        Creates atoms_plus/atoms_minus once and evaluates all internals
        on both geometries for efficiency.

        Returns list of (label, component, idx) tuples.
        """
        atoms = ASE_Atoms(
            symbols="".join(self.atoms),
            positions=self.geom.copy(),
        )
        
        # CORREZIONE FISICA: Rimuove il peso della massa per avere gli spostamenti geometrici puri
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
        """Convert a raw component to arc-length in Å."""
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
        """Project a normal mode onto internal coordinates.

        Uses central finite differences to compute the component of the mode
        along each internal coordinate (bond, angle, dihedral).
        Also returns the percentage contribution based on arc-length normalisation
        (angles/dihedrals converted to Å via mean bond lengths).

        When normalization_internals is provided, percentages are computed over
        the full molecular internal coordinate space, making bond and angle
        contributions directly comparable.

        Args:
            mode: Index of the normal mode to project.
            internals: List of internal coordinate definitions to report.
            normalization_internals: Full set of internals for percentage
                normalisation. If None, uses internals only.
            step_size: Small displacement for finite difference in Å.

        Returns:
            List of (label, raw_component, percentage) tuples.
        """
        norm_set = normalization_internals if normalization_internals is not None else internals
        raw_norm = self._compute_raw_components(mode, norm_set, step_size)
        if not raw_norm:
            return []

        arc_weights = [self._arc_weight(comp, idx) for _, comp, idx in raw_norm]
        total = sum(arc_weights)
        if total < 1e-30:
            return [(self._format_internal_label(idx), 0.0, 0.0) for idx in internals]

        # Build map: label → percentage
        pct_map: dict[str, float] = {}
        for (label, _, _), arc in zip(raw_norm, arc_weights):
            pct_map[label] = float(arc / total * 100.0)

        # Return results only for requested internals
        result: list[tuple[str, float, float]] = []
        for idx in internals:
            label = self._format_internal_label(idx)
            comp = 0.0
            # Find raw component
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