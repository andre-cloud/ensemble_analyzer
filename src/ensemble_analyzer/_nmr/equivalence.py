"""Chemically equivalent nuclei for NMR.

Copied from FROGChem (paclab-lib: ``frogchem/symmetry.py``, ``graph_theory.py``,
``math_geometry.py``, ``nmr.get_nmr_equivalencies``); keep in sync. Changes: atomic radii
and masses from ``ase.data`` instead of ``periodictable``; the input is ``(symbols, coords)``.
"""

from collections import defaultdict
from typing import Any, Dict, List, Tuple

import numpy as np
from ase.data import atomic_masses, atomic_numbers, covalent_radii
from scipy.spatial.distance import cdist


def get_atomic_weights(symbol: str) -> float:
    return float(atomic_masses[atomic_numbers[symbol]])


def get_atomic_radius(symbol: str) -> float:
    return float(covalent_radii[atomic_numbers[symbol]])


def normalize(v: np.ndarray) -> np.ndarray:
    """Normalize a vector to unit length, returning the original if norm < 1e-10."""
    norm = np.linalg.norm(v)
    return v if norm < 1e-10 else v / norm


def get_rotation_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
    """Generate a 3x3 rotation matrix for rotation around an axis by theta radians."""
    axis = normalize(axis)
    a, b, c = axis
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    m = 1 - cos_t
    return np.array([
        [a*a*m + cos_t, a*b*m - c*sin_t, a*c*m + b*sin_t],
        [b*a*m + c*sin_t, b*b*m + cos_t, b*c*m - a*sin_t],
        [c*a*m - b*sin_t, c*b*m + a*sin_t, c*c*m + cos_t]
    ])


def get_reflection_matrix(normal: np.ndarray) -> np.ndarray:
    """Generate a 3x3 Householder reflection matrix across the plane perpendicular to normal."""
    normal = normalize(normal)
    return np.eye(3) - 2 * np.outer(normal, normal)


def calculate_center_of_mass(symbols: np.ndarray, positions: np.ndarray) -> Tuple[np.ndarray, float]:
    """Calculate center of mass and total molecular mass."""
    masses = np.array([get_atomic_weights(s) for s in symbols])
    total_mass = np.sum(masses)
    com = np.sum(positions * masses[:, np.newaxis], axis=0) / total_mass
    return com, total_mass


def _inertia_tensor_from_masses(masses: np.ndarray, coords: np.ndarray) -> np.ndarray:
    """Compute 3x3 moment of inertia tensor from masses and centered coordinates."""
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    Ixx = np.sum(masses * (y**2 + z**2))
    Iyy = np.sum(masses * (x**2 + z**2))
    Izz = np.sum(masses * (x**2 + y**2))
    Ixy = -np.sum(masses * x * y)
    Ixz = -np.sum(masses * x * z)
    Iyz = -np.sum(masses * y * z)
    return np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])



def create_graph_matrix(
    symbols: List[str],
    positions: np.ndarray,
    radii: Dict[str, float],
    buffer: float = 0.3
) -> np.ndarray:
    """Generates the connectivity graph matrix of the system.

    Two atoms are considered connected if their distance is less than the sum
    of their radii plus the buffer.

    Args:
        symbols (List[str]): List of atomic symbols.
        positions (np.ndarray): Array of shape (N, 3) with atomic coordinates.
        radii (Dict[str, float]): Dictionary mapping symbols to atomic radii.
        buffer (float, optional): Buffer distance in Angstrom. Defaults to 0.3.

    Returns:
        np.ndarray: Adjacency matrix (N, N) where 1 implies connection.
    """
    n_atoms = len(symbols)
    # ponytail: cdist replaces N*(N-1)/2 individual norm calls
    dist_matrix = cdist(positions, positions)

    r = np.array([radii.get(s, 1.5) for s in symbols])
    threshold = r[:, None] + r[None, :] + buffer

    graph = (dist_matrix < threshold).astype(int)
    np.fill_diagonal(graph, 0)
    return graph



class DisjointSet:
    """Disjoint-set (Union-Find) data structure with path compression and union by rank."""

    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, i: int) -> int:
        """Find the representative element of the set containing i with path compression."""
        if self.parent[i] == i:
            return i
        self.parent[i] = self.find(self.parent[i])
        return self.parent[i]

    def union(self, i: int, j: int) -> bool:
        """Merge the sets containing i and j. Returns True if merged, False if already in the same set."""
        root_i = self.find(i)
        root_j = self.find(j)
        if root_i != root_j:
            if self.rank[root_i] < self.rank[root_j]:
                root_i, root_j = root_j, root_i
            self.parent[root_j] = root_i
            if self.rank[root_i] == self.rank[root_j]:
                self.rank[root_i] += 1
            return True
        return False

    def get_groups(self) -> Dict[int, List[int]]:
        """Return a dictionary mapping root representatives to member indices."""
        groups = defaultdict(list)
        for i in range(len(self.parent)):
            groups[self.find(i)].append(i)
        return dict(groups)


class PointGroupDetector:
    """Point-group detector.

    exclude_h=False (default): the true molecular point group (all atoms).
    exclude_h=True: the point group of the heavy-atom framework, centred on the
    framework's own centre of mass (useful to tolerate noisy/rotameric H). If the
    framework has fewer than two atoms, all atoms are used instead.
    """
    def __init__(self, symbols: List[str], positions: np.ndarray, tolerance: float = 0.1, exclude_h: bool = False):
        self.symbols = np.array(symbols)
        positions = np.asarray(positions, dtype=float)
        self.tol = tolerance
        self.n_atoms = len(symbols)
        self.exclude_h = exclude_h

        if self.exclude_h:
            self.mask = np.array([s.upper() != 'H' for s in self.symbols])
            if np.sum(self.mask) < 2:
                self.mask = np.ones(self.n_atoms, dtype=bool)
        else:
            self.mask = np.ones(self.n_atoms, dtype=bool)

        # Symmetry elements pass through the centre of mass of the analysed atoms
        com, _ = calculate_center_of_mass(self.symbols[self.mask], positions[self.mask])
        self.coords = positions - com

        self.valid_ops = []
        self.found_axes = []
        self.found_planes = []
        self.found_improper = []
        self.has_inversion = False
        self.max_op_rmsd = 0.0

        self._analyze_inertia()
        self.candidate_axes = self._generate_candidate_axes()

    def _analyze_inertia(self):
        syms_active = self.symbols[self.mask]
        coords_active = self.coords[self.mask]
        masses = np.array([get_atomic_weights(s) for s in syms_active])
        tensor = _inertia_tensor_from_masses(masses, coords_active)
        evals, evecs = np.linalg.eigh(tensor)
        idx = np.argsort(evals)
        self.inertia_evals = evals[idx]
        self.inertia_evecs = evecs[:, idx]
        if np.linalg.det(self.inertia_evecs) < 0:
            self.inertia_evecs[:, 0] *= -1

    def _generate_candidate_axes(self) -> List[np.ndarray]:
        """Candidate rotation axes / plane normals.

        Atom-derived candidates come first (exact for symmetric geometries): atom
        directions, pair midpoints (C2 axes) and, for same-element pairs, the
        difference vectors r_i - r_j (mirror normals, perpendicular C2 axes) and
        the normal of the plane spanned by the pair and the origin.
        Principal inertia axes are appended last. Near-duplicates are removed with a
        tight threshold (|cos| > 0.9999, ~0.8 deg) so that an approximate axis can
        never shadow an exact one.
        """
        candidates = []
        active_indices = np.where(self.mask)[0]
        coords = self.coords
        syms = self.symbols

        for i in active_indices:
            if np.linalg.norm(coords[i]) > 0.1:
                candidates.append(normalize(coords[i]))

        if len(active_indices) < 50:
            for idx_i in range(len(active_indices)):
                i = active_indices[idx_i]
                for idx_j in range(idx_i + 1, len(active_indices)):
                    j = active_indices[idx_j]
                    mid = coords[i] + coords[j]
                    if np.linalg.norm(mid) > 0.1:
                        candidates.append(normalize(mid))
                    if syms[i] == syms[j]:
                        diff = coords[i] - coords[j]
                        if np.linalg.norm(diff) > 0.1:
                            candidates.append(normalize(diff))
                        # normal of the plane through the origin and the pair (C3/C4 axes of cubic groups)
                        cr = np.cross(coords[i], coords[j])
                        if np.linalg.norm(cr) > 0.1:
                            candidates.append(normalize(cr))

        for i in range(3):
            candidates.append(self.inertia_evecs[:, i])

        unique = []
        for c in candidates:
            if not any(abs(np.dot(c, u)) > 0.9999 for u in unique):
                unique.append(c)
        return unique

    def _check_op_mat(self, op_mat: np.ndarray) -> Tuple[bool, float]:
        coords_sub = self.coords[self.mask]
        syms_sub = self.symbols[self.mask]
        
        transformed = np.dot(coords_sub, op_mat.T)
        
        total_max_dist = 0.0
        
        for el in np.unique(syms_sub):
            mask_el = syms_sub == el
            orig_el = coords_sub[mask_el]
            trans_el = transformed[mask_el]
            
            dists = cdist(trans_el, orig_el)
            min_dists = np.min(dists, axis=1)
            current_max = np.max(min_dists)
            
            if current_max > total_max_dist:
                total_max_dist = current_max
                
            if current_max > self.tol:
                return False, current_max
                
        return True, total_max_dist

    def analyze(self):
        self.valid_ops = [np.eye(3)]
        self.max_op_rmsd = 0.0 
        
        valid, rmsd = self._check_op_mat(-np.eye(3))
        if valid:
            self.has_inversion = True
            self.valid_ops.append(-np.eye(3))
            self.max_op_rmsd = max(self.max_op_rmsd, rmsd)

        # Linear framework: every C_n about the molecular axis is valid; keep all orders
        # so that atoms mapped by C3/C5... (e.g. methyl H) are found equivalent.
        linear = self._is_linear()
        for axis in self.candidate_axes:
            for n in [8, 6, 5, 4, 3, 2]:
                mat = get_rotation_matrix(axis, 2*np.pi/n)
                valid, rmsd = self._check_op_mat(mat)
                if valid:
                    self.found_axes.append((axis, n))
                    self.max_op_rmsd = max(self.max_op_rmsd, rmsd)
                    for k in range(1, n):
                        mat_k = get_rotation_matrix(axis, k*2*np.pi/n)
                        self.valid_ops.append(mat_k)
                    if not linear:
                        break

        # Improper rotations S_n = sigma_h * C_n (S1 = sigma, S2 = i are handled separately;
        # S3/S5 imply C3h/C5h and are covered by C_n + sigma_h)
        for axis in self.candidate_axes:
            refl = get_reflection_matrix(axis)
            for n in [8, 6, 4]:
                mat = refl @ get_rotation_matrix(axis, 2*np.pi/n)
                valid, rmsd = self._check_op_mat(mat)
                if valid:
                    self.found_improper.append((axis, n))
                    self.max_op_rmsd = max(self.max_op_rmsd, rmsd)
                    op = np.eye(3)
                    for _ in range(1, n):
                        op = mat @ op
                        self.valid_ops.append(op.copy())
                    break

        for axis in self.candidate_axes:
            mat = get_reflection_matrix(axis)
            valid, rmsd = self._check_op_mat(mat)
            if valid:
                self.found_planes.append(axis)
                self.valid_ops.append(mat)
                self.max_op_rmsd = max(self.max_op_rmsd, rmsd)

    def _is_linear(self) -> bool:
        return self.inertia_evals[0] < 1e-4 and self.inertia_evals[1] >= 1e-4

    def _identify_group(self) -> str:
        if self.inertia_evals[1] < 1e-4:
            return "Kh"  # single atom
        if self._is_linear():
            return "D_inf_h" if self.has_inversion else "C_inf_v"

        if not self.found_axes:
            if self.found_planes: return "Cs"
            if self.has_inversion: return "Ci"
            return "C1"

        self.found_axes.sort(key=lambda x: x[1], reverse=True)
        n = self.found_axes[0][1]

        high_order = [x for x in self.found_axes if x[1] >= 3]
        unique_high = []
        for ax, order in high_order:
            if not any(abs(np.dot(ax, u[0])) > 0.9 for u in unique_high):
                unique_high.append((ax, order))

        if len(unique_high) > 1:
            orders = {order for _, order in unique_high}
            if 5 in orders: return "Ih" if self.has_inversion else "I"
            if 4 in orders: return "Oh" if self.has_inversion else "O"
            if self.has_inversion: return "Th"
            return "Td" if self.found_planes else "T"

        def coaxial_improper(ax, order):
            return any(m == order and abs(np.dot(ax, a)) > 0.9 for a, m in self.found_improper)

        # Principal axis: highest order; among equal orders (e.g. three C2 in D2d)
        # prefer the one carrying a coaxial S_2n.
        top = [ax for ax, order in self.found_axes if order == n]
        principal_axis = next((ax for ax in top if coaxial_improper(ax, 2 * n)), top[0])

        perp_c2s = 0
        for ax, order in self.found_axes:
            if order == 2 and abs(np.dot(ax, principal_axis)) < 0.1:
                perp_c2s += 1

        has_sigma_h = False
        has_sigma_v = False
        for normal in self.found_planes:
            dot = abs(np.dot(normal, principal_axis))
            if dot > 0.9: has_sigma_h = True
            if dot < 0.1: has_sigma_v = True

        if perp_c2s >= 1:
            if has_sigma_h: return f"D{n}h"
            if has_sigma_v: return f"D{n}d"
            return f"D{n}"
        else:
            if has_sigma_h: return f"C{n}h"
            if has_sigma_v: return f"C{n}v"
            if coaxial_improper(principal_axis, 2 * n): return f"S{2*n}"
            if self.has_inversion and n%2!=0: return f"S{2*n}"
            return f"C{n}"

    def get_equivalent_atoms(self, map_tolerance: float = 0.3) -> Dict[str, List[List[int]]]:
        if not self.valid_ops:
            self.analyze()

        dset = DisjointSet(self.n_atoms)

        for mat in self.valid_ops:
            if isinstance(mat, tuple):
                mat = mat[0]
            transformed = np.dot(self.coords, mat.T)

            for el in np.unique(self.symbols):
                indices = np.where(self.symbols == el)[0]
                if len(indices) == 0:
                    continue

                sub_orig = self.coords[indices]
                sub_trans = transformed[indices]

                dists = cdist(sub_orig, sub_trans)
                rows, cols = np.where(dists < map_tolerance)

                for r, c in zip(rows, cols):
                    dset.union(indices[r], indices[c])

        groups = dset.get_groups()
        output = defaultdict(list)
        for g in groups.values():
            el = self.symbols[g[0]]
            output[el].append(sorted(g))

        for k in output:
            output[k].sort(key=lambda x: x[0])

        return dict(output)

def detect_point_group(symbols: List[str], positions: np.ndarray, tolerance: float = 0.1, exclude_h: bool = False) -> Dict[str, Any]:
    """Schoenflies symbol and symmetry-equivalent atoms.

    exclude_h=False (default) gives the true molecular point group. exclude_h=True
    gives the group of the heavy-atom framework (centred on its own COM); hydrogens
    are then mapped onto equivalence classes with a 3x looser tolerance.
    """
    detector = PointGroupDetector(symbols, positions, tolerance, exclude_h=exclude_h)
    detector.analyze()
    label = detector._identify_group()
    map_tol = tolerance * 3 if exclude_h else tolerance
    equiv = detector.get_equivalent_atoms(map_tolerance=map_tol)
    
    return {
        'is_valid': True,
        'sch_symbol': label,
        'equivalencies': equiv,
        'error_rmsd': detector.max_op_rmsd
    }




def get_nmr_equivalencies(
    atoms: Any,
    coords: Any,
    buffer: float = 0.3,
) -> Dict[str, List[List[int]]]:
    """Groups of chemically equivalent atoms: heavy-atom point group filtered by
    Weisfeiler-Lehman colours on the bonded graph, plus rotors (>= 3 H on one atom).

    Returns {element: [[0-based indices], ...]}, every atom in exactly one group.
    """
    atoms = [str(a) for a in atoms]
    coords = np.array(coords, dtype=float)
    n_atoms = len(atoms)
    if n_atoms == 0:
        return {}

    # 2. Rigid global symmetry detection
    # exclude_h=True detects group on heavy framework, then maps all atoms including H
    pg_data = detect_point_group(atoms, coords, tolerance=0.1, exclude_h=True)
    global_equiv = pg_data.get('equivalencies', {})

    # 3. Local topology for rotators
    radii = {k: get_atomic_radius(k) for k in set(atoms)}
    adj = create_graph_matrix(atoms, coords, radii, buffer).astype(bool)

    # Topological colours (Weisfeiler-Lehman refinement on the full graph, H included).
    # Symmetry-equivalent atoms must share the same colour; this rejects spurious
    # equivalencies from a heavy-atom-only point group that the H atoms do not obey
    # (e.g. D_inf_h of a C-C framework carrying CH and CH3 ends).
    colours = list(atoms)
    for _ in range(min(n_atoms, 6)):
        colours = [colours[i] + '(' + ','.join(sorted(colours[n] for n in np.where(adj[i])[0])) + ')'
                   for i in range(n_atoms)]
        uniq = {c: str(k) for k, c in enumerate(sorted(set(colours)))}
        colours = [atoms[i] + uniq[colours[i]] for i in range(n_atoms)]

    # 4. Disjoint Set Union (DSU) merge
    dsu = DisjointSet(n_atoms)
    # Apply global symmetry equivalencies (only between topologically identical atoms)
    for el, groups in global_equiv.items():
        for group in groups:
            if len(group) > 1:
                by_colour: Dict[str, List[int]] = defaultdict(list)
                for member in group:
                    by_colour[colours[int(member)]].append(int(member))
                for sub in by_colour.values():
                    for member in sub[1:]:
                        dsu.union(sub[0], member)

    # Apply rotor equivalencies (e.g. -CH3, -NH3+)
    for i, atom in enumerate(atoms):
        if atom != 'H':
            neighbors = np.where(adj[i])[0]
            h_neighbors = [int(n) for n in neighbors if atoms[n] == 'H']
            if len(h_neighbors) >= 3:
                root = h_neighbors[0]
                for h in h_neighbors[1:]:
                    dsu.union(root, h)

    # Build final groupings grouped by element
    output: Dict[str, List[List[int]]] = defaultdict(list)
    for g in dsu.get_groups().values():
        el = atoms[g[0]]
        output[el].append(sorted(g))

    for k in output:
        output[k].sort(key=lambda x: x[0])

    return dict(output)
