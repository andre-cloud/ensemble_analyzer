import numpy as np
from ase import Atoms as ASE_Atoms
from ase import units


_ATOMIC_MASSES: dict[str, float] = {}

# Converts sqrt(eigenvalue) of a mass-weighted Hessian in eV/(Å² amu) to eV,
# the same factor ase.vibrations uses.
_HESSIAN_TO_EV = units._hbar * 1e10 / np.sqrt(units._e * units._amu)


def rigid_body_basis(
    masses: np.ndarray, positions: np.ndarray, tol: float = 1e-6
) -> np.ndarray:
    """Orthonormal basis (3N, k) of the mass-weighted rigid-body motions.

    Spans the three translations and the three rotations about the centre of
    mass. k is 6 for a non-linear molecule, 5 for a linear one and 3 for a
    single atom: null rotations are dropped through the SVD, so no linearity
    flag is needed.
    """
    m = np.asarray(masses, dtype=float)
    x = np.asarray(positions, dtype=float)
    n = len(m)
    sqm = np.sqrt(m)
    r = x - (m[:, None] * x).sum(axis=0) / m.sum()

    vecs = []
    for a in range(3):
        t = np.zeros((n, 3))
        t[:, a] = sqm
        vecs.append(t.ravel())
    for a in range(3):
        e = np.zeros(3)
        e[a] = 1.0
        vecs.append((sqm[:, None] * np.cross(e, r)).ravel())

    U, s, _ = np.linalg.svd(np.array(vecs).T, full_matrices=False)
    return U[:, s > tol * s.max()]


def project_hessian(
    hessian: np.ndarray, masses: np.ndarray, positions: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Vibrational frequencies and modes with translations/rotations removed.

    The mass-weighted Hessian is diagonalised only inside the internal
    subspace orthogonal to the rigid-body motions, so exactly 3N-6 (3N-5 for
    linear molecules) modes are returned and no near-zero spurious mode has to
    be recognised and discarded afterwards.

    Args:
        hessian: (3N, 3N) Cartesian Hessian in eV/Å² (ASE units).
        masses: (N,) atomic masses in amu.
        positions: (N, 3) Cartesian coordinates in Å.

    Returns:
        freqs: (n_modes,) frequencies in cm⁻¹, imaginary ones as negative values,
            sorted by ascending eigenvalue.
        modes: (n_modes, N, 3) mass-weighted, orthonormal normal modes (same
            convention as ``Vibrations.get_mode(i) * sqrt(m)``).
    """
    m = np.asarray(masses, dtype=float)
    n = len(m)
    H = np.asarray(hessian, dtype=float)
    H = 0.5 * (H + H.T)
    inv_sqm = np.repeat(1.0 / np.sqrt(m), 3)
    Hmw = H * inv_sqm[:, None] * inv_sqm[None, :]

    D = rigid_body_basis(m, positions)
    k = D.shape[1]
    if 3 * n - k <= 0:
        return np.zeros(0), np.zeros((0, n, 3))

    # Orthonormal complement of D: the internal (vibrational) subspace
    Q, _ = np.linalg.qr(np.hstack([D, np.eye(3 * n)]))
    B = Q[:, k:3 * n]

    w2, L = np.linalg.eigh(B.T @ Hmw @ B)
    energies = _HESSIAN_TO_EV * np.sqrt(np.abs(w2)) * np.sign(w2)
    freqs = energies / units.invcm
    modes = (B @ L).T.reshape(-1, n, 3)
    return freqs, modes


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
        real_displacements = self.normal_modes[mode]
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

    def delta_positions(
        self, mode: int, scale: float = 0.3
    ) -> dict:
        mode_vec = self.normal_modes[mode]
        delta_mag = 2 * scale * np.linalg.norm(mode_vec, axis=1)
        total = np.sum(delta_mag)
        if total == 0:
            return {
                "delta_mag": np.zeros(self.n_atoms),
                "percent": np.zeros(self.n_atoms),
                "total_delta": 0.0,
            }
        percent = (delta_mag / total) * 100.0
        return {
            "delta_mag": delta_mag,
            "percent": percent,
            "total_delta": float(total),
        }

    def displace_geometry(
        self, mode: int, scale: float = 0.3
    ) -> np.ndarray:
        real_displacements = self.normal_modes[mode]
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