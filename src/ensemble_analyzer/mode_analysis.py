import numpy as np
from ase import Atoms as ASE_Atoms


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