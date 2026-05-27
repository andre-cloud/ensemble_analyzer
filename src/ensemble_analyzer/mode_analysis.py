import numpy as np
from typing import Optional


THRESHOLD_IMAG = 20.0


class NormalModeAnalyzer:
    def __init__(
        self,
        normal_modes: np.ndarray,
        geom: np.ndarray,
        atoms: tuple[str, ...],
    ):
        self.normal_modes = np.asarray(normal_modes, dtype=float)
        self.geom = np.asarray(geom, dtype=float)
        self.atoms = atoms
        self.n_atoms = len(atoms)
        self.n_modes = self.normal_modes.shape[0]

    def localize_mode(self, mode: int) -> np.ndarray:
        displ = self.normal_modes[mode]
        sq = np.sum(displ ** 2, axis=1)
        total = np.sum(sq)
        if total < 1e-14:
            return np.zeros(self.n_atoms)
        return (sq / total) * 100.0

    def localize_mode_fragment(
        self,
        mode: int,
        fragments: dict[str, list[int]],
    ) -> dict[str, float]:
        atomic = self.localize_mode(mode)
        return {
            name: float(np.sum(atomic[indices]))
            for name, indices in fragments.items()
        }

    def displace_geometry(
        self, mode: int, scale: float = 0.3
    ) -> np.ndarray:
        return self.geom + scale * self.normal_modes[mode]

    @staticmethod
    def classify_negative_freqs(
        freqs: np.ndarray,
        threshold: float = THRESHOLD_IMAG,
    ) -> tuple[list[int], list[int]]:
        neg = np.where(freqs < 0)[0]
        significant = [int(i) for i in neg if abs(freqs[i]) > threshold]
        noise = [int(i) for i in neg if abs(freqs[i]) <= threshold]
        return significant, noise

    def imag_mode_summary(
        self,
        mode: int,
        fragments: Optional[dict[str, list[int]]] = None,
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
