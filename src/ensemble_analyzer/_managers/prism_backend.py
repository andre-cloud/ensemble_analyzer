"""
Similarity-pruning backend based on PRISM Pruner (selected with prune_method='prism').

PRISM (Pruning Interface for Similar Molecules, N. Tampellini,
https://github.com/ntampellini/prism_pruner) removes duplicate conformers by
up to three metrics:

* ``moi``     - relative deviation of the principal moments of inertia
                (mirror-image conformers are treated as duplicates);
* ``rmsd``    - Kabsch RMSD + maximum atomic deviation (handedness-aware);
* ``rotcorr`` - RMSD corrected for degenerate rotamers (slowest).

Behaviour of ``prism_pruner`` 0.1.2 that this module compensates for:

1. When energies are supplied, the sequential code path (always used by the
   MOI filter, and by the RMSD filter on large ensembles) re-sorts structures
   by energy internally, which mis-aligns precomputed descriptors with the
   structures and can discard a *unique* conformer. Input is therefore always
   passed pre-sorted by strictly increasing energy, which makes PRISM's
   internal sort the identity permutation.
2. Which member of a duplicate pair survives is not energy-aware (RMSD keeps
   the first, MOI keeps the last). After every stage the representative is
   re-assigned to the lowest-energy member of the pair.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np

PRISM_METHODS = ("prism", "prism_rotcorr")
ALL_PRUNE_METHODS = ("rotconst",) + PRISM_METHODS

# Strictly increasing offset added to sorted energies [kcal/mol]; far below any
# meaningful threshold, only used to break ties so that sorting is the identity.
_TIE_BREAK = 1e-9


@dataclass
class PrismRecord:
    """One conformer removed by a PRISM stage."""
    removed: int          # index in the (energy-sorted) input
    reference: int        # index of the retained duplicate partner
    stage: str            # 'moi' | 'rmsd' | 'rotcorr'
    metric: float         # MOI max relative deviation, or Kabsch RMSD [Å]


def import_prism():
    """Import ``prism_pruner.pruner`` or raise an informative ImportError."""
    try:
        from prism_pruner import pruner  # noqa: F401
        return pruner
    except ImportError as exc:  # pragma: no cover - message path tested via mock
        raise ImportError(
            "prune_method='prism'/'prism_rotcorr' requires 'prism_pruner', a core "
            "dependency of ensemble-analyzer that is missing from this environment. "
            "Reinstall the package (`pip install -U ensemble-analyzer`) or run "
            "`pip install prism_pruner`."
        ) from exc


# ---------------------------------------------------------------------------
# Lightweight descriptors used only to identify the duplicate partner
# ---------------------------------------------------------------------------

def _masses(atoms: Sequence[str]) -> np.ndarray:
    from ase.data import atomic_masses, atomic_numbers
    return np.array([atomic_masses[atomic_numbers[a]] for a in atoms])


def principal_moments(coords: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Sorted principal moments of inertia of one structure."""
    com = (coords * masses[:, None]).sum(0) / masses.sum()
    r = coords - com
    inertia = np.einsum("i,ij,ik->jk", masses, r, r)
    inertia = np.eye(3) * np.trace(inertia) - inertia
    return np.linalg.eigvalsh(inertia)


def kabsch_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Proper-rotation (no reflection) Kabsch RMSD between two structures."""
    p = a - a.mean(0)
    q = b - b.mean(0)
    u, s, vt = np.linalg.svd(p.T @ q)
    if np.linalg.det(u @ vt) < 0:
        s[-1] = -s[-1]
    msd = (np.sum(p * p) + np.sum(q * q) - 2.0 * s.sum()) / len(a)
    return float(np.sqrt(max(msd, 0.0)))


def _partner(
    r: int, candidates: np.ndarray, coords: np.ndarray, energies: np.ndarray,
    max_dE: float, stage: str, moments: Optional[np.ndarray], sel: np.ndarray,
) -> Tuple[Optional[int], float]:
    """Most similar retained structure to ``r`` within the energy window."""
    cands = candidates[np.abs(energies[candidates] - energies[r]) < max_dE]
    if cands.size == 0:
        return None, np.nan
    if stage == "moi":
        dev = np.max(np.abs(moments[cands] - moments[r]) / moments[r], axis=1)
    else:
        dev = np.array([kabsch_rmsd(coords[r][sel], coords[c][sel]) for c in cands])
    best = int(np.argmin(dev))
    return int(cands[best]), float(dev[best])


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_prism_pruning(
    coords: np.ndarray,
    atoms: Sequence[str],
    energies: np.ndarray,
    *,
    max_dE: float,
    method: str = "prism",
    use_moi: bool = False,
    moi_dev: float = 1e-2,
    max_rmsd: float = 0.25,
    max_dev: Optional[float] = None,
    heavy_atoms_only: bool = True,
    timeout_s: int = 60,
    debugfunction: Optional[Callable[[str], None]] = None,
) -> Tuple[np.ndarray, List[PrismRecord]]:
    """
    Run the PRISM stages on an ensemble and return which structures to keep.

    Args:
        coords: (n, N, 3) coordinates [Å], **sorted by ascending energy**.
        atoms: N atomic symbols.
        energies: (n,) relative energies [kcal/mol], ascending.
        max_dE: only pairs with |ΔE| < max_dE are compared [kcal/mol].
        method: 'prism' (rmsd) or 'prism_rotcorr' (rmsd + rotamer-corrected rmsd).
        use_moi: also run the moment-of-inertia stage first (opt-in; at 1 % it
            merges distinct conformers of large molecules, see docs).
        moi_dev: relative MOI deviation threshold (0.01 = 1 %).
        max_rmsd: RMSD threshold [Å].
        max_dev: max atomic deviation threshold [Å] (default 2*max_rmsd).
        heavy_atoms_only: exclude H from RMSD stages.
        timeout_s: per-stage PRISM timeout [s].

    Returns:
        keep: (n,) boolean mask over the input order.
        records: list of PrismRecord, one per removed structure.
    """
    if method not in PRISM_METHODS:
        raise ValueError(f"Unknown PRISM method {method!r}; expected one of {PRISM_METHODS}")

    pruner = import_prism()
    coords = np.asarray(coords, dtype=float)
    atoms = np.asarray(atoms)
    energies = np.asarray(energies, dtype=float)
    n = len(coords)
    keep = np.ones(n, dtype=bool)
    records: List[PrismRecord] = []
    if n < 2:
        return keep, records

    if np.any(np.diff(energies) < 0):
        raise ValueError("run_prism_pruning expects structures sorted by ascending energy")

    # Strictly increasing energies -> PRISM's internal argsort is the identity.
    e_in = energies - energies.min() + np.arange(n) * _TIE_BREAK
    max_dev = max_dev if max_dev is not None else 2.0 * max_rmsd
    sel = (atoms != "H") if heavy_atoms_only else np.ones(len(atoms), dtype=bool)
    moments = np.array([principal_moments(c, _masses(atoms)) for c in coords]) if use_moi else None

    stages = (["moi"] if use_moi else []) + ["rmsd"] + (["rotcorr"] if method == "prism_rotcorr" else [])

    for stage in stages:
        idx = np.flatnonzero(keep)
        if idx.size < 2:
            break
        S, E = coords[idx], e_in[idx]

        if stage == "moi":
            _, mask = pruner.prune_by_moment_of_inertia(
                S, atoms, max_deviation=moi_dev, energies=E, max_dE=max_dE,
                timeout_s=timeout_s, debugfunction=debugfunction)
        elif stage == "rmsd":
            _, mask = pruner.prune_by_rmsd(
                S, atoms, max_rmsd=max_rmsd, max_dev=max_dev, energies=E, max_dE=max_dE,
                timeout_s=timeout_s, heavy_atoms_only=heavy_atoms_only,
                debugfunction=debugfunction)
        else:
            from prism_pruner.graph_manipulations import graphize
            graph = graphize(atoms, S[0])
            _, mask = pruner.prune_by_rmsd_rot_corr(
                S, atoms, graph, max_rmsd=max_rmsd, max_dev=max_dev, energies=E,
                max_dE=max_dE, timeout_s=timeout_s, heavy_atoms_only=heavy_atoms_only,
                debugfunction=debugfunction)

        mask = np.asarray(mask, dtype=bool)
        if mask.shape != (idx.size,):
            raise RuntimeError(f"PRISM stage '{stage}' returned a mask of shape {mask.shape}")

        removed = idx[~mask]
        retained = idx[mask]
        keep[removed] = False

        # Identify the partner of each removed structure and make sure the
        # lowest-energy member of the pair is the one that survives.
        for r in removed:
            p, metric = _partner(r, retained, coords, energies, max_dE, stage, moments, sel)
            if p is None:
                # No partner found (should not happen): keep PRISM's decision.
                records.append(PrismRecord(int(r), -1, stage, np.nan))
                continue
            if energies[r] < energies[p]:
                keep[r], keep[p] = True, False
                retained = np.where(retained == p, r, retained)
                # re-point earlier records that referenced p
                for rec in records:
                    if rec.reference == p:
                        rec.reference = int(r)
                records.append(PrismRecord(int(p), int(r), stage, metric))
            else:
                records.append(PrismRecord(int(r), int(p), stage, metric))

    return keep, records
