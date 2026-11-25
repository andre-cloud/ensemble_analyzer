"""1D spectrum of one observed nucleus, with homonuclear and heteronuclear partners.

Homonuclear couplings go into the exact spin Hamiltonian (:mod:`.spin_sim`). Heteronuclear
partners are weakly coupled (only the Iz·Sz term survives), so for every state of the partner
groups (their total M, with degeneracy and natural abundance) the observed resonances move by
Σ J·M and the homonuclear problem is solved again: exact to second order in the homonuclear part.

- A partner whose isotope has natural abundance < ``min_abundance`` is ignored (13C satellites
  in a 1H spectrum); one above it contributes with its abundance (e.g. 119Sn satellites).
- Homonuclear couplings are ignored when the observed isotope itself is below ``min_abundance``
  (no 13C–13C splitting at natural abundance).
- Partners in ``decouple`` are ignored.
"""

from collections import defaultdict
from dataclasses import dataclass, field
from itertools import product
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.sparse.csgraph import connected_components

from ensemble_analyzer._nmr.nuclei import (
    abundance,
    get_nucleus_spin,
    larmor_mhz,
    normalize_nucleus_name,
)
from ensemble_analyzer._nmr.spin_sim import SpinTransition, simulate_quantum_transitions

HETERO_MIN_ABUNDANCE = 0.05


@dataclass
class SpinGroup:
    """A set of equivalent nuclei: isotope, shift (ppm), number of nuclei, atom indices."""

    nucleus: str
    delta_ppm: float
    count: int = 1
    members: Tuple[int, ...] = field(default_factory=tuple)


def _partner_distribution(spin: float, n: int, p: float, prune: float = 1e-3) -> Dict[float, float]:
    """Probability of each total M of n partner nuclei, each magnetic with probability p.

    States below ``prune`` are dropped and the rest renormalised (no 1H-free isotopologues).
    """
    ms = [-spin + i for i in range(int(round(2 * spin + 1)))]
    single: Dict[float, float] = defaultdict(float)
    single[0.0] += 1.0 - p
    for m in ms:
        single[m] += p / len(ms)
    dist: Dict[float, float] = {0.0: 1.0}
    for _ in range(n):
        new: Dict[float, float] = defaultdict(float)
        for m_tot, w in dist.items():
            for m, ws in single.items():
                if ws > 0:
                    new[round(m_tot + m, 6)] += w * ws
        dist = dict(new)
    kept = {m: w for m, w in dist.items() if w > prune}
    total = sum(kept.values())
    return {m: w / total for m, w in kept.items()}


def _components(n: int, J: np.ndarray, cutoff: float) -> List[List[int]]:
    """Groups connected by |J| >= cutoff, ordered by their lowest index."""
    _, labels = connected_components(np.abs(J[:n, :n]) >= cutoff, directed=False)
    return [np.flatnonzero(labels == c).tolist() for c in range(labels.max(initial=-1) + 1)]


def simulate_nucleus(
    observed: str,
    groups: Sequence[SpinGroup],
    j_hz: np.ndarray,
    frequency_1h: float,
    *,
    decouple: Sequence[str] = (),
    j_cutoff: float = 0.3,
    min_abundance: float = HETERO_MIN_ABUNDANCE,
    max_configurations: int = 4096,
    min_intensity: float = 1e-4,
    max_cluster_size: int = 8,
    max_hilbert_dim: int = 256,
) -> List[SpinTransition]:
    """Transitions of the ``observed`` isotope.

    Args:
        observed: isotope whose spectrum is simulated (e.g. ``"1H"``).
        groups: all groups of the molecule (any isotope); only coupled partners matter.
        j_hz: group × group coupling matrix in Hz.
        frequency_1h: spectrometer frequency for 1H in MHz.
        decouple: partner isotopes to ignore.

    Returns:
        Transitions with ``freq_hz`` measured from 0 ppm at the observed Larmor frequency;
        the intensities of a group sum to its ``count``.
    """
    obs = normalize_nucleus_name(observed)
    nu0 = larmor_mhz(obs, frequency_1h)
    J = np.asarray(j_hz, dtype=float)
    iso = [normalize_nucleus_name(g.nucleus) for g in groups]
    obs_idx = [i for i, n in enumerate(iso) if n == obs]
    if not obs_idx:
        return []

    dec = {normalize_nucleus_name(d) for d in decouple}
    partners = [
        i for i, n in enumerate(iso)
        if n != obs and n not in dec and abundance(n) >= min_abundance
    ]
    J_obs = J[np.ix_(obs_idx, obs_idx)].copy()
    if abundance(obs) < min_abundance:
        J_obs[:] = 0.0

    transitions: List[SpinTransition] = []
    for cid, comp in enumerate(_components(len(obs_idx), J_obs, j_cutoff)):
        sub = [obs_idx[c] for c in comp]
        coupled = [k for k in partners if any(abs(J[g, k]) >= j_cutoff for g in sub)]
        coupled.sort(key=lambda k: -max(abs(J[g, k]) for g in sub))
        dists = {k: _partner_distribution(get_nucleus_spin(iso[k]), groups[k].count, abundance(iso[k]))
                 for k in coupled}

        exact, first_order, n_conf = [], [], 1
        for k in coupled:
            if n_conf * len(dists[k]) <= max_configurations:
                exact.append(k)
                n_conf *= len(dists[k])
            else:
                first_order.append(k)

        base = np.array([groups[g].delta_ppm for g in sub], dtype=float)
        counts = [groups[g].count for g in sub]
        J_sub = J_obs[np.ix_(comp, comp)]
        lines: Dict[float, float] = defaultdict(float)
        for combo in product(*[list(dists[k].items()) for k in exact]):
            weight = float(np.prod([w for _, w in combo])) if combo else 1.0
            shift_hz = np.array([sum(J[g, k] * m for k, (m, _) in zip(exact, combo)) for g in sub])
            for t in simulate_quantum_transitions(
                base + shift_hz / nu0, J_sub, group_counts=counts, freq_mhz=nu0,
                j_cutoff=j_cutoff, min_intensity=min_intensity, nucleus=obs,
                max_cluster_size=max_cluster_size, max_hilbert_dim=max_hilbert_dim,
            ):
                lines[round(t.freq_hz, 6)] += t.intensity * weight

        # Partners beyond the configuration budget: first order, line → closest group
        base_hz = base * nu0
        for k in first_order:
            new: Dict[float, float] = defaultdict(float)
            for f, w in lines.items():
                g = sub[int(np.argmin(np.abs(base_hz - f)))]
                for m, pm in dists[k].items():
                    new[round(f + J[g, k] * m, 6)] += w * pm
            lines = new

        transitions.extend(
            SpinTransition(freq_hz=f, intensity=w, cluster_id=cid, delta_ppm=f / nu0)
            for f, w in lines.items() if w >= min_intensity
        )
    return transitions


def stick_spectrum(transitions: Sequence[SpinTransition], observed: str,
                   frequency_1h: float) -> Tuple[np.ndarray, np.ndarray]:
    """Line positions (ppm) and intensities, sorted by position."""
    nu0 = larmor_mhz(observed, frequency_1h)
    if not transitions:
        return np.zeros(0), np.zeros(0)
    ppm = np.array([t.freq_hz / nu0 for t in transitions])
    inten = np.array([t.intensity for t in transitions])
    order = np.argsort(ppm)
    return ppm[order], inten[order]
