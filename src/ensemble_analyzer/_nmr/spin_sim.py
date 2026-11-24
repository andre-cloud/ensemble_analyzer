"""Spin simulation of 1D NMR spectra (full isotropic spin Hamiltonian per cluster).

Copied from FROGChem (paclab-lib ``frogchem/nmr_sim.py``: the 1D core only, no 2D/FID/plots);
keep in sync. Changes: nuclear data from :mod:`ensemble_analyzer._nmr.nuclei`.
"""

import math
from collections import defaultdict
from dataclasses import dataclass
from itertools import product
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from ensemble_analyzer._nmr.nuclei import get_nucleus_spin


@dataclass
class SpinTransition:
    """Represents an allowed NMR dipole transition."""
    freq_hz: float
    intensity: float
    cluster_id: int = 0
    delta_ppm: float = 0.0


def solve_spin_cluster(
    nu: np.ndarray,
    J: np.ndarray,
    spins: Optional[Union[float, Sequence[float]]] = None,
    min_intensity: float = 1e-5,
) -> List[Tuple[float, float]]:
    """Solve the isotropic spin Hamiltonian for a cluster of n coupled spins with arbitrary
    spin quantum numbers S_k (e.g. S=0.5 for 1H/13C, S=1 for 2H/14N, S=1.5 for 33S/11B, S=2.5 for 17O)
    using block diagonalization by total magnetic quantum number M = sum(m_k).

    Parameters:
        nu: (n,) array of chemical shift offsets / resonance frequencies in Hz.
        J: (n, n) symmetric scalar coupling matrix in Hz.
        spins: Single spin quantum number (float, e.g. 2.5) or sequence of spin quantum numbers
               for each nucleus in the cluster. If None, defaults to 0.5 (spin-1/2).
        min_intensity: Minimum transition probability to keep.

    Returns:
        List of (frequency_hz, transition_probability) tuples, one per Delta M = +1
        transition (frequencies keep their sign: negative for delta < 0 ppm).
    """
    n = len(nu)
    if n == 0:
        return []

    if spins is None:
        spin_list = [0.5] * n
    elif isinstance(spins, (int, float)):
        spin_list = [float(spins)] * n
    else:
        spin_list = [float(s) for s in spins]

    if len(spin_list) != n:
        raise ValueError(f"Length of spins ({len(spin_list)}) must match number of frequencies ({n})")

    if n == 1:
        return [(float(nu[0]), 1.0)]

    # Product Zeeman basis states: for each spin k, m_k in [-S_k, -S_k+1, ..., +S_k]
    m_ranges = []
    for s in spin_list:
        d = int(round(2.0 * s + 1.0))
        m_ranges.append([-s + i for i in range(d)])

    all_states = [np.array(st, dtype=float) for st in product(*m_ranges)]

    # 1. Group Zeeman product basis states by total magnetic quantum number M = sum(m_k)
    blocks: Dict[float, List[int]] = defaultdict(list)
    for idx, st in enumerate(all_states):
        M = round(float(np.sum(st)), 4)
        blocks[M].append(idx)

    # 2. Construct and diagonalize each M block independently
    eigvals: Dict[float, np.ndarray] = {}
    eigvecs: Dict[float, np.ndarray] = {}

    for M, b_list in blocks.items():
        dim = len(b_list)
        H_block = np.zeros((dim, dim), dtype=float)

        for i, bi in enumerate(b_list):
            st_i = all_states[bi]
            # Diagonal element: sum nu_k m_k + sum_{k<l} J_kl m_k m_l
            diag = float(np.sum(nu * st_i))
            for k in range(n):
                for l in range(k + 1, n):
                    diag += J[k, l] * st_i[k] * st_i[l]
            H_block[i, i] = diag

            # Off-diagonal element: 1/2 J_kl (I+_k I-_l + I-_k I+_l)
            for j in range(i + 1, dim):
                bj = b_list[j]
                st_j = all_states[bj]
                diff = st_j - st_i
                non_zero = np.where(np.abs(diff) > 1e-4)[0]
                if len(non_zero) == 2:
                    k, l = non_zero[0], non_zero[1]
                    dk, dl = diff[k], diff[l]
                    if np.isclose(dk, 1.0) and np.isclose(dl, -1.0):
                        # state j has m_k = m_i + 1 and m_l = m_i - 1
                        val_k = math.sqrt(spin_list[k] * (spin_list[k] + 1) - st_i[k] * (st_i[k] + 1))
                        val_l = math.sqrt(spin_list[l] * (spin_list[l] + 1) - st_i[l] * (st_i[l] - 1))
                        val = 0.5 * J[k, l] * val_k * val_l
                        H_block[i, j] = val
                        H_block[j, i] = val
                    elif np.isclose(dk, -1.0) and np.isclose(dl, 1.0):
                        # state j has m_k = m_i - 1 and m_l = m_i + 1
                        val_k = math.sqrt(spin_list[k] * (spin_list[k] + 1) - st_i[k] * (st_i[k] - 1))
                        val_l = math.sqrt(spin_list[l] * (spin_list[l] + 1) - st_i[l] * (st_i[l] + 1))
                        val = 0.5 * J[k, l] * val_k * val_l
                        H_block[i, j] = val
                        H_block[j, i] = val

        w, v = np.linalg.eigh(H_block)
        eigvals[M] = w
        eigvecs[M] = v

    # 3. Calculate dipole transitions between adjacent blocks: Delta M = +1
    transitions: List[Tuple[float, float]] = []
    sorted_M = sorted(blocks.keys())

    for idx_M, M in enumerate(sorted_M[:-1]):
        M_next = sorted_M[idx_M + 1]
        if not np.isclose(M_next - M, 1.0):
            continue

        b_M = blocks[M]
        b_M1 = blocks[M_next]
        dim_M = len(b_M)
        dim_M1 = len(b_M1)

        # Build raising operator I+ matrix: connects |st_j> in M to |st_i> in M+1 if exactly one spin k has m_i = m_j + 1
        I_plus = np.zeros((dim_M1, dim_M), dtype=float)
        for j, bj in enumerate(b_M):
            st_j = all_states[bj]
            for i, bi in enumerate(b_M1):
                st_i = all_states[bi]
                diff = st_i - st_j
                non_zero = np.where(np.abs(diff) > 1e-4)[0]
                if len(non_zero) == 1:
                    k = non_zero[0]
                    if np.isclose(diff[k], 1.0):
                        I_plus[i, j] = math.sqrt(spin_list[k] * (spin_list[k] + 1) - st_j[k] * (st_j[k] + 1))

        # Transform to eigenbasis: <psi_b^(M+1) | I+ | psi_a^(M)>
        M_trans = eigvecs[M_next].T @ (I_plus @ eigvecs[M])
        w_M = eigvals[M]
        w_M1 = eigvals[M_next]

        for i in range(dim_M1):
            for j in range(dim_M):
                freq = float(w_M1[i] - w_M[j])
                prob = float(M_trans[i, j] ** 2)
                # Only Delta M = +1 is enumerated, so no line is double counted: keep
                # negative frequencies too (nuclei with delta < 0 ppm).
                if prob >= min_intensity:
                    transitions.append((freq, prob))

    return transitions


def partition_spin_network(
    n_spins: int,
    J: np.ndarray,
    cutoff: float = 0.3,
    max_cluster_size: int = 8,
    spins: Optional[Sequence[float]] = None,
    max_dim: int = 256,
    counts: Optional[Sequence[int]] = None,
    fallback: str = 'singles',
) -> List[List[int]]:
    """Partition spin network into connected components (clusters) based on J threshold.
    If a cluster exceeds max_cluster_size or total Hilbert space dimension exceeds max_dim,
    its internal threshold is iteratively increased to prune weaker couplings, preventing
    combinatorial Hilbert space explosion (crucial for high-spin nuclei like 17O where (2S+1)=6).

    Parameters:
        counts: Optional number of equivalent spins represented by each node (group). The
            cluster size is then sum(counts) and its dimension prod((2S+1)**count).
        fallback: What to do with a cluster that cannot be subdivided into valid clusters:
            'singles' (default) splits it into isolated nodes (couplings lost);
            'none' returns it unchanged, so the caller can treat it differently
            (e.g. first order, as simulate_quantum_transitions does).
    """
    spin_list = [0.5] * n_spins if spins is None else [float(s) for s in spins]
    count_list = [1] * n_spins if counts is None else [int(c) for c in counts]

    def cluster_dim(c_indices: List[int]) -> int:
        d = 1
        for idx in c_indices:
            d *= int(round(2.0 * spin_list[idx] + 1.0)) ** count_list[idx]
        return d

    def is_valid_cluster(c_indices: List[int]) -> bool:
        size = sum(count_list[idx] for idx in c_indices)
        return size <= max_cluster_size and cluster_dim(c_indices) <= max_dim

    adj = np.abs(J) >= cutoff
    visited = [False] * n_spins
    raw_clusters: List[List[int]] = []

    for i in range(n_spins):
        if not visited[i]:
            cluster: List[int] = []
            queue = [i]
            visited[i] = True
            while queue:
                curr = queue.pop(0)
                cluster.append(curr)
                for neighbor in range(n_spins):
                    if adj[curr, neighbor] and not visited[neighbor]:
                        visited[neighbor] = True
                        queue.append(neighbor)
            raw_clusters.append(sorted(cluster))

    final_clusters: List[List[int]] = []
    for c in raw_clusters:
        if is_valid_cluster(c):
            final_clusters.append(c)
        else:
            # Subdivide large cluster by progressively raising threshold
            c_J = np.zeros((len(c), len(c)), dtype=float)
            for i, a in enumerate(c):
                for j, b in enumerate(c):
                    c_J[i, j] = J[a, b]
            c_spins = [spin_list[idx] for idx in c]
            c_counts = [count_list[idx] for idx in c]

            upper_vals = sorted(list(set(np.abs(c_J[np.triu_indices(len(c), k=1)]))))
            upper_vals = [v for v in upper_vals if v >= cutoff]

            sub_assigned = False
            for thresh in upper_vals:
                sub_c = partition_spin_network(
                    len(c),
                    c_J,
                    cutoff=thresh + 1e-4,
                    max_cluster_size=max_cluster_size,
                    spins=c_spins,
                    max_dim=max_dim,
                    counts=c_counts,
                    fallback='none',
                )
                if all(is_valid_cluster([c[idx] for idx in sc]) for sc in sub_c):
                    for sc in sub_c:
                        final_clusters.append(sorted([c[idx] for idx in sc]))
                    sub_assigned = True
                    break
            if not sub_assigned:
                if fallback == 'none':
                    final_clusters.append(c)
                else:
                    # Fallback: single spins
                    for idx in c:
                        final_clusters.append([idx])

    return final_clusters


def _m_distribution(spin: float, n: int) -> Dict[float, float]:
    """Degeneracy of the total M of n uncoupled spins I (n spin-1/2 -> binomial / Pascal row)."""
    single = [-spin + i for i in range(int(round(2.0 * spin + 1.0)))]
    dist: Dict[float, float] = {0.0: 1.0}
    for _ in range(n):
        new: Dict[float, float] = defaultdict(float)
        for m_tot, w in dist.items():
            for m in single:
                new[round(m_tot + m, 6)] += w
        dist = dict(new)
    return dist


def _first_order_multiplets(
    cluster: Sequence[int],
    nu_hz: np.ndarray,
    J_arr: np.ndarray,
    counts: Sequence[int],
    spins: Sequence[float],
    j_cutoff: float,
    min_intensity: float,
) -> List[Tuple[float, float]]:
    """First-order (weak coupling, |J| << |dnu|) multiplets of every group in a cluster.

    Group k (n_k nuclei) gives lines at nu_k + sum_j J_kj * M_j, where M_j runs over the
    total magnetic quantum number of the n_j spins of each coupled group j, weighted by its
    degeneracy (Pascal triangle for spin 1/2). The lines of group k sum to n_k.
    """
    out: List[Tuple[float, float]] = []
    for k in cluster:
        lines: Dict[float, float] = {float(nu_hz[k]): 1.0}
        for j in cluster:
            if j == k or abs(J_arr[k, j]) < j_cutoff:
                continue
            dist = _m_distribution(spins[j], counts[j])
            new: Dict[float, float] = defaultdict(float)
            for f, w in lines.items():
                for m_tot, deg in dist.items():
                    new[round(f + J_arr[k, j] * m_tot, 6)] += w * deg
            lines = dict(new)
        total = sum(lines.values())
        for f, w in lines.items():
            inten = counts[k] * w / total
            if inten >= min_intensity:
                out.append((f, inten))
    return out


def _add_external_couplings(
    lines: List[Tuple[float, float]],
    cluster: Sequence[int],
    nu_hz: np.ndarray,
    J_arr: np.ndarray,
    counts: Sequence[int],
    spins: Sequence[float],
    j_cutoff: float,
    min_intensity: float,
) -> List[Tuple[float, float]]:
    """Apply to first order the couplings (|J| >= j_cutoff) between a cluster and the groups
    outside it (couplings dropped by the partitioning). Each line is attributed to the cluster
    group with the closest resonance frequency and split by the total-M distribution of every
    external group coupled to it; intensities are conserved."""
    outside = [j for j in range(len(nu_hz)) if j not in set(cluster)]
    if not outside or not lines:
        return lines
    pattern: Dict[int, Dict[float, float]] = {}
    for k in cluster:
        offs: Dict[float, float] = {0.0: 1.0}
        for j in outside:
            if abs(J_arr[k, j]) < j_cutoff:
                continue
            dist = _m_distribution(spins[j], counts[j])
            tot = sum(dist.values())
            new: Dict[float, float] = defaultdict(float)
            for o, w in offs.items():
                for m_tot, deg in dist.items():
                    new[round(o + J_arr[k, j] * m_tot, 6)] += w * deg / tot
            offs = dict(new)
        pattern[k] = offs
    if all(len(p) == 1 for p in pattern.values()):
        return lines
    out: List[Tuple[float, float]] = []
    for f, inten in lines:
        k = min(cluster, key=lambda c: abs(nu_hz[c] - f))
        for o, w in pattern[k].items():
            if inten * w >= min_intensity:
                out.append((f + o, inten * w))
    return out


def _auto_range(values: Sequence[float], pad: float, clamp: Optional[float] = None) -> Tuple[float, float]:
    """(min - pad, max + pad); the lower bound is raised to `clamp` only if no value
    lies within half a pad of it (so lines below the usual window are never cut off)."""
    lo = float(np.min(values)) - pad
    hi = float(np.max(values)) + pad
    if clamp is not None and float(np.min(values)) - 0.5 * pad >= clamp:
        lo = max(clamp, lo)
    return lo, hi


def simulate_quantum_transitions(
    shifts_ppm: Sequence[float],
    j_matrix: np.ndarray,
    group_counts: Optional[Sequence[int]] = None,
    freq_mhz: float = 500.0,
    j_cutoff: float = 0.3,
    max_cluster_size: int = 8,
    min_intensity: float = 1e-4,
    spins: Optional[Union[float, Sequence[float]]] = None,
    nucleus: Optional[str] = None,
    max_hilbert_dim: int = 256,
) -> List[SpinTransition]:
    """Simulate exact quantum mechanical transitions for a set of spins with arbitrary nuclear spin I.

    Parameters:
        shifts_ppm: Sequence of chemical shifts in ppm.
        j_matrix: (N, N) coupling matrix in Hz.
        group_counts: Optional sequence of number of equivalent nuclei per spin/group.
        freq_mhz: Spectrometer frequency in MHz (default: 500.0).
        j_cutoff: Minimum J (Hz) to treat two spins as coupled in the Hamiltonian.
        max_cluster_size: Maximum spin cluster size to diagonalize simultaneously.
        min_intensity: Minimum transition probability to keep.
        spins: Single spin quantum number (float) or sequence of spin quantum numbers.
        nucleus: Optional nucleus symbol (e.g. '17O', '1H', '33S') used to look up spin if spins is None.
        max_hilbert_dim: Maximum Hilbert space dimension per cluster (default: 256).

    Equivalent groups (group_counts[k] = n > 1, e.g. CH3) are expanded into n individual
    spins with the group shift and the group-group couplings (intra-group J is omitted: it
    does not affect the spectrum of magnetically equivalent spins). Clusters are partitioned
    at group level, so a group is never split. A cluster whose expanded size exceeds
    max_cluster_size or whose dimension prod((2I+1)**n) exceeds max_hilbert_dim, and that
    cannot be subdivided by raising the J threshold, is treated to FIRST ORDER: each group
    gives a multiplet nu_k + sum_j J_kj * M_j over the binomial (Pascal) distribution of
    the total M_j of the n_j spins of every coupled group j. Couplings between different
    sub-clusters (removed when a large network is split by raising the J threshold) are
    added to first order in the same way, each line being attributed to the group of its
    cluster with the closest resonance frequency.

    Intensities: each cluster is normalised to its number of nuclei (sum of counts); for an
    exactly solved cluster of like spins this gives each nucleus the same integrated weight,
    so the integrals of the groups are in the ratio of their counts.

    Returns:
        List of SpinTransition objects.
    """
    n_spins = len(shifts_ppm)
    if n_spins == 0:
        return []

    # Determine spin array
    if spins is not None:
        if isinstance(spins, (int, float)):
            spin_arr = [float(spins)] * n_spins
        else:
            spin_arr = [float(s) for s in spins]
    elif nucleus is not None:
        s_val = get_nucleus_spin(nucleus)
        spin_arr = [s_val] * n_spins
    else:
        spin_arr = [0.5] * n_spins

    count_list = [max(1, int(round(c))) for c in (group_counts if group_counts is not None else [1] * n_spins)]
    counts = np.array(count_list, dtype=float)
    nu_hz = np.array(shifts_ppm, dtype=float) * freq_mhz
    J_arr = np.array(j_matrix, dtype=float)

    def _dim(idx_list: Sequence[int]) -> int:
        d = 1
        for idx in idx_list:
            d *= int(round(2.0 * spin_arr[idx] + 1.0)) ** count_list[idx]
        return d

    def _exact_ok(idx_list: Sequence[int]) -> bool:
        return sum(count_list[i] for i in idx_list) <= max_cluster_size and _dim(idx_list) <= max_hilbert_dim

    clusters = partition_spin_network(
        n_spins,
        J_arr,
        cutoff=j_cutoff,
        max_cluster_size=max_cluster_size,
        spins=spin_arr,
        max_dim=max_hilbert_dim,
        counts=count_list,
        fallback='none',
    )
    all_transitions: List[SpinTransition] = []

    for cluster_id, cluster in enumerate(clusters):
        if len(cluster) == 1:
            # Isolated group: a single line carrying all its nuclei
            k = cluster[0]
            lines = [(float(nu_hz[k]), float(counts[k]))]
        elif not _exact_ok(cluster):
            # First-order fallback (cluster too large for exact diagonalization)
            lines = _first_order_multiplets(cluster, nu_hz, J_arr, count_list, spin_arr,
                                            j_cutoff, min_intensity)
        else:
            # Expand each group into its individual spins
            members = [k for k in cluster for _ in range(count_list[k])]
            c_size = len(members)
            c_nu = nu_hz[members]
            c_spins = [spin_arr[k] for k in members]
            c_J = np.zeros((c_size, c_size), dtype=float)
            for i, a in enumerate(members):
                for j, b in enumerate(members):
                    if a != b:
                        c_J[i, j] = J_arr[a, b]

            raw_trans = solve_spin_cluster(c_nu, c_J, spins=c_spins, min_intensity=min_intensity)

            # Scale transitions so that the sum of intensities in the cluster equals total spins (atoms)
            total_cluster_spins = float(np.sum(counts[cluster]))
            cluster_prob_sum = float(sum(prob for _, prob in raw_trans))
            c_scale = (total_cluster_spins / cluster_prob_sum) if cluster_prob_sum > 0 else 1.0
            lines = [(freq, prob * c_scale) for freq, prob in raw_trans]

        # Couplings to groups outside this cluster (dropped by the partitioning): first order
        lines = _add_external_couplings(lines, cluster, nu_hz, J_arr, count_list, spin_arr,
                                        j_cutoff, min_intensity)

        # Merge coincident lines (degenerate transitions of expanded equivalent spins)
        merged: Dict[float, float] = defaultdict(float)
        for freq, inten in lines:
            merged[round(freq, 6)] += inten
        for freq, inten in merged.items():
            all_transitions.append(
                SpinTransition(
                    freq_hz=freq,
                    intensity=inten,
                    cluster_id=cluster_id,
                    delta_ppm=freq / freq_mhz,
                )
            )

    return all_transitions


def _lorentzian_profiles(axis: np.ndarray, centers: np.ndarray, fwhm: float) -> np.ndarray:
    """Compute Lorentzian line profiles for multiple center frequencies/shifts.

    Parameters:
        axis: (P,) grid axis.
        centers: (N,) center positions.
        fwhm: Full width at half maximum.

    Returns:
        (N, P) array where profiles[i] is the normalized Lorentzian for centers[i].
    """
    if len(centers) == 0:
        return np.zeros((0, len(axis)), dtype=float)
    gamma = fwhm / 2.0
    gamma2 = gamma * gamma
    diff = axis[np.newaxis, :] - centers[:, np.newaxis]
    return (gamma / math.pi) / (diff * diff + gamma2)


def simulate_1d_spectrum(
    transitions: Sequence[SpinTransition],
    ppm_range: Optional[Tuple[float, float]] = None,
    n_points: int = 8192,
    freq_mhz: float = 500.0,
    fwhm_hz: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Convolve quantum dipole transitions with a Lorentzian lineshape to generate
    a continuous 1D NMR spectrum.

    Parameters:
        transitions: List of SpinTransition objects.
        ppm_range: (min_ppm, max_ppm) tuple. If None, auto-determined from transitions.
        n_points: Number of points on the spectrum grid.
        freq_mhz: Spectrometer frequency in MHz.
        fwhm_hz: Full width at half maximum of Lorentzian line in Hz.

    Returns:
        (ppm_axis, intensity_array)
    """
    if not transitions:
        if ppm_range is None:
            ppm_range = (-0.5, 12.0)
        axis = np.linspace(ppm_range[0], ppm_range[1], n_points)
        return axis, np.zeros(n_points)

    deltas = np.array([t.delta_ppm for t in transitions], dtype=float)
    intensities = np.array([t.intensity for t in transitions], dtype=float)

    if ppm_range is None:
        ppm_range = _auto_range(deltas, pad=1.0, clamp=-1.0)

    axis = np.linspace(ppm_range[0], ppm_range[1], n_points)
    fwhm_ppm = fwhm_hz / freq_mhz
    profs = _lorentzian_profiles(axis, deltas, fwhm_ppm)
    spectrum = intensities @ profs

    return axis, spectrum


def simulate_1d_lines(
    shifts_ppm: Sequence[float],
    intensities: Optional[Sequence[float]] = None,
    ppm_range: Optional[Tuple[float, float]] = None,
    n_points: int = 16384,
    fwhm_ppm: float = 0.04,
) -> Tuple[np.ndarray, np.ndarray]:
    """Simulate a 1D decoupled spectrum (e.g. 13C{1H} or single lines) with Lorentzian lineshapes."""
    n_lines = len(shifts_ppm)
    shifts = np.array(shifts_ppm, dtype=float)
    ints = np.array(intensities if intensities is not None else [1.0] * n_lines, dtype=float)

    if ppm_range is None:
        if n_lines > 0:
            pad = max(5.0, 0.1 * (float(np.max(shifts)) - float(np.min(shifts))))
            ppm_range = _auto_range(shifts, pad=pad, clamp=-5.0)
        else:
            ppm_range = (0.0, 220.0)

    ppm_axis = np.linspace(ppm_range[0], ppm_range[1], n_points)
    profs = _lorentzian_profiles(ppm_axis, shifts, fwhm_ppm)
    spectrum = ints @ profs if n_lines > 0 else np.zeros(n_points, dtype=float)
    return ppm_axis, spectrum
