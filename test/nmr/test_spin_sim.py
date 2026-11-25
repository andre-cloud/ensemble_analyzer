"""Spin simulation tests, ported from FROGChem (tests/unit/test_nmr_sim.py, test_math_review_nmr.py)."""

import math
from collections import defaultdict

import numpy as np
import pytest

from ensemble_analyzer._nmr.spin_sim import (
    SpinTransition,
    _auto_range,
    partition_spin_network,
    simulate_1d_lines,
    simulate_1d_spectrum,
    simulate_quantum_transitions,
    solve_spin_cluster,
)


def _aggregate(transitions, decimals=1):
    """Sum intensities of lines that coincide to `decimals` Hz."""
    agg = defaultdict(float)
    for t in transitions:
        agg[round(t.freq_hz, decimals)] += t.intensity
    return dict(agg)


def _bins(transitions, gap=0.5):
    """Merge lines closer than `gap` Hz; returns [(intensity-weighted centre, total intensity)]."""
    out = []
    for f, i in sorted((t.freq_hz, t.intensity) for t in transitions):
        if out and f - out[-1][2] < gap:
            c, s, _ = out[-1]
            out[-1] = ((c * s + f * i) / (s + i), s + i, f)
        else:
            out.append((f, i, f))
    return [(c, s) for c, s, _ in out]


def _ab_analytic(nu_a, nu_b, J):
    d = nu_a - nu_b
    C = math.sqrt(d * d + J * J)
    c = 0.5 * (nu_a + nu_b)
    s = J / C
    return sorted([(c + J / 2 + C / 2, 1 - s), (c - J / 2 + C / 2, 1 + s),
                   (c + J / 2 - C / 2, 1 + s), (c - J / 2 - C / 2, 1 - s)])



def test_spin_hamiltonian_ab_system_analytical():
    """Verify that exact quantum diagonalization matches analytical AB system formulas:
    Positions: nu0 +/- C/2 +/- J/2
    Intensities: 1 +/- J/C (roof effect)
    """
    nu_A = 200.0
    nu_B = 220.0
    delta_nu = nu_B - nu_A  # 20 Hz
    J_val = 10.0            # 10 Hz
    nu_0 = (nu_A + nu_B) / 2.0  # 210 Hz

    C = math.sqrt(delta_nu**2 + J_val**2)  # sqrt(500) = 22.36068 Hz

    expected_freqs = sorted([
        nu_0 - C / 2.0 - J_val / 2.0,
        nu_0 - C / 2.0 + J_val / 2.0,
        nu_0 + C / 2.0 - J_val / 2.0,
        nu_0 + C / 2.0 + J_val / 2.0,
    ])
    expected_intensities = [
        1.0 - J_val / C,
        1.0 + J_val / C,
        1.0 + J_val / C,
        1.0 - J_val / C,
    ]

    nu = np.array([nu_A, nu_B])
    J = np.array([[0.0, J_val], [J_val, 0.0]])

    transitions = solve_spin_cluster(nu, J)
    transitions.sort(key=lambda x: x[0])

    assert len(transitions) == 4
    calc_freqs = [t[0] for t in transitions]
    calc_intensities = [t[1] for t in transitions]

    np.testing.assert_allclose(calc_freqs, expected_freqs, atol=1e-4)
    np.testing.assert_allclose(calc_intensities, expected_intensities, atol=1e-4)


def test_spin_hamiltonian_identical_spins():
    """For two magnetically equivalent spins (A2 system), mutual J coupling
    produces zero splitting and the spectrum collapses to a single singlet.
    """
    nu = np.array([150.0, 150.0])
    J = np.array([[0.0, 12.0], [12.0, 0.0]])

    transitions = solve_spin_cluster(nu, J)
    # Both allowed transitions in the triplet manifold (-1 -> 0 and 0 -> +1)
    # occur at the exact same frequency 150.0 Hz (zero splitting between equivalent spins).
    assert len(transitions) == 2
    for freq, prob in transitions:
        assert np.isclose(freq, 150.0)
        assert prob > 0.0


def test_1d_simulation_spectrum():
    nu = np.array([100.0, 200.0])  # in Hz, corresponding to 0.2 and 0.4 ppm at 500 MHz
    J = np.array([[0.0, 0.0], [0.0, 0.0]])
    transitions = simulate_quantum_transitions(
        shifts_ppm=[0.2, 0.4],
        j_matrix=J,
        freq_mhz=500.0,
    )
    ppm_axis, intensity = simulate_1d_spectrum(
        transitions,
        ppm_range=(0.0, 0.6),
        n_points=1000,
        freq_mhz=500.0,
        fwhm_hz=2.0,
    )
    assert len(ppm_axis) == 1000
    assert len(intensity) == 1000
    assert np.all(intensity >= 0.0)

    # Peaks should be centered near 0.2 and 0.4 ppm
    peak1_idx = np.argmin(np.abs(ppm_axis - 0.2))
    peak2_idx = np.argmin(np.abs(ppm_axis - 0.4))
    assert intensity[peak1_idx] > 10.0
    assert intensity[peak2_idx] > 10.0


def test_simulate_1d_lines():
    c_shifts = [25.0, 55.0, 130.0, 175.0]
    ppm, intensity = simulate_1d_lines(c_shifts, ppm_range=(0.0, 200.0), n_points=1000, fwhm_ppm=1.0)
    assert len(ppm) == 1000
    assert len(intensity) == 1000
    assert np.all(intensity >= 0.0)

    for cs in c_shifts:
        idx = np.argmin(np.abs(ppm - cs))
        assert intensity[idx] > 0.1


def test_solve_spin_cluster_oxygen17_equivalent():
    """Two magnetically equivalent Oxygen-17 spins (I = 5/2) should collapse
    into a single unsplit resonance frequency at nu0.
    """
    nu = np.array([250.0, 250.0])
    J = np.array([[0.0, 15.0], [15.0, 0.0]])
    trans = solve_spin_cluster(nu, J, spins=2.5)
    assert len(trans) > 0
    # Every allowed transition must occur at exactly 250.0 Hz
    for freq, prob in trans:
        assert np.isclose(freq, 250.0, atol=1e-3)
        assert prob > 0.0


def test_solve_spin_cluster_oxygen17_coupled_ab():
    """Two non-equivalent Oxygen-17 spins (I = 5/2) coupled together."""
    nu = np.array([200.0, 300.0])
    J = np.array([[0.0, 20.0], [20.0, 0.0]])
    trans = solve_spin_cluster(nu, J, spins=2.5)
    assert len(trans) > 10
    # Center of gravity of transitions should match average frequency
    weighted_freq = sum(f * p for f, p in trans) / sum(p for _, p in trans)
    assert np.isclose(weighted_freq, 250.0, atol=1e-2)


def test_solve_spin_cluster_mixed_spins():
    """Solve spin cluster with mixed spins: spin 1/2 (e.g. 1H) coupled to spin 5/2 (e.g. 17O)."""
    nu = np.array([500.0, 100.0])
    J = np.array([[0.0, 12.0], [12.0, 0.0]])
    trans = solve_spin_cluster(nu, J, spins=[0.5, 2.5])
    assert len(trans) > 0
    # Hilbert space is 2 * 6 = 12 states
    for freq, prob in trans:
        assert freq > 0.0
        assert prob > 0.0


def test_partition_spin_network_high_spin():
    """Verify that partition_spin_network respects max_dim for high-spin nuclei."""
    # 4 spins of spin 5/2: Hilbert space is 6^4 = 1296 > 256.
    # The partitioner must subdivide this into smaller clusters.
    J = np.array([
        [0.0, 10.0, 1.0, 0.5],
        [10.0, 0.0, 0.5, 1.0],
        [1.0, 0.5, 0.0, 10.0],
        [0.5, 1.0, 10.0, 0.0],
    ])
    spins = [2.5, 2.5, 2.5, 2.5]
    clusters = partition_spin_network(4, J, cutoff=0.3, max_dim=256, spins=spins)
    for c in clusters:
        dim = 1
        for idx in c:
            dim *= int(round(2.0 * spins[idx] + 1.0))
        assert dim <= 256, f"Cluster {c} exceeded max_dim 256 with dim {dim}"


def test_simulate_quantum_transitions_oxygen():
    """Verify simulate_quantum_transitions works automatically with nucleus='17O'."""
    shifts = [10.0, 20.0]
    J = np.array([[0.0, 15.0], [15.0, 0.0]])
    transitions = simulate_quantum_transitions(
        shifts_ppm=shifts,
        j_matrix=J,
        freq_mhz=500.0,
        nucleus='17O',
    )
    assert len(transitions) > 0
    # Total intensity must integrate to exactly 2.0 (2 atoms)
    total_intensity = sum(t.intensity for t in transitions)
    assert np.isclose(total_intensity, 2.0, atol=1e-3)


def test_group_expansion_ethyl_matches_explicit_spins():
    J = np.array([[0.0, 7.0], [7.0, 0.0]])
    grouped = simulate_quantum_transitions([1.2, 3.6], J, group_counts=[3, 2], freq_mhz=500.0)

    J5 = np.zeros((5, 5))
    for i in range(3):
        for j in range(3, 5):
            J5[i, j] = J5[j, i] = 7.0
    explicit = simulate_quantum_transitions([1.2] * 3 + [3.6] * 2, J5, freq_mhz=500.0)

    a, b = _aggregate(grouped, 3), _aggregate(explicit, 3)
    assert set(a) == set(b)
    for k in a:
        assert abs(a[k] - b[k]) < 1e-3

    # CH3: 1:2:1 triplet (J = 7 Hz) around 600 Hz; CH2: 1:3:3:1 quartet around 1800 Hz
    # (second order at J/dnu = 0.006 shifts lines by < 0.1 Hz and intensities by < 2 %)
    ch3 = _bins([t for t in grouped if t.freq_hz < 1000])
    ch2 = _bins([t for t in grouped if t.freq_hz > 1000])
    np.testing.assert_allclose([c for c, _ in ch3], [593.0, 600.0, 607.0], atol=0.1)
    np.testing.assert_allclose([s for _, s in ch3], [0.75, 1.5, 0.75], atol=0.02)
    np.testing.assert_allclose([c for c, _ in ch2], [1789.5, 1796.5, 1803.5, 1810.5], atol=0.1)
    np.testing.assert_allclose([s for _, s in ch2], [0.25, 0.75, 0.75, 0.25], atol=0.02)


def test_group_integrals_follow_counts_when_coupled():
    for J_val in (0.0, 7.0):
        tr = simulate_quantum_transitions([1.0, 3.0], np.array([[0.0, J_val], [J_val, 0.0]]),
                                          group_counts=[3, 1], freq_mhz=500.0)
        s1 = sum(t.intensity for t in tr if t.delta_ppm < 2.0)
        s2 = sum(t.intensity for t in tr if t.delta_ppm > 2.0)
        assert abs(s1 - 3.0) < 1e-3 and abs(s2 - 1.0) < 1e-3


def test_first_order_fallback_pascal_multiplets():
    J = np.array([[0.0, 7.0], [7.0, 0.0]])
    tr = simulate_quantum_transitions([1.2, 3.6], J, group_counts=[3, 2], freq_mhz=500.0,
                                      max_cluster_size=4)   # 5 spins > 4 -> first order
    ch3 = _aggregate([t for t in tr if t.freq_hz < 1000], decimals=3)
    ch2 = _aggregate([t for t in tr if t.freq_hz > 1000], decimals=3)
    assert ch3 == {593.0: 0.75, 600.0: 1.5, 607.0: 0.75}
    assert ch2 == {1789.5: 0.25, 1796.5: 0.75, 1803.5: 0.75, 1810.5: 0.25}


def test_first_order_fallback_oversize_group_binomial():
    # tBu-like group of 9 equivalent spins (dim 512 > 256) coupled to one spin: not exactly
    # solvable and not splittable -> first order: X is a 10-line binomial C(9,k)/2^9 multiplet.
    J = np.array([[0.0, 1.0], [1.0, 0.0]])
    tr = simulate_quantum_transitions([1.0, 5.0], J, group_counts=[9, 1], freq_mhz=500.0)
    x = sorted((t.freq_hz, t.intensity) for t in tr if t.freq_hz > 1000)
    t9 = sorted((t.freq_hz, t.intensity) for t in tr if t.freq_hz < 1000)
    np.testing.assert_allclose([f for f, _ in x], 2500.0 + np.arange(-4.5, 4.6, 1.0), atol=1e-9)
    np.testing.assert_allclose([i for _, i in x], [math.comb(9, k) / 512 for k in range(10)], atol=1e-12)
    np.testing.assert_allclose(t9, [(499.5, 4.5), (500.5, 4.5)], atol=1e-9)


def test_negative_frequency_transitions_kept():
    nu = np.array([-1000.0, -750.0])
    J = np.array([[0.0, 7.0], [7.0, 0.0]])
    tr = sorted(solve_spin_cluster(nu, J))
    ref = _ab_analytic(-1000.0, -750.0, 7.0)
    assert len(tr) == 4
    np.testing.assert_allclose([t[0] for t in tr], [r[0] for r in ref], atol=1e-8)
    np.testing.assert_allclose([t[1] for t in tr], [r[1] for r in ref], atol=1e-8)

    mixed = simulate_quantum_transitions([-2.0, 1.5], J, freq_mhz=500.0)
    neg = sum(t.intensity for t in mixed if t.delta_ppm < 0)
    pos = sum(t.intensity for t in mixed if t.delta_ppm > 0)
    assert abs(neg - 1.0) < 1e-3 and abs(pos - 1.0) < 1e-3


def test_auto_range_includes_negative_shifts():
    tr = [SpinTransition(freq_hz=-1500.0, intensity=1.0, delta_ppm=-3.0)]
    ax, inten = simulate_1d_spectrum(tr, freq_mhz=500.0, fwhm_hz=1.0)
    assert ax[0] < -3.0 < ax[-1]
    assert inten.max() > 0.9 * 1.0 / (math.pi * 0.001)   # Lorentzian peak height 1/(pi*gamma)

    ax, inten = simulate_1d_lines([-120.0, -60.0])
    assert ax[0] < -120.0 and ax[-1] > -60.0

    # Usual 1H window is unchanged when nothing is below the clamp
    ax, _ = simulate_1d_spectrum([SpinTransition(freq_hz=3500.0, intensity=1.0, delta_ppm=7.0),
                                  SpinTransition(freq_hz=500.0, intensity=1.0, delta_ppm=1.0)], freq_mhz=500.0)
    assert ax[0] == 0.0 and ax[-1] == 8.0
