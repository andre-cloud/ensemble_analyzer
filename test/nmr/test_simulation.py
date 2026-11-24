"""Spectrum of one observed nucleus with homo- and heteronuclear partners."""

import numpy as np
import pytest

from ensemble_analyzer._nmr.nuclei import abundance, larmor_mhz
from ensemble_analyzer._nmr.simulation import SpinGroup, simulate_nucleus, stick_spectrum


def _lines(transitions, decimals=2):
    out = {}
    for t in transitions:
        key = round(t.freq_hz, decimals)
        out[key] = out.get(key, 0.0) + t.intensity
    return dict(sorted(out.items()))


def _j(n, pairs):
    J = np.zeros((n, n))
    for (a, b), v in pairs.items():
        J[a, b] = J[b, a] = v
    return J


def test_ch3f_proton_doublet_from_fluorine():
    groups = [SpinGroup("1H", 4.0, 3), SpinGroup("19F", -270.0, 1)]
    tr = simulate_nucleus("1H", groups, _j(2, {(0, 1): 47.0}), frequency_1h=400)
    lines = _lines(tr)
    nu = 4.0 * 400
    assert list(lines) == [round(nu - 23.5, 2), round(nu + 23.5, 2)]
    assert np.allclose(list(lines.values()), [1.5, 1.5])


def test_carbon_quartet_from_protons_and_no_proton_satellites():
    groups = [SpinGroup("13C", 20.0, 1), SpinGroup("1H", 1.0, 3)]
    J = _j(2, {(0, 1): 125.0})
    c = _lines(simulate_nucleus("13C", groups, J, frequency_1h=400))
    assert np.allclose(np.array(list(c.values())) / sum(c.values()), [1 / 8, 3 / 8, 3 / 8, 1 / 8])
    h = _lines(simulate_nucleus("1H", groups, J, frequency_1h=400))
    assert len(h) == 1 and np.isclose(sum(h.values()), 3.0)


def test_low_abundance_homonuclear_couplings_ignored():
    groups = [SpinGroup("13C", 20.0, 1), SpinGroup("13C", 60.0, 1)]
    tr = simulate_nucleus("13C", groups, _j(2, {(0, 1): 35.0}), frequency_1h=400)
    assert len(_lines(tr)) == 2


def test_decoupled_partner_gives_singlet():
    groups = [SpinGroup("13C", 20.0, 1), SpinGroup("1H", 1.0, 3)]
    tr = simulate_nucleus("13C", groups, _j(2, {(0, 1): 125.0}), frequency_1h=400, decouple=["1H"])
    assert len(_lines(tr)) == 1


def test_partial_abundance_gives_satellites():
    groups = [SpinGroup("1H", 0.0, 1), SpinGroup("119Sn", 0.0, 1)]
    tr = simulate_nucleus("1H", groups, _j(2, {(0, 1): 50.0}), frequency_1h=400)
    lines = _lines(tr)
    p = abundance("119Sn")
    assert np.allclose(list(lines.values()), [p / 2, 1 - p, p / 2])
    assert list(lines) == [-25.0, 0.0, 25.0]


def test_second_order_ab_kept_with_heteronuclear_partner():
    # AB proton pair, one of them coupled to 19F: each F state shifts the AB system
    groups = [SpinGroup("1H", 4.00, 1), SpinGroup("1H", 4.02, 1), SpinGroup("19F", -200.0, 1)]
    J = _j(3, {(0, 1): 15.0, (0, 2): 10.0})
    lines = _lines(simulate_nucleus("1H", groups, J, frequency_1h=400))
    assert len(lines) == 8                                   # two AB quartets
    assert np.isclose(sum(lines.values()), 2.0)
    assert not np.allclose(list(lines.values()), 0.25)       # roof effect, not first order


def test_stick_spectrum_in_ppm():
    tr = simulate_nucleus("13C", [SpinGroup("13C", 77.0, 1)], np.zeros((1, 1)), frequency_1h=400)
    ppm, inten = stick_spectrum(tr, "13C", 400)
    assert np.allclose(ppm, [77.0]) and np.allclose(inten, [1.0])


def test_larmor_frequency_scales_with_gamma():
    assert larmor_mhz("13C", 400) == pytest.approx(100.6, abs=0.1)
