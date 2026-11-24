"""Ensemble NMR analysis: Boltzmann averages, equivalent groups, shifts, tables and spectra."""

import json
from unittest.mock import MagicMock

import numpy as np
import pytest

from ensemble_analyzer._nmr.analysis import average_nmr, process_nmr_step
from ensemble_analyzer._nmr.nuclei import larmor_mhz, shielding_to_shift
from ensemble_analyzer.conformer.conformer import Conformer
from ensemble_analyzer.conformer.energy_data import EnergyRecord
from ensemble_analyzer.conformer.nmr_data import NMRRecord
from ensemble_analyzer.protocol.nmr import check_protocols
from ensemble_analyzer.protocol.protocol import Protocol

# CH3F: C, F, 3 H (methyl rotor)
ATOMS = ["C", "F", "H", "H", "H"]
GEOM = np.array([[0, 0, 0], [0, 0, 1.38], [1.03, 0, -0.36], [-0.515, 0.892, -0.36], [-0.515, -0.892, -0.36]])
REFS = {"1H": {"sigma": 31.0, "reference": "TMS", "delta": 0.0, "n_atoms": 12, "opt_step": 0},
        "13C": {"sigma": 180.0, "reference": "TMS", "delta": 0.0, "n_atoms": 4, "opt_step": 0},
        "19F": {"sigma": 190.0, "reference": "CFCl3", "delta": 0.0, "n_atoms": 1, "opt_step": 0}}


def _p(number, **nmr):
    s = {"shielding": True, "coupling": True, "nuclei": ["1H", "13C", "19F"],
         "references": [{"name": "TMS", "nuclei": ["1H", "13C"]}, {"name": "CFCl3", "nuclei": ["19F"]}]}
    s.update(nmr)
    if not s["shielding"]:
        s["references"] = []
    return Protocol(number=number, functional="PBE0", basis="def2-SVP", calculator="orca", nmr=s)


def _conf(number, pop, step, shield, coupling):
    c = Conformer(number=number, geom=GEOM.copy(), atoms=ATOMS, raw=True)
    c.energies.add(step, EnergyRecord(E=-100.0, Pop=pop))
    c.nmr_data.add(step, NMRRecord(shielding=shield, coupling=coupling))
    return c


def _ensemble(step=1):
    a = _conf(1, 75.0, step, {0: 100.0, 1: 400.0, 2: 27.0, 3: 27.2, 4: 27.4},
              {(1, 2): 46.0, (1, 3): 46.0, (1, 4): 46.0, (0, 1): -160.0, (0, 2): 150.0, (0, 3): 150.0, (0, 4): 150.0})
    b = _conf(2, 25.0, step, {0: 104.0, 1: 404.0, 2: 27.4, 3: 27.6, 4: 27.8},
              {(1, 2): 50.0, (1, 3): 50.0, (1, 4): 50.0, (0, 1): -164.0, (0, 2): 154.0, (0, 3): 154.0, (0, 4): 154.0})
    return [a, b]


def _peaks(x, y, floor=0.5):
    inner = (y[1:-1] > y[:-2]) & (y[1:-1] > y[2:]) & (y[1:-1] > floor)
    return x[1:-1][inner]


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_boltzmann_averages_and_groups():
    r = average_nmr(_ensemble(), _p(1), REFS)
    assert r.sigma[0] == pytest.approx(0.75 * 100 + 0.25 * 104)
    assert r.j[(1, 2)] == pytest.approx(0.75 * 46 + 0.25 * 50)
    h = [g for g in r.groups if g.nucleus == "1H"]
    assert len(h) == 1 and h[0].members == (2, 3, 4)
    sig_h = np.mean([0.75 * 27.0 + 0.25 * 27.4, 0.75 * 27.2 + 0.25 * 27.6, 0.75 * 27.4 + 0.25 * 27.8])
    assert h[0].delta == pytest.approx(shielding_to_shift(sig_h, 31.0))
    assert r.weights == {1: 0.75, 2: 0.25}


def test_inactive_conformers_are_skipped():
    confs = _ensemble()
    confs[1].active = False
    r = average_nmr(confs, _p(1), REFS)
    assert r.sigma[0] == pytest.approx(100.0) and r.weights == {1: 1.0}


def test_user_equivalence_never_mixes_elements():
    r = average_nmr(_ensemble(), _p(1, equivalent=[[0, 1]]), REFS)
    assert any(g.members == (0,) for g in r.groups) and any(g.members == (1,) for g in r.groups)


def test_coupling_from_previous_step():
    steps = [_p(0, shielding=False), _p(1, coupling=False)]
    check_protocols(steps, MagicMock())
    confs = _ensemble(step=0)
    for c in confs:                                           # step 1: shielding only, 50/50
        c.nmr_data.add(1, NMRRecord(shielding={0: 100.0, 1: 400.0, 2: 27.0, 3: 27.0, 4: 27.0}))
        c.energies.add(1, EnergyRecord(E=-100.0, Pop=50.0))
    r = average_nmr(confs, steps[1], REFS)
    assert r.j[(1, 2)] == pytest.approx(48.0)
    assert r.coupling_step == 0


def test_process_writes_tables_and_spectra(workdir):
    process_nmr_step(_ensemble(), _p(1), REFS, MagicMock())
    table = (workdir / "nmr_p1.out").read_text()
    assert "TMS" in table and "δ" in table and "J" in table
    data = json.loads((workdir / "nmr_p1.json").read_text())
    assert set(data) >= {"atoms", "groups", "couplings", "references", "weights", "spectra"}
    for nuc in ("1H", "13C", "19F"):
        assert (workdir / f"NMR_{nuc}_p1_comp.xy").exists()
    x, y = np.loadtxt(workdir / "NMR_1H_p1_comp.xy", unpack=True)
    assert y.max() == pytest.approx(1.0)
    peaks = _peaks(x, y)                                      # 1H: doublet from 19F
    assert len(peaks) == 2
    assert abs(peaks[0] - peaks[1]) * larmor_mhz("1H", 400) == pytest.approx(47.0, abs=0.5)


def test_shielding_only_gives_singlets(workdir):
    steps = [_p(1, coupling=False, coupling_from=False)]
    check_protocols(steps, MagicMock())
    process_nmr_step(_ensemble(), steps[0], REFS, MagicMock())
    x, y = np.loadtxt(workdir / "NMR_1H_p1_comp.xy", unpack=True)
    assert len(_peaks(x, y)) == 1


def test_coupling_only_gives_tables_only(workdir):
    process_nmr_step(_ensemble(), _p(1, shielding=False), {}, MagicMock())
    assert (workdir / "nmr_p1.out").exists()
    assert not list(workdir.glob("NMR_*_comp.xy"))


def test_experimental_reference_shift_fitted(workdir):
    process_nmr_step(_ensemble(), _p(1), REFS, MagicMock())
    x, y = np.loadtxt(workdir / "NMR_19F_p1_comp.xy", unpack=True)
    np.savetxt(workdir / "nmr_19f_ref.dat", np.column_stack((x + 0.8, y)))
    process_nmr_step(_ensemble(step=2), _p(2), REFS, MagicMock())
    result = json.loads((workdir / "nmr_p2.json").read_text())["spectra"]["19F"]
    assert result["shift_ppm"] == pytest.approx(0.8, abs=0.05) and result["similarity"] > 90
