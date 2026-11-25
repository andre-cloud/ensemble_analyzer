"""Validation of the ``nmr`` protocol key and cross-step resolution."""

from unittest.mock import MagicMock

import pytest

from ensemble_analyzer.protocol.nmr import check_protocols
from ensemble_analyzer.protocol.protocol import Protocol

TMS = {"name": "TMS", "nuclei": ["1H", "13C"]}


def _p(number, **kw):
    base = dict(functional="PBE0", basis="pcSseg-1", calculator="orca")
    base.update(kw)
    return Protocol(number=number, **base)


def _nmr(**kw):
    base = {"shielding": True, "nuclei": ["1H", "13C"], "references": [TMS]}
    base.update(kw)
    return base


# --- normalize_settings ---------------------------------------------------

def test_defaults_and_normalisation():
    s = _p(1, nmr=_nmr(nuclei=["1h", "13-C"], references=[{"name": "tms", "nuclei": ["1H", "13c"]}])).nmr
    assert s["nuclei"] == ["1H", "13C"]
    assert s["references"][0]["nuclei"] == ["1H", "13C"]
    assert s["coupling"] is False and s["coupling_terms"] == "all"
    assert s["frequency"] == 400.0 and s["fwhm"] == 1.0
    assert s["atoms"] is None and s["equivalent"] == [] and s["decouple"] == {}
    assert s["reference_opt_step"] is None and s["coupling_step"] is None


def test_off_values():
    assert _p(1, nmr=None).nmr is None
    assert _p(1, nmr=False).nmr is None


def test_coupling_only_needs_no_reference():
    s = _p(1, nmr={"coupling": True, "nuclei": ["1H", "19F"]}).nmr
    assert s["shielding"] is False and s["references"] == []


@pytest.mark.parametrize("nmr,match", [
    (True, "dictionary"),
    ({"nuclei": ["1H"]}, "shielding.*coupling"),
    (_nmr(nuclei=[]), "nuclei"),
    (_nmr(nuclei=["H"]), "isotope"),
    (_nmr(nuclei=["1H", "2H"], references=[{"name": "TMS", "nuclei": ["1H", "2H"]}]), "same element"),
    (_nmr(foo=1), "unknown"),
    (_nmr(references=[]), "reference"),
    (_nmr(references=[{"name": "TMS", "nuclei": ["1H"]}]), "13C.*no reference"),
    (_nmr(references=[TMS, {"name": "PhF", "smiles": "Fc1ccccc1", "nuclei": ["13C"]}]), "13C.*more than one"),
    (_nmr(references=[{"name": "TMS", "nuclei": ["1H", "13C", "29Si"]}]), "29Si.*not requested"),
    (_nmr(references=[{"name": "nothing", "nuclei": ["1H", "13C"]}]), "xyz.*smiles"),
    (_nmr(references=[{"name": "TMS", "nuclei": ["1H", "13C"], "bad": 1}]), "unknown"),
    (_nmr(references=[{"nuclei": ["1H", "13C"], "smiles": "C"}]), "name"),
    (_nmr(atoms=[0, -1]), "atoms"),
    (_nmr(equivalent=[[1]]), "equivalent"),
    (_nmr(decouple={"19F": ["1H"]}), "decouple"),
    (_nmr(coupling_terms="psO"), "coupling_terms"),
    (_nmr(frequency=0), "frequency"),
    (_nmr(fwhm={"1H": -1}), "fwhm"),
])
def test_rejected(nmr, match):
    with pytest.raises(ValueError, match=match):
        _p(1, nmr=nmr)


@pytest.mark.parametrize("kw,match", [
    (dict(opt=True), "single point"),
    (dict(freq=True), "single point"),
    (dict(nroots=5), "TD-DFT"),
    (dict(deltascf={}), "deltascf"),
    (dict(calculator="tblite", functional="GFN2-xTB", basis=""), "orca"),
    (dict(add_input="NMR"), "add_input"),
    (dict(add_input="%eprnmr nuclei = all H {shift} end"), "add_input"),
])
def test_step_constraints(kw, match):
    with pytest.raises(ValueError, match=match):
        _p(1, nmr=_nmr(), **kw)


def test_fwhm_per_nucleus_and_decouple():
    s = _p(1, nmr=_nmr(fwhm={"1h": 0.5, "13C": 2}, decouple={"13C": ["1h"]})).nmr
    assert s["fwhm"] == {"1H": 0.5, "13C": 2.0}
    assert s["decouple"] == {"13C": ["1H"]}


def test_settings_survive_reload():
    p = _p(1, nmr=_nmr(atoms=[0, 3], equivalent=[[4, 5, 6]]))
    q = Protocol(**p.__dict__)
    assert q.nmr == p.nmr


# --- check_protocols ------------------------------------------------------

def test_reference_opt_step_is_last_ground_state_opt():
    steps = [
        _p(0, calculator="tblite", functional="GFN2-xTB", basis="", opt=True),
        _p(1, opt=True, freq=True),
        _p(2, nroots=5, opt=True),                 # TD-DFT opt: excluded
        _p(3, ts=True),                            # TS: excluded
        _p(4, opt=True, deltascf={}),              # ΔSCF: excluded
        _p(5, nmr=_nmr()),
    ]
    check_protocols(steps, MagicMock())
    assert steps[5].nmr["reference_opt_step"] == 1


def test_ml_optimisation_counts():
    steps = [_p(0, calculator="tblite", functional="GFN2-xTB", basis="", opt=True), _p(1, nmr=_nmr())]
    check_protocols(steps, MagicMock())
    assert steps[1].nmr["reference_opt_step"] == 0


def test_no_optimisation_warns():
    log = MagicMock()
    steps = [_p(0), _p(1, nmr=_nmr())]
    check_protocols(steps, log)
    assert steps[1].nmr["reference_opt_step"] is None
    assert "level of the NMR step" in log.warning.call_args[0][0]


def test_opt_from_override():
    steps = [_p(0, opt=True), _p(1, opt=True), _p(2, nmr=_nmr(opt_from=0))]
    check_protocols(steps, MagicMock())
    assert steps[2].nmr["reference_opt_step"] == 0


@pytest.mark.parametrize("opt_from", [3, 1, 5])
def test_opt_from_invalid(opt_from):
    steps = [_p(0, opt=True), _p(1, ts=True), _p(2, opt=True), _p(3, nmr=_nmr(opt_from=opt_from))]
    with pytest.raises(ValueError, match="opt_from"):
        check_protocols(steps, MagicMock())


def test_coupling_step_auto_and_disabled():
    steps = [
        _p(0, opt=True),
        _p(1, nmr={"coupling": True, "nuclei": ["1H", "19F"]}),
        _p(2, nmr=_nmr(nuclei=["1H"], references=[{"name": "TMS", "nuclei": ["1H"]}])),
        _p(3, nmr=_nmr(nuclei=["1H"], references=[{"name": "TMS", "nuclei": ["1H"]}], coupling_from=False)),
    ]
    check_protocols(steps, MagicMock())
    assert steps[2].nmr["coupling_step"] == 1
    assert steps[3].nmr["coupling_step"] is None
    assert steps[1].nmr["coupling_step"] is None


def test_coupling_step_explicit_must_have_coupling():
    steps = [_p(0, nmr=_nmr()), _p(1, nmr=_nmr(coupling_from=0))]
    with pytest.raises(ValueError, match="coupling_from"):
        check_protocols(steps, MagicMock())


def test_same_name_must_be_same_molecule():
    steps = [
        _p(0, nmr=_nmr(nuclei=["13C"], references=[{"name": "ref", "smiles": "C", "nuclei": ["13C"]}])),
        _p(1, nmr=_nmr(nuclei=["13C"], references=[{"name": "ref", "smiles": "CC", "nuclei": ["13C"]}])),
    ]
    with pytest.raises(ValueError, match="ref"):
        check_protocols(steps, MagicMock())


def test_missing_rdkit_only_for_smiles(monkeypatch):
    import ensemble_analyzer.protocol.nmr as mod
    monkeypatch.setattr(mod, "rdkit_available", lambda: False)
    check_protocols([_p(0, nmr=_nmr())], MagicMock())          # bundled TMS: fine
    steps = [_p(0, nmr=_nmr(nuclei=["13C"], references=[{"name": "PhF", "smiles": "Fc1ccccc1", "nuclei": ["13C"]}]))]
    with pytest.raises(ImportError, match="rdkit"):
        check_protocols(steps, MagicMock())
