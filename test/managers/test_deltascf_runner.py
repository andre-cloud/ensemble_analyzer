"""
ΔSCF runner: input of the sub-calculations built from the step calculator.
The full workflow is validated end-to-end with ORCA 6.1.1 and NWChem 7.3.1
(see docs/source/user_guide/deltascf.md); here only the generated keywords.
"""

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from ensemble_analyzer._managers.deltascf_runner import (
    DeltaSCFRunner, _NWChemEngine, _OrcaEngine, _Reference, _Target,
)
from ensemble_analyzer.protocol import deltascf as ds
from ensemble_analyzer.protocol.deltascf import Excitation
from ensemble_analyzer.protocol.protocol import Protocol


@pytest.fixture
def conf():
    return SimpleNamespace(folder="conf_1", number=1, atoms=["C", "O", "H", "H"])


def _orca(protocol, conf):
    return _OrcaEngine(protocol, conf, 1)


def _orca_build(engine, **kw):
    with patch("enan_calculators._orca.shutil.which", return_value="orca"), \
         patch("ensemble_analyzer.calculators.orca.OrcaCalc.VERSION", 6):
        calc = engine.build(Path("run"), **kw)
    return calc.parameters["orcasimpleinput"], calc.parameters["orcablocks"], calc


def test_orca_attempt_keywords(conf, tmp_path):
    p = Protocol(number=2, functional="PBE", basis="def2-SVP", read_orbitals="1", deltascf={})
    simple, blocks, calc = _orca_build(
        _orca(p, conf), mode="protocol", guess=tmp_path / "ref.gbw",
        target=_Target([Excitation("alpha", 1, 0)]), attempt=ds.ORCA_LADDER[0])
    assert "UKS DeltaSCF GMF MORead" in simple
    assert f'%moinp "{(tmp_path / "ref.gbw").resolve()}"' in blocks
    assert "%scf alphaconf 0,1,1 ConvForced true SOSCFUpdateSPOEst false MaxIter 200 end" in blocks
    assert blocks.count("%moinp") == 1          # the step read_orbitals is not added again
    assert Path(calc.directory) == Path("run")


def test_orca_reference_is_tddft_with_iroot0(conf):
    p = Protocol(number=2, functional="PBE", basis="def2-SVP", deltascf={}, opt=True)
    simple, blocks, _ = _orca_build(_orca(p, conf), mode="sp", nroots=5)
    assert "%tddft nroots 5 tda true iroot 0 end" in blocks
    assert "opt" not in simple.lower()
    assert "DeltaSCF" not in simple


def test_orca_triplet_and_quartet_are_plain_uks(conf, tmp_path):
    p = Protocol(number=2, functional="PBE", basis="def2-SVP", deltascf={})
    engine = _orca(p, conf)
    simple, _, calc = _orca_build(engine, mult=3, guess=tmp_path / "mix.gbw")
    assert "UKS" in simple and "DeltaSCF" not in simple and calc.parameters["mult"] == 3
    simple, blocks, _ = _orca_build(engine, mult=4, guess=tmp_path / "ref.gbw", target=_Target([]),
                                    attempt=ds.ORCA_LADDER[0])
    assert "DeltaSCF" not in simple and "%scf" not in blocks


def test_orca_two_spin_determinant(conf):
    p = Protocol(number=2, functional="PBE", basis="def2-SVP", charge=1, mult=2, deltascf={})
    target = _Target([Excitation("alpha", 0, 0), Excitation("beta", 0, 0)])   # recoupling M2
    _, blocks, _ = _orca_build(_orca(p, conf), target=target, attempt=ds.ORCA_LADDER[0])
    assert "alphaconf 0,1 betaconf 0,1" in blocks


def test_orca_raw_override(conf):
    p = Protocol(number=2, functional="PBE", basis="def2-SVP", deltascf={"alphaconf": "0,1,1"})
    _, blocks, _ = _orca_build(_orca(p, conf), target=_Target(alphaconf="0,1,1"), attempt=ds.ORCA_LADDER[2])
    assert "%scf alphaconf 0,1,1 ConvForced true MaxIter 300 end" in blocks


@patch("ensemble_analyzer.calculators.nwchem.NWCHEM_COMMAND", "nwchem")
def test_nwchem_attempt_vectors_and_swaps(conf, tmp_path):
    p = Protocol(number=2, functional="PBE", basis="def2-SVP", calculator="nwchem", charge=1, mult=2,
                 read_orbitals="1", deltascf={})
    target = _Target([Excitation("alpha", 0, 0), Excitation("beta", 0, 0)])
    calc = _NWChemEngine(p, conf, 1).build(tmp_path / "try1", guess=tmp_path / "ref.movecs", target=target,
                                           attempt=ds.NWCHEM_LADDER[2], n_occ=(8, 7))
    dft = calc.parameters["dft"]
    assert dft["vectors"] == (f'input "{(tmp_path / "ref.movecs").resolve()}" swap alpha 8 9 swap beta 7 8 '
                              f'output "{(tmp_path / "try1").resolve() / "deltascf.movecs"}"')
    assert "odft" in dft and dft["mult"] == 2 and "max_ovl" in dft and dft["iterations"] == 300
    assert dft["convergence"] == {"lshift": 0.0, "damp": 70, "ncydp": 50}
    assert calc.parameters["charge"] == 1


@patch("ensemble_analyzer.calculators.nwchem.NWCHEM_COMMAND", "nwchem")
def test_nwchem_quartet_keeps_occupation_with_mom(conf, tmp_path):
    # without MOM the quartet drifts to a higher solution (validated on the formaldehyde cation)
    p = Protocol(number=2, functional="PBE", basis="def2-SVP", calculator="nwchem", charge=1, mult=2, deltascf={})
    calc = _NWChemEngine(p, conf, 1).build(tmp_path / "q", mult=4, guess=tmp_path / "ref.movecs",
                                           target=_Target([]), attempt=ds.NWCHEM_LADDER[0], n_occ=(9, 6))
    dft = calc.parameters["dft"]
    assert "swap" not in dft["vectors"] and "max_ovl" in dft and dft["mult"] == 4


def _ref(states, n=(8, 8)):
    return _Reference(Path("r.out"), Path("r.gbw"), -114.0, *n, states=states)


def test_target_auto_and_explicit():
    states = [ds.ExcitedState(1, 3.9, [(7, 8, "alpha", 0.99)])]
    settings = {**ds.DEFAULTS, "excitation": "auto"}
    target, exc, kind = DeltaSCFRunner._target(settings, _ref(states), 1, SimpleNamespace(index_base=0))
    assert exc.label == "HOMO->LUMO (alpha)" and kind == "singlet_ap" and target.excitations == [exc]

    settings = {**ds.DEFAULTS, "excitation": "HOMO-1->LUMO", "spin": "alpha"}
    _, exc, kind = DeltaSCFRunner._target(settings, _ref([], (8, 7)), 2, SimpleNamespace(index_base=0))
    assert kind == "doublet_recoupling" and (exc.n, exc.m) == (1, 0)

    settings = {**ds.DEFAULTS, "alphaconf": "0,1"}
    target, exc, kind = DeltaSCFRunner._target(settings, _ref([]), 1, SimpleNamespace(index_base=0))
    assert exc is None and target.alphaconf == "0,1" and kind == "singlet_ap"


# ---------------------------------------------------------------------------
# "inherit"
# ---------------------------------------------------------------------------

def test_excitation_file_roundtrip_and_inherited_target(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    conf = SimpleNamespace(folder="conf_1", number=1)
    a = Protocol(number=0, functional="PBE", basis="def2-SVP", opt=True, deltascf={})
    base = Path("conf_1/protocol_0/deltascf")
    DeltaSCFRunner._save_excitation(base, Excitation("beta", 1, 0, weight=0.9, source="auto"), 3.2, a, -114.1)
    data = DeltaSCFRunner._load_inherited(conf, 0)
    assert (data["spin"], data["n"], data["m"], data["opt"], data["step"]) == ("beta", 1, 0, True, 0)

    settings = {**ds.DEFAULTS, "excitation": "inherit"}
    target, exc, kind = DeltaSCFRunner._target(settings, _ref([], (8, 7)), 2, SimpleNamespace(index_base=0), data)
    assert exc == Excitation("beta", 1, 0, source="inherit") and target.excitations == [exc]
    assert kind == "pure"            # beta into the SOMO of the doublet


def test_inherit_missing_source_or_raw(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    conf = SimpleNamespace(folder="conf_1", number=1)
    with pytest.raises(RuntimeError, match="no accepted ΔSCF state"):
        DeltaSCFRunner._load_inherited(conf, 3)
    DeltaSCFRunner._save_excitation(Path("conf_1/protocol_3/deltascf"), None, 3.0,
                                    SimpleNamespace(opt=False, number=3), -1.0)
    with pytest.raises(RuntimeError, match="raw alphaconf"):
        DeltaSCFRunner._load_inherited(conf, 3)


def test_inherited_gap_comparison():
    runner = DeltaSCFRunner(SimpleNamespace(logger=SimpleNamespace(info=lambda m: None, warning=_Collect())))
    settings = {**ds.DEFAULTS, "inherit_gap_tol_ev": 1.0}
    runner._compare_inherited_gap(5.0, {"gap_ev": 3.5, "step": 1, "opt": False, "label": "HOMO->LUMO"}, settings, "t")
    assert len(runner.log.warning.calls) == 1
    runner._compare_inherited_gap(5.0, {"gap_ev": 3.5, "step": 0, "opt": True, "label": "HOMO->LUMO"}, settings, "t")
    runner._compare_inherited_gap(3.9, {"gap_ev": 3.5, "step": 1, "opt": False, "label": "HOMO->LUMO"}, settings, "t")
    assert len(runner.log.warning.calls) == 1


class _Collect:
    def __init__(self):
        self.calls = []

    def __call__(self, msg):
        self.calls.append(msg)


def test_inherit_reference_has_no_tddft(conf):
    p = Protocol(number=2, functional="PBE0", basis="def2-SVP", deltascf={"excitation": "inherit"})
    _, blocks, _ = _orca_build(_orca(p, conf), mode="sp", nroots=None)
    assert "%tddft" not in blocks
