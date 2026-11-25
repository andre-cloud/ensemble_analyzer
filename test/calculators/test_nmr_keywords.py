"""NMR input generation for ORCA, Gaussian and NWChem (syntax checked on real runs, see plan phase 0)."""

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest
from ase import Atoms

from ensemble_analyzer._nmr.inputs import selected_atoms
from ensemble_analyzer.calculators.gaussian import GaussianCalc
from ensemble_analyzer.calculators.nwchem import NWChemCalc
from ensemble_analyzer.calculators.orca import OrcaCalc
from ensemble_analyzer.protocol.protocol import Protocol

FLUOROETHANE = ["C", "C", "F", "H", "H", "H", "H", "H"]
TMS = {"name": "TMS", "nuclei": ["1H", "13C"]}
CFCL3 = {"name": "CFCl3", "nuclei": ["19F"]}


@pytest.fixture
def conf():
    return SimpleNamespace(folder="conf_1", number=1, atoms=FLUOROETHANE)


def _p(calculator="orca", **nmr):
    settings = {"shielding": True, "coupling": True, "nuclei": ["1H", "13C", "19F"], "references": [TMS, CFCL3]}
    settings.update(nmr)
    if not settings["shielding"]:
        settings["references"] = []
    return Protocol(number=2, functional="B3LYP", basis="def2-SVP", calculator=calculator, nmr=settings)


def _orca_input(protocol, conf, tmp_path) -> str:
    with patch("enan_calculators._orca.shutil.which", return_value="orca"), \
         patch("ensemble_analyzer.calculators.orca.OrcaCalc.VERSION", 6):
        calc, _ = OrcaCalc(protocol, 4, conf).single_point()
    calc.directory = tmp_path
    Path(tmp_path).mkdir(parents=True, exist_ok=True)
    calc.write_inputfiles(Atoms(conf.atoms, positions=[[i, 0, 0] for i in range(len(conf.atoms))]), ["energy"])
    return (Path(tmp_path) / calc.template.inputname).read_text()


# --- atom selection ---------------------------------------------------------

def test_selected_atoms_all_and_subset():
    assert selected_atoms(_p().nmr, FLUOROETHANE) == {"1H": [3, 4, 5, 6, 7], "13C": [0, 1], "19F": [2]}
    assert selected_atoms(_p(atoms=[1, 2, 6, 7]).nmr, FLUOROETHANE) == {"1H": [6, 7], "13C": [1], "19F": [2]}


def test_selected_atoms_errors():
    with pytest.raises(ValueError, match="out of range"):
        selected_atoms(_p(atoms=[8]).nmr, FLUOROETHANE)
    with pytest.raises(ValueError, match="no atom"):
        selected_atoms(_p(nuclei=["31P"], references=[{"name": "H3PO4", "nuclei": ["31P"]}]).nmr, FLUOROETHANE)


# --- ORCA -------------------------------------------------------------------

def test_orca_all_atoms(conf, tmp_path):
    text = _orca_input(_p(), conf, tmp_path)
    assert "! NMR" in text
    assert text.index("%eprnmr") > text.rindex("*")          # after the coordinates
    block = text[text.index("%eprnmr"):]
    assert "Nuclei = all H {shift, ssall, ist = 1}" in block
    assert "Nuclei = all C {shift, ssall, ist = 13}" in block
    assert "Nuclei = all F {shift, ssall, ist = 19}" in block


def test_orca_subset_is_one_based_with_isotopes(conf, tmp_path):
    block = _orca_input(_p(atoms=[1, 2, 6, 7]), conf, tmp_path).split("%eprnmr")[1]
    assert "Nuclei = 7, 8 {shift, ssall, ist = 1}" in block
    assert "Nuclei = 2 {shift, ssall, ist = 13}" in block
    assert "Nuclei = 3 {shift, ssall, ist = 19}" in block


def test_orca_shielding_only_and_fc_only(conf, tmp_path):
    assert "Nuclei = all H {shift, ist = 1}" in _orca_input(_p(coupling=False), conf, tmp_path / "a")
    text = _orca_input(_p(shielding=False, coupling_terms="fc"), conf, tmp_path / "b")
    assert "Nuclei = all H {ssfc, ist = 1}" in text


def test_orca_keeps_user_post_blocks(conf, tmp_path):
    p = _p()
    p.add_input = "%frag 1 {0:1} end end"
    text = _orca_input(p, conf, tmp_path)
    assert "%frag" in text and "%eprnmr" in text


# --- Gaussian ---------------------------------------------------------------

def test_gaussian_route_and_readatoms(conf):
    calc, _ = GaussianCalc(_p("gaussian"), 4, conf).single_point()
    assert calc.parameters["extra"].endswith(" NMR=(GIAO,SpinSpin)")
    assert not calc.parameters.get("addsec")
    calc, _ = GaussianCalc(_p("gaussian", atoms=[1, 2, 6, 7]), 4, conf).single_point()
    assert calc.parameters["extra"].endswith(" NMR=(GIAO,SpinSpin,ReadAtoms)")
    assert calc.parameters["addsec"] == "atoms=2,3,7,8"
    calc, _ = GaussianCalc(_p("gaussian", coupling=False), 4, conf).single_point()
    assert calc.parameters["extra"].endswith(" NMR=GIAO")


def test_gaussian_readatoms_when_nuclei_are_a_subset(conf):
    p = _p("gaussian", nuclei=["1H"], references=[TMS | {"nuclei": ["1H"]}])
    calc, _ = GaussianCalc(p, 4, conf).single_point()
    assert calc.parameters["addsec"] == "atoms=4,5,6,7,8"


# --- NWChem -----------------------------------------------------------------

@patch("ensemble_analyzer.calculators.nwchem.NWCHEM_COMMAND", "nwchem")
def test_nwchem_property_block(conf):
    calc, _ = NWChemCalc(_p("nwchem"), 4, conf).single_point()
    assert calc.parameters["task"] == "property"
    assert calc.parameters["property"] == {"shielding": None, "spinspin": None}
    assert "odft" in calc.parameters["dft"]


@patch("ensemble_analyzer.calculators.nwchem.NWCHEM_COMMAND", "nwchem")
def test_nwchem_subset_and_shielding_only(conf):
    calc, _ = NWChemCalc(_p("nwchem", atoms=[1, 2, 6, 7]), 4, conf).single_point()
    assert calc.parameters["property"] == {"shielding": "4 2 3 7 8", "spinspin": "6 2 3 2 7 2 8 3 7 3 8 7 8"}
    calc, _ = NWChemCalc(_p("nwchem", coupling=False), 4, conf).single_point()
    assert calc.parameters["property"] == {"shielding": None}
    assert "odft" not in calc.parameters["dft"]
