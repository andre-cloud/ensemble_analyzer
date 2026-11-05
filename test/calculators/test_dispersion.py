"""
Tests for the standalone `dispersion` protocol keyword.
Verifies normalisation/validation in Protocol and the engine-specific
translation in ORCA, Gaussian and NWChem input generation.
"""

import pytest
from unittest.mock import patch
from ase.io.nwchem.nwwriter import _format_block

from ensemble_analyzer.protocol.protocol import Protocol
from ensemble_analyzer.protocol.dispersion import normalize_dispersion
from ensemble_analyzer.calculators.gaussian import GaussianCalc
from ensemble_analyzer.calculators.nwchem import NWChemCalc
from ensemble_analyzer.calculators.orca import OrcaCalc


class TestProtocolDispersion:

    @pytest.mark.parametrize("value, expected", [
        (None, None), ("", None), ("none", None),
        ("D3BJ", "D3BJ"), ("d3(bj)", "D3BJ"), ("GD3BJ", "D3BJ"),
        ("d3zero", "D3ZERO"), ("GD3", "D3ZERO"), ("d4", "D4"), ("pfd", "PFD"),
    ])
    def test_normalize(self, value, expected):
        assert normalize_dispersion(value) == expected

    def test_ambiguous_d3_rejected(self):
        with pytest.raises(ValueError, match="ambiguous"):
            Protocol(number=0, functional="B3LYP", basis="def2-SVP", dispersion="D3")

    def test_unknown_rejected(self):
        with pytest.raises(ValueError, match="not recognised"):
            Protocol(number=0, functional="B3LYP", basis="def2-SVP", dispersion="D5")

    @pytest.mark.parametrize("calculator, disp", [
        ("gaussian", "D4"), ("nwchem", "D4"), ("nwchem", "PFD"), ("orca", "PFD"),
    ])
    def test_unsupported_by_engine(self, calculator, disp):
        with pytest.raises(ValueError, match="not available"):
            Protocol(number=0, functional="B3LYP", basis="def2-SVP",
                     calculator=calculator, dispersion=disp)

    def test_ml_calculator_rejected(self):
        with pytest.raises(ValueError, match="does not support"):
            Protocol(number=0, functional="model.pt", calculator="mace", dispersion="D3BJ")

    @pytest.mark.parametrize("functional", ["wB97X-D4", "B3LYP-D3BJ", "wB97XD", "B97D3", "r2SCAN-3c"])
    def test_builtin_dispersion_rejected(self, functional):
        with pytest.raises(ValueError, match="already includes"):
            Protocol(number=0, functional=functional, basis="def2-SVP", dispersion="D3BJ")

    def test_double_request_in_add_input_rejected(self):
        with pytest.raises(ValueError, match="add_input"):
            Protocol(number=0, functional="B3LYP", basis="def2-SVP", calculator="gaussian",
                     dispersion="D3BJ", add_input="EmpiricalDispersion=GD3BJ")

    def test_repr(self):
        p = Protocol(number=0, functional="B3LYP", basis="def2-SVP", dispersion="d3bj")
        assert repr(p) == "B3LYP-D3BJ/def2-SVP"


class TestCalculatorDispersion:

    @pytest.fixture
    def setup_calc(self, mock_conformer, mock_protocol):
        mock_conformer.folder = "conf_1"
        mock_protocol.charge = 0
        mock_protocol.mult = 1
        mock_protocol.solvent = None
        mock_protocol.add_input = ""
        mock_protocol.read_orbitals = False
        mock_protocol.constraints = []
        mock_protocol.opt = False
        mock_protocol.freq = False
        mock_protocol.ts = False
        mock_protocol.nroots = None
        mock_protocol.dispersion = "D3BJ"
        return mock_conformer, mock_protocol

    def test_gaussian_route(self, setup_calc):
        conf, proto = setup_calc
        calc = GaussianCalc(proto, 4, conf)
        assert "# B3LYP/6-31G* EmpiricalDispersion=GD3BJ" in calc.common_str()
        ase_calc, _ = calc.single_point()
        assert "EmpiricalDispersion=GD3BJ" in ase_calc.parameters["extra"]

    def test_gaussian_d3zero(self, setup_calc):
        conf, proto = setup_calc
        proto.dispersion = "D3ZERO"
        ase_calc, _ = GaussianCalc(proto, 4, conf).optimisation()
        assert "EmpiricalDispersion=GD3 " in ase_calc.parameters["extra"]

    def test_orca_simple_input(self, setup_calc):
        conf, proto = setup_calc
        proto.dispersion = "D4"
        with patch("enan_calculators._orca.shutil.which", return_value="orca"), \
             patch("ensemble_analyzer.calculators.orca.orca_profile"):
            calc = OrcaCalc(proto, 4, conf)
            si, _, _ = calc.common_str()
            assert si.split()[-1] == "D4"
            ase_calc, _ = calc.single_point()
            assert " D4" in ase_calc.parameters["orcasimpleinput"]

    @patch("ensemble_analyzer.calculators.nwchem.NWCHEM_COMMAND", "nwchem")
    def test_nwchem_dft_block(self, setup_calc):
        conf, proto = setup_calc
        calc = NWChemCalc(proto, 4, conf)
        assert calc.common_str()["dft"]["disp"] == {"vdw": 4}
        ase_calc, _ = calc.single_point()
        block = _format_block("dft", ase_calc.parameters["dft"])
        assert "  disp vdw 4" in block
        assert block.index("  disp vdw 4") < block.index("end")

    @patch("ensemble_analyzer.calculators.nwchem.NWCHEM_COMMAND", "nwchem")
    def test_no_dispersion_leaves_input_untouched(self, setup_calc):
        conf, proto = setup_calc
        proto.dispersion = None
        ase_calc, _ = NWChemCalc(proto, 4, conf).single_point()
        assert "disp" not in ase_calc.parameters["dft"]
        g_calc, _ = GaussianCalc(proto, 4, conf).single_point()
        assert "EmpiricalDispersion" not in g_calc.parameters["extra"]
