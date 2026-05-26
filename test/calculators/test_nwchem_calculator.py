import pytest
from unittest.mock import MagicMock, patch
from ase.calculators.nwchem import NWChem
from ensemble_analyzer.calculators.nwchem import NWChemCalc


class TestNWChemCalc:

    @pytest.fixture
    def setup_calc(self, mock_conformer, mock_protocol):
        mock_conformer.folder = "conf_1"
        mock_protocol.charge = 0
        mock_protocol.mult = 1
        mock_protocol.solvent = MagicMock()
        mock_protocol.solvent.solvent = "Water"
        mock_protocol.solvent.smd = False
        mock_protocol.add_input = ""
        mock_protocol.read_orbitals = False
        mock_protocol.constrains = []
        mock_protocol.freq = False
        return mock_conformer, mock_protocol

    def test_common_str_basic(self, setup_calc):
        conf, proto = setup_calc
        calc = NWChemCalc(proto, 4, conf)
        kw = calc.common_str()
        assert kw["theory"] == "dft"
        assert kw["xc"] == "B3LYP"
        assert kw["basis"] == "6-31G*"
        assert "charge" not in kw
        assert kw["dft"]["mult"] == 1

    def test_common_str_charge_nonzero(self, setup_calc):
        conf, proto = setup_calc
        proto.charge = 2
        calc = NWChemCalc(proto, 4, conf)
        kw = calc.common_str()
        assert kw["charge"] == 2
        assert kw["dft"]["mult"] == 1

    def test_common_str_with_solvent(self, setup_calc):
        conf, proto = setup_calc
        calc = NWChemCalc(proto, 4, conf)
        kw = calc.common_str()
        assert "cosmo" in kw
        assert kw["cosmo"]["solvent"] == "water"

    def test_common_str_no_solvent(self, setup_calc):
        conf, proto = setup_calc
        proto.solvent = None
        calc = NWChemCalc(proto, 4, conf)
        kw = calc.common_str()
        assert "cosmo" not in kw

    def test_common_memory(self, setup_calc):
        conf, proto = setup_calc
        calc = NWChemCalc(proto, 8, conf)
        kw = calc.common_str()
        assert kw["memory"] == "40000 mb"

    @patch("ensemble_analyzer.calculators.nwchem.NWCHEM_COMMAND", "nwchem")
    def test_single_point(self, setup_calc):
        conf, proto = setup_calc
        calc = NWChemCalc(proto, 4, conf)
        ase_calc, label = calc.single_point()
        assert label == "nwchem"

    @patch("ensemble_analyzer.calculators.nwchem.NWCHEM_COMMAND", "nwchem")
    def test_optimisation(self, setup_calc):
        conf, proto = setup_calc
        proto.freq = False
        calc = NWChemCalc(proto, 4, conf)
        ase_calc, label = calc.optimisation()
        assert ase_calc.parameters["task"] == "optimize"

    @patch("ensemble_analyzer.calculators.nwchem.NWCHEM_COMMAND", "nwchem")
    def test_optimisation_with_freq(self, setup_calc):
        conf, proto = setup_calc
        proto.freq = True
        proto.opt = True
        calc = NWChemCalc(proto, 4, conf)
        ase_calc, label = calc.optimisation()
        assert 'optimize' in ase_calc.parameters["task"].split()
        assert 'freq' in ase_calc.parameters["task"].split()

    def test_frequency(self, setup_calc):
        conf, proto = setup_calc
        calc = NWChemCalc(proto, 4, conf)
        ase_calc, label = calc.frequency()
        assert ase_calc.parameters["task"] == "freq"

    @patch("ensemble_analyzer.calculators.nwchem.NWCHEM_COMMAND", "nwchem")
    def test_single_point_with_add_input(self, setup_calc):
        conf, proto = setup_calc
        proto.add_input = "scf\n  thresh 1e-8\nend"
        calc = NWChemCalc(proto, 4, conf)
        ase_calc, label = calc.single_point()
        assert label == "nwchem"
        assert ase_calc.write_input.__name__ != NWChem.write_input.__name__
