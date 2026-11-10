"""
NWChem TD-DFT input generation and MO-vector restart (``read_orbitals``).

``nroots`` generates ``task tddft energy`` and ``read_orbitals`` reads the MO
vectors written by the referenced protocol step.
"""

from unittest.mock import patch

import pytest

from ensemble_analyzer.calculators.nwchem import NWChemCalc


@pytest.fixture
def tddft_protocol(mock_protocol):
    mock_protocol.charge = 0
    mock_protocol.mult = 1
    mock_protocol.solvent = None
    mock_protocol.add_input = ""
    mock_protocol.read_orbitals = ""
    mock_protocol.constraints = []
    mock_protocol.ts = False
    mock_protocol.opt = False
    mock_protocol.freq = False
    mock_protocol.nroots = 10
    mock_protocol.tda = True
    return mock_protocol


@pytest.fixture
def conf(mock_conformer):
    mock_conformer.folder = "conf_1"
    mock_conformer.number = 1
    return mock_conformer



@patch("ensemble_analyzer.calculators.nwchem.NWCHEM_COMMAND", "nwchem")
class TestNWChemTDDFT:

    def test_single_point_tddft_task(self, tddft_protocol, conf):
        tddft_protocol.calculator = "nwchem"
        ase_calc, _ = NWChemCalc(tddft_protocol, 4, conf).single_point()
        assert ase_calc.parameters["theory"] == "tddft"
        assert ase_calc.parameters["task"] == "energy"
        assert ase_calc.parameters["tddft"] == {"nroots": 10, "cis": None}

    def test_rpa_when_not_tda(self, tddft_protocol, conf):
        tddft_protocol.tda = False
        ase_calc, _ = NWChemCalc(tddft_protocol, 4, conf).single_point()
        assert ase_calc.parameters["tddft"] == {"nroots": 10, "rpa": None}

    def test_no_tddft_without_nroots(self, tddft_protocol, conf):
        tddft_protocol.nroots = None
        ase_calc, _ = NWChemCalc(tddft_protocol, 4, conf).single_point()
        assert ase_calc.parameters.get("theory", "dft") == "dft"
        assert "tddft" not in ase_calc.parameters

    def test_optimisation_task_untouched(self, tddft_protocol, conf):
        tddft_protocol.opt = True
        ase_calc, _ = NWChemCalc(tddft_protocol, 4, conf).optimisation()
        assert ase_calc.parameters["task"] == "optimize"
        assert "tddft" not in ase_calc.parameters


@patch("ensemble_analyzer.calculators.nwchem.NWCHEM_COMMAND", "nwchem")
class TestNWChemReadOrbitals:

    def test_vectors_output_only(self, tddft_protocol, conf):
        tddft_protocol.number = 3
        ase_calc, _ = NWChemCalc(tddft_protocol, 4, conf).single_point()
        vectors = ase_calc.parameters["dft"]["vectors"]
        assert vectors.startswith('output "')
        assert vectors.endswith('conf_1/protocol_3/1_p3_nwchem/1_p3_nwchem.movecs"')

    def test_vectors_read_from_referenced_step(self, tddft_protocol, conf):
        tddft_protocol.number = 3
        tddft_protocol.read_orbitals = "1"
        ase_calc, _ = NWChemCalc(tddft_protocol, 4, conf).single_point()
        vectors = ase_calc.parameters["dft"]["vectors"]
        assert 'conf_1/protocol_1/1_p1_nwchem/1_p1_nwchem.movecs" output "' in vectors
        assert vectors.endswith('conf_1/protocol_3/1_p3_nwchem/1_p3_nwchem.movecs"')

    @pytest.mark.parametrize("step", ["0", 0])
    def test_step_zero_is_a_valid_reference(self, tddft_protocol, conf, step):
        tddft_protocol.number = 2
        tddft_protocol.read_orbitals = step
        ase_calc, _ = NWChemCalc(tddft_protocol, 4, conf).single_point()
        assert "protocol_0/1_p0_nwchem/1_p0_nwchem.movecs" in ase_calc.parameters["dft"]["vectors"]

    def test_common_str_uses_same_paths(self, tddft_protocol, conf):
        tddft_protocol.number = 3
        tddft_protocol.read_orbitals = "1"
        calc = NWChemCalc(tddft_protocol, 4, conf)
        assert calc.common_str()["dft"]["vectors"] == calc._vectors_directive()
