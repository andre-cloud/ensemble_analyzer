"""
Tests for QM Calculators (Gaussian, ORCA).
Verifies input file generation for different calculation types (OPT, FREQ, SP).
"""

import pytest
from unittest.mock import MagicMock, patch
from ensemble_analyzer.calculators.gaussian import GaussianCalc
from ensemble_analyzer.calculators.orca import OrcaCalc

class TestCalculators:
    
    @pytest.fixture
    def setup_calc(self, mock_conformer, mock_protocol):
        """Shared setup for calculator instantiation."""
        # Setup specific protocol details
        mock_protocol.charge = 0
        mock_protocol.mult = 1
        
        # Setup Solvent Mock
        mock_protocol.solvent = MagicMock()
        mock_protocol.solvent.solvent = "Water"
        mock_protocol.solvent.smd = False # Default CPCM

        mock_protocol.solvent.__str__.return_value = "CPCM"
        
        mock_protocol.add_input = ""
        mock_protocol.read_orbitals = False
        mock_protocol.constraints = []
        mock_protocol.ts = False
        
        return mock_conformer, mock_protocol

    def test_gaussian_opt_constraints(self, setup_calc):
        conf, proto = setup_calc
        proto.constraints = [[1], [2]] # Atom indices
        
        calc = GaussianCalc(proto, 4, conf)
        ase_calc, label = calc.optimisation()
        
        assert "opt=(modredundant)" in ase_calc.parameters["extra"]
        assert "X 1 F" in ase_calc.parameters["addsec"]
        assert "X 2 F" in ase_calc.parameters["addsec"]

    def test_orca_freq_block(self, setup_calc):
        conf, proto = setup_calc
        proto.freq = True
        with patch("enan_calculators._orca.shutil.which", return_value="orca"), \
             patch("ensemble_analyzer.calculators.orca.OrcaCalc.VERSION", 6):
            calc = OrcaCalc(proto, 4, conf)
            ase_calc, label = calc.frequency()

            assert "freq" in ase_calc.parameters["orcasimpleinput"]
            assert "%freq vcd true end" in ase_calc.parameters["orcablocks"]

    def test_orca_constraints(self, setup_calc):
        conf, proto = setup_calc
        proto.constraints = [[1]]
        with patch("enan_calculators._orca.shutil.which", return_value="orca"):
            calc = OrcaCalc(proto, 4, conf)
            ase_calc, label = calc.optimisation()
        
            assert "%geom Constraints {C 1 C} end end" in ase_calc.parameters["orcasimpleinput"]

    def test_constraints_all_types(self, setup_calc):
        conf, proto = setup_calc
        proto.constraints = [[1, 2], [1, 2, 3], [1, 2, 3, 4], [1]]

        with patch("enan_calculators._orca.shutil.which", return_value="orca"):
            calc = OrcaCalc(proto, 4, conf)
            ase_calc, label = calc.optimisation()
            assert "%geom Constraints {B 1 2 C} {A 1 2 3 C} {D 1 2 3 4 C} {C 1 C} end end" in ase_calc.parameters["orcasimpleinput"]

        calc = GaussianCalc(proto, 4, conf)
        ase_calc, label = calc.optimisation()
        assert "B 1 2 F" in ase_calc.parameters["addsec"]
        assert "A 1 2 3 F" in ase_calc.parameters["addsec"]
        assert "D 1 2 3 4 F" in ase_calc.parameters["addsec"]
        assert "X 1 F" in ase_calc.parameters["addsec"]



    @pytest.mark.parametrize("smd, model", [(True, "SMD"), (False, "CPCM")])
    def test_gaussian_solvent_model(self, setup_calc, smd, model):
        from ensemble_analyzer.protocol.solvent import Solvent
        conf, proto = setup_calc
        proto.solvent = Solvent("Water", smd=smd)
        ase_calc, _ = GaussianCalc(proto, 4, conf).single_point()
        assert f"SCRF=({model},Solvent=Water)" in ase_calc.parameters["extra"]
