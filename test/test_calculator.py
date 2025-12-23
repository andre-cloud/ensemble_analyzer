"""
Tests for QM Calculators (Gaussian, ORCA).
Verifies input file generation for different calculation types (OPT, FREQ, SP).
"""

import pytest
from unittest.mock import MagicMock, patch
from ensemble_analyzer._calculators._gaussian import GaussianCalc
from ensemble_analyzer._calculators._orca import OrcaCalc

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
        mock_protocol.constrains = []
        
        return mock_conformer, mock_protocol

    def test_gaussian_common_string(self, setup_calc):
        conf, proto = setup_calc
        # Signature: (protocol, cpu, conf)
        calc = GaussianCalc(proto, 4, conf)
        
        route = calc.common_str()
        assert "# B3LYP/6-31G* SCRF=(CPCM,Solvent=Water)" in route

    def test_gaussian_opt_constraints(self, setup_calc):
        conf, proto = setup_calc
        proto.constrains = [0, 1] # Atom indices
        
        calc = GaussianCalc(proto, 4, conf)
        ase_calc, label = calc.optimisation()
        
        assert "opt=(modredudant)" in ase_calc.parameters["extra"]
        # Gaussian indices start at 1
        assert "X 1 F" in ase_calc.parameters["addsec"]
        assert "X 2 F" in ase_calc.parameters["addsec"]

    def test_gaussian_smd_solvent(self, setup_calc):
        conf, proto = setup_calc
        proto.solvent.smd = True
        
        calc = GaussianCalc(proto, 4, conf)
        route = calc.common_str()
        assert "SCRF=(SMD,Solvent=Water)" in route

    def test_orca_common_string(self, setup_calc):
        conf, proto = setup_calc
        # Mock ORCA profile availability
        with patch("ensemble_analyzer._calculators._orca.orca_profile"):
            calc = OrcaCalc(proto, 4, conf)
            si, ob = calc.common_str()
            
            assert "B3LYP 6-31G*" in si
            assert "CPCM" in si 
            assert "nopop" in si
            assert "%pal nprocs 4 end" in ob

    def test_orca_freq_block(self, setup_calc):
        conf, proto = setup_calc
        proto.freq = True
        with patch("ensemble_analyzer._calculators._orca.orca_profile"):
            with patch("ensemble_analyzer._calculators._orca.OrcaCalc.VERSION", 6):
                calc = OrcaCalc(proto, 4, conf)
                ase_calc, label = calc.frequency()

                assert "freq" in ase_calc.parameters["orcasimpleinput"]
                assert "%freq vcd true end" in ase_calc.parameters["orcablocks"]

    def test_orca_constraints(self, setup_calc):
        conf, proto = setup_calc
        proto.constrains = [0]
        with patch("ensemble_analyzer._calculators._orca.orca_profile"):
            calc = OrcaCalc(proto, 4, conf)
            ase_calc, label = calc.optimisation()
        
            assert "%geom Constraints  {C 0 C}end end" in ase_calc.parameters["orcasimpleinput"]