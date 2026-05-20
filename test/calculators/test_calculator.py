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
        proto.constrains = [[1], [2]] # Atom indices
        
        calc = GaussianCalc(proto, 4, conf)
        ase_calc, label = calc.optimisation()
        
        assert "opt=(modredudant)" in ase_calc.parameters["extra"]
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
            si, ob, post = calc.common_str()
            
            assert "B3LYP 6-31G*" in si
            assert "CPCM" in si 
            assert "nopop" in si
            assert "%pal nprocs 4 end" in ob
            assert post == ""

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
        proto.constrains = [[1]]
        with patch("ensemble_analyzer._calculators._orca.orca_profile"):
            calc = OrcaCalc(proto, 4, conf)
            ase_calc, label = calc.optimisation()
        
            assert "%geom Constraints {C 1 C} end end" in ase_calc.parameters["orcasimpleinput"]

    def test_constraints_all_types(self, setup_calc):
        conf, proto = setup_calc
        proto.constrains = [[1, 2], [1, 2, 3], [1, 2, 3, 4], [1]]

        with patch("ensemble_analyzer._calculators._orca.orca_profile"):
            calc = OrcaCalc(proto, 4, conf)
            ase_calc, label = calc.optimisation()
            assert "%geom Constraints {B 1 2 C} {A 1 2 3 C} {D 1 2 3 4 C} {C 1 C} end end" in ase_calc.parameters["orcasimpleinput"]

        calc = GaussianCalc(proto, 4, conf)
        ase_calc, label = calc.optimisation()
        assert "B 1 2 F" in ase_calc.parameters["addsec"]
        assert "A 1 2 3 F" in ase_calc.parameters["addsec"]
        assert "D 1 2 3 4 F" in ase_calc.parameters["addsec"]
        assert "X 1 F" in ase_calc.parameters["addsec"]


class TestSplitPostBlocks:
    """Tests for _split_post_blocks edge cases (indented vs compact)."""

    def test_indented_frag(self):
        from ensemble_analyzer._calculators._orca import _split_post_blocks
        t = "%frag\n Definition\n  1 {18:27} end\n  2 {0:17} end\n end\nend"
        pre, post = _split_post_blocks(t)
        assert pre == ""
        assert post == "%frag\n Definition\n  1 {18:27} end\n  2 {0:17} end\n end\nend"

    def test_compact_frag(self):
        from ensemble_analyzer._calculators._orca import _split_post_blocks
        t = "%frag\nDefinition\n1 {18:27} end\n2 {0:17} end\nend\nend"
        pre, post = _split_post_blocks(t)
        assert pre == ""
        assert post == "%frag\nDefinition\n1 {18:27} end\n2 {0:17} end\nend\nend"

    def test_newline_prefix_frag(self):
        from ensemble_analyzer._calculators._orca import _split_post_blocks
        t = "\n%frag\nDefinition\n1 {18:27} end\n2 {0:17} end\nend\nend"
        pre, post = _split_post_blocks(t)
        assert post == "%frag\nDefinition\n1 {18:27} end\n2 {0:17} end\nend\nend"

    def test_mixed_pre_and_post(self):
        from ensemble_analyzer._calculators._orca import _split_post_blocks
        t = "%scf maxiter 500 end\n%frag\nDefinition\n1 {18:27} end\nend\nend"
        pre, post = _split_post_blocks(t)
        assert pre == "%scf maxiter 500 end"
        assert post == "%frag\nDefinition\n1 {18:27} end\nend\nend"

    def test_multiple_post_blocks(self):
        from ensemble_analyzer._calculators._orca import _split_post_blocks
        t = "%frag\n1 {0:5} end\nend\n%eprnmr\ngtensor 1\nend"
        pre, post = _split_post_blocks(t)
        assert pre == ""
        assert post == "%frag\n1 {0:5} end\nend\n%eprnmr\ngtensor 1\nend"

    def test_pre_non_post_block_and_post(self):
        from ensemble_analyzer._calculators._orca import _split_post_blocks
        t = "%output\n print[ P_Mulliken 1 ] end\n%frag\n1 {0:5} end\nend"
        pre, post = _split_post_blocks(t)
        assert pre == "%output\n print[ P_Mulliken 1 ] end"
        assert post == "%frag\n1 {0:5} end\nend"

    def test_empty(self):
        from ensemble_analyzer._calculators._orca import _split_post_blocks
        pre, post = _split_post_blocks("")
        assert pre == ""
        assert post == ""

    def test_no_post_blocks(self):
        from ensemble_analyzer._calculators._orca import _split_post_blocks
        t = "%scf maxiter 500 end\n%maxcore 8000"
        pre, post = _split_post_blocks(t)
        assert pre == t
        assert post == ""