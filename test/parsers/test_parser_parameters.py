"""
Tests for Parser Parameter extraction.
Covers get_conf_parameters logic, thermodynamics integration, and error handling.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from datetime import datetime
from ensemble_analyzer.parser_parameter import get_conf_parameters

class TestParserParameter:

    @pytest.fixture
    def mock_args(self, mock_conformer, mock_protocol, mock_logger):
        """Setup common arguments for get_conf_parameters."""
        mock_conformer.folder = "conf_1"
        mock_conformer.number = 1
        
        # Setup energies/spectral storage mocks
        mock_conformer.energies.add = MagicMock()
        mock_conformer.graphs_data.add = MagicMock()
        
        # Setup Protocol defaults
        mock_protocol.calculator = "mock_calc"
        mock_protocol.opt = False
        mock_protocol.freq = False
        mock_protocol.add_input = ""
        mock_protocol.functional = "B3LYP"
        mock_protocol.freq_fact = 1.0
        mock_protocol.mult = 1
        mock_protocol.skip_opt_fail = False
        
        return {
            "conf": mock_conformer,
            "number": 1,
            "output": "out.log",
            "p": mock_protocol,
            "time": datetime.now(),
            "temp": 298.15,
            "log": mock_logger
        }

    @patch("ensemble_analyzer.parser_parameter.PARSER_REGISTRY")
    def test_successful_parsing_full(self, mock_registry, mock_args):
        """Test parsing of a complete calculation (OPT+FREQ) with valid results."""
        # Setup Mock Parser
        MockParserClass = MagicMock()
        mock_parser = MockParserClass.return_value
        mock_registry.__getitem__.return_value = MockParserClass
        
        # Configure Parser Returns
        mock_parser.correct_exiting = True
        mock_parser.parse_energy.return_value = -100.0
        mock_parser.opt_done.return_value = True
        mock_parser.parse_geom.return_value = np.zeros((3,3))
        # Freq, IR, VCD
        mock_parser.parse_freq.return_value = (np.array([100.0, 200.0]), np.ones((2,2)), np.ones((2,2)))
        # B, M vectors
        mock_parser.parse_B_m.return_value = (np.array([1,0,0]), np.array([0,1,0]))
        # TD-DFT (UV, ECD)
        mock_parser.parse_tddft.return_value = (np.ones((2,2)), np.ones((2,2)))

        # Enable OPT + FREQ flags
        mock_args["p"].opt = True
        mock_args["p"].freq = True
        
        # Execute
        with patch("ensemble_analyzer.parser_parameter.free_gibbs_energy") as mock_rrho:
            # Mock thermodynamics calculation
            mock_rrho.return_value = (-100.1, 0.1, -100.05, 0.05) # G, ZPVE, H, S
            
            success = get_conf_parameters(**mock_args)
            
        assert success is True
        # Verify storage
        mock_args["conf"].energies.add.assert_called_once()
        # Verify geometry update
        assert np.array_equal(mock_args["conf"].last_geometry, np.zeros((3,3)))

    @patch("ensemble_analyzer.parser_parameter.PARSER_REGISTRY")
    def test_optimization_failure(self, mock_registry, mock_args):
        """Test handling of failed optimization."""
        MockParserClass = MagicMock()
        mock_parser = MockParserClass.return_value
        mock_registry.__getitem__.return_value = MockParserClass
        
        mock_parser.correct_exiting = True
        mock_parser.parse_energy.return_value = -100.0
        mock_parser.opt_done.return_value = False # OPT FAILED
        
        mock_args["p"].opt = True
        
        # Case 1: Strict (Default) -> Should raise Exception -> caught -> return False
        mock_args["p"].skip_opt_fail = False
        success = get_conf_parameters(**mock_args)
        assert success is False
        assert mock_args["conf"].active is False
        
        # Case 2: Permissive -> Should warn and return True but deactivate
        mock_args["conf"].active = True # Reset
        mock_args["p"].skip_opt_fail = True
        success = get_conf_parameters(**mock_args)
        assert success is True
        assert mock_args["conf"].active is False # Still deactivated, but handled gracefully

    @patch("ensemble_analyzer.parser_parameter.PARSER_REGISTRY")
    def test_missing_frequencies_inheritance(self, mock_registry, mock_args):
        """Test fetching thermo data from previous step when frequencies are missing."""
        MockParserClass = MagicMock()
        mock_parser = MockParserClass.return_value
        mock_registry.__getitem__.return_value = MockParserClass
        
        mock_parser.correct_exiting = True
        mock_parser.parse_energy.return_value = -100.0
        mock_parser.parse_freq.return_value = (np.array([]), None, None) # No freq
        mock_parser.parse_B_m.return_value = (None, None)
        
        mock_args["p"].opt = False
        mock_args["p"].freq = False
        
        # Mock previous protocol data
        prev_data = MagicMock()
        prev_data.G_E = -0.1 # Correction from prev step
        prev_data.zpve = 0.05
        prev_data.H = -100.05
        prev_data.S = 0.2
        
        mock_args["conf"].energies.__getitem__.return_value = prev_data
        
        success = get_conf_parameters(**mock_args)
        
        assert success is True
        # Verify that stored G = E + prev_G_E
        # args[0][1] is the EnergyRecord passed to add
        stored_record = mock_args["conf"].energies.add.call_args[0][1]
        assert stored_record.G == -100.1 # -100.0 + (-0.1)

    @patch("ensemble_analyzer.parser_parameter.PARSER_REGISTRY")
    def test_crash_detection(self, mock_registry, mock_args):
        """Test detection of abnormal termination."""
        MockParserClass = MagicMock()
        mock_parser = MockParserClass.return_value
        mock_registry.__getitem__.return_value = MockParserClass
        
        mock_parser.correct_exiting = False # Gaussian crashed
        
        success = get_conf_parameters(**mock_args)
        
        assert success is True # Handled "gracefully" in logic flow
        assert mock_args["conf"].active is False
        mock_args["log"].warning.assert_called_with(
            "Parser detected abnormal termination for Conf 1. Deactivating."
        )