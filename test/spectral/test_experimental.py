"""
Tests for Experimental Data Loading.
Covers ExperimentalGraph class.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock, patch, mock_open
from ensemble_analyzer._spectral.experimental import ExperimentalGraph

class TestExperimental:
    
    @pytest.fixture
    def exp_graph(self, mock_logger, mock_protocol):
        # ExperimentalGraph inherits from BaseGraph, so it needs:
        # confs, protocol, graph_type, log
        # It DOES NOT accept 'filename' in __init__.
        
        # We pass an empty list of confs since we are testing experimental loading
        graph = ExperimentalGraph(
            confs=[],
            protocol=mock_protocol,
            graph_type="IR",
            log=mock_logger
        )
        # Manually set the filename attribute which is usually set in __post_init__ via defaults
        graph.fname = "spectrum.dat"
        return graph

    def test_load_data_success(self, exp_graph):
        """Test parsing of experimental data file."""
        # Mock file content: X Y pairs
        content = "1000.0  0.5\n1100.0  0.8\n"
        
        with patch("builtins.open", mock_open(read_data=content)):
            # Mock loadtxt since it might check file existence or open it itself
            with patch("numpy.loadtxt") as mock_load:
                mock_load.return_value = np.array([[1000.0, 0.5], [1100.0, 0.8]])
                
                # The method to call is load_file_experimental, not load
                exp_graph.load_file_experimental()
                
                assert exp_graph.x is not None
                assert len(exp_graph.x) == 2
                # Check normalization (max is 0.8 -> becomes 1.0)
                assert exp_graph.y[1] == 1.0 

    def test_load_data_file_not_found(self, exp_graph):
        """Test handling of missing file."""
        # We patch loadtxt to raise OSError/FileNotFoundError
        # Note: experimental.py typically lets the exception propagate or handles it
        # depending on implementation. Looking at your code, it seems to use loadtxt directly.
        # If it crashes, we expect the crash (or need to wrap it).
        
        with patch("numpy.loadtxt", side_effect=FileNotFoundError):
            with pytest.raises(FileNotFoundError):
                exp_graph.load_file_experimental()

    def test_normalization_logic(self, exp_graph):
        """Test normalization of experimental signals."""
        y_raw = np.array([0.5, 1.0, 0.2])
        # BaseGraph provides normalize, inherited by ExperimentalGraph
        y_norm = exp_graph.normalize(y_raw)
        assert np.max(y_norm) == 1.0
        
        # Test negative peaks (e.g. VCD)
        y_raw_vcd = np.array([0.5, -2.0, 0.2])
        y_norm_vcd = exp_graph.normalize(y_raw_vcd)
        # Should normalize by max abs value (2.0)
        assert np.min(y_norm_vcd) == -1.0
        assert y_norm_vcd[0] == 0.25