"""
Tests for Spectral Comparison and Plotting.
Covers ComparedGraph logic and visualization output.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock, patch, mock_open
from ensemble_analyzer._spectral.compare import ComparedGraph

class TestSpectralCompare:
    
    @pytest.fixture
    def setup_compare(self, mock_logger):
        """Setup ComparedGraph instance with mock data."""
        # Mock internal loaders to avoid file system access during init
        with patch("ensemble_analyzer._spectral.compare.ComparedGraph._load_experimental", return_value=(None, None, None, None)):
            with patch("ensemble_analyzer._spectral.compare.ComparedGraph._load_computed", return_value={}):
                
                comp = ComparedGraph(
                    graph_type="IR",
                    log=mock_logger,
                    experimental_file="exp.dat"
                )
                return comp

    def test_initialization_validation(self):
        """Test graph type validation during init."""
        with patch("ensemble_analyzer._spectral.compare.ComparedGraph._load_experimental", return_value=(None, None, None, None)):
            with patch("ensemble_analyzer._spectral.compare.ComparedGraph._load_computed", return_value={}):
                # Valid
                g = ComparedGraph(graph_type="UV")
                assert g.graph_type == "UV"
                
                # Invalid
                with pytest.raises(ValueError):
                    ComparedGraph(graph_type="INVALID_TYPE")

    @patch("ensemble_analyzer._spectral.compare.Path")
    @patch("numpy.loadtxt")
    def test_load_computed(self, mock_loadtxt, mock_path, setup_compare):
        """Test loading of computed spectra files."""
        # Setup mock file system finding
        mock_file = MagicMock()
        mock_file.name = "IR_p1_comp.xy" # Must match pattern {TYPE}_p
        mock_file.stem = "IR_p1_comp"
        # glob returns iterator
        mock_path.return_value.glob.return_value = [mock_file]
        
        # Setup mock file content (X, Y)
        mock_loadtxt.return_value = (np.array([1000., 1100.]), np.array([0.5, 1.0]))
        
        # Execute private method
        data = setup_compare._load_computed()
        
        # Verify
        assert "1" in data
        X, Y = data["1"]
        assert np.allclose(X, [1000., 1100.])
        assert np.max(Y) == 1.0 # Checks normalization if max > 1 (mocked max is 1.0)

    @patch("ensemble_analyzer._spectral.compare.plt")
    def test_plot_generation(self, mock_plt, setup_compare):
        """Test plotting logic."""
        # Inject data manually
        setup_compare.data = {"1": (np.linspace(100, 200, 10), np.random.rand(10))}
        
        # Inject experimental data
        setup_compare.Xr = np.linspace(100, 200, 10)
        setup_compare.Yr = np.random.rand(10)
        setup_compare.bounders = np.array([0, 9]) # indices
        
        # Mock subplots
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)
        
        # Execute
        # We patch pickle dump because _save_or_show tries to save a pickle of the figure
        # We also patch open to prevent creating real files
        with patch("pickle.dump"):
            with patch("builtins.open", mock_open()):
                setup_compare.plot(save=True, show=False)
        
        # Verify
        # The source code uses plt.savefig(), not fig.savefig()
        assert mock_plt.savefig.called
        assert mock_plt.close.called
        # Check title set
        mock_ax.set_title.assert_called()