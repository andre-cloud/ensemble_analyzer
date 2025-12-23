import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from ensemble_analyzer._spectral.comp_electronic import ComputedElectronic
from ensemble_analyzer._spectral.comp_vibronic import ComputedVibronic
from ensemble_analyzer._spectral.base import BaseGraph

class TestSpectral:
    
    @pytest.fixture
    def mock_defaults(self):
        with patch("ensemble_analyzer._spectral.base.GraphDefault") as MockDef:
            MockDef.return_value.start = 0
            MockDef.return_value.end = 100
            yield MockDef

    @pytest.fixture
    def base_setup(self, mock_defaults, mock_protocol, mock_logger):
        return ComputedElectronic(
            confs=[],
            protocol=mock_protocol,
            graph_type="UV",
            log=mock_logger,
            definition=2 
        )

    def test_retrieve_data_filtering(self, base_setup, mock_conformer, mock_protocol):
        # Create a special mock class to handle the custom dunder method
        class MockStore:
            def __getitem__(self, *args, **kwargs):
                return MagicMock(X=[10, 20], Y=[0.1, 0.2])
            def __contains__(self, item):
                return True
            def __has_graph_type__(self, *args):
                return True

        c1 = mock_conformer
        c1.graphs_data = MockStore() 
        c1.energies.__getitem__.return_value.Pop = 0.5
        
        c2 = MagicMock(active=False)

        c3 = MagicMock(active=True)
        c3.graphs_data.__contains__.return_value = False
        
        base_setup.confs = [c1, c2, c3]
        
        base_setup.retrieve_data(mock_protocol)
        
        assert len(base_setup.energies) == 2
        assert np.allclose(base_setup.impulse, [0.05, 0.1])

    def test_normalize(self, base_setup):
        y = np.array([10.0, 5.0, -20.0])
        norm = base_setup.normalize(y)
        assert np.allclose(norm, [0.5, 0.25, -1.0])

    def test_electronic_convolution_shift(self, mock_defaults, mock_protocol, mock_logger):
        spec = ComputedElectronic(confs=[], protocol=mock_protocol, graph_type="UV", log=mock_logger)
        energies = np.array([50.0])
        impulses = np.array([1.0])
        shift = 10.0
        fwhm = 1.0
        
        with patch.object(spec, 'gaussian') as mock_gauss:
            spec.convolute(energies, impulses, shift, fwhm)
            passed_energies = mock_gauss.call_args[0][0]
            assert passed_energies[0] == 60.0 

    def test_vibronic_convolution_shift(self, mock_defaults, mock_protocol, mock_logger):
        spec = ComputedVibronic(confs=[], protocol=mock_protocol, graph_type="IR", log=mock_logger)
        energies = np.array([1000.0])
        impulses = np.array([1.0])
        shift = 0.98 
        fwhm = 4.0
        
        with patch.object(spec, 'lorentzian') as mock_lor:
            spec.convolute(energies, impulses, shift, fwhm)
            passed_energies = mock_lor.call_args[0][0]
            assert passed_energies[0] == 980.0 

    def test_autoconvolution_logic(self, base_setup):
        base_setup.ref = MagicMock()
        base_setup.ref.Y = np.zeros(100)
        base_setup.ref.x_min_idx = 0
        base_setup.ref.x_max_idx = 100
        # MUST use real numpy array for Numba @njit compatibility
        base_setup.ref.weight = np.ones(100)
        
        base_setup.energies = np.array([10.0])
        base_setup.impulse = np.array([1.0])
        
        # Boundaries must be set before optimization
        base_setup.set_boundaries()
        
        with patch("ensemble_analyzer._spectral.base.minimize") as mock_min:
            mock_min.return_value.success = True
            mock_min.return_value.x = [0.5, 0.5] 
            
            base_setup.autoconvolution()
            
            assert base_setup.SHIFT == 0.5
            assert base_setup.FWHM == 0.5
            base_setup.log.spectra_result.assert_called()