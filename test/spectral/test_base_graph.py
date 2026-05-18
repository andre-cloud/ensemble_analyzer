import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from ensemble_analyzer._spectral.base import BaseGraph
from ensemble_analyzer._spectral.base import gaussian_njit, lorentzian_njit, diversity_function_njit


class TestNumbaFunctions:

    def test_gaussian_njit(self):
        X = np.linspace(0, 100, 1001)
        x0 = np.array([50.0])
        I = np.array([1.0])
        Y = gaussian_njit(X, x0, I, 10.0)
        assert np.isclose(np.max(Y), 1.0 / (10.0 / (2 * np.sqrt(2 * np.log(2))) * np.sqrt(2 * np.pi)), atol=0.1)
        assert Y.shape == X.shape

    def test_gaussian_njit_empty_peaks(self):
        X = np.linspace(0, 100, 100)
        Y = gaussian_njit(X, np.array([]), np.array([]), 10.0)
        assert np.all(Y == 0.0)

    def test_lorentzian_njit(self):
        X = np.linspace(0, 100, 1001)
        x0 = np.array([50.0])
        I = np.array([1.0])
        Y = lorentzian_njit(X, x0, I, 10.0)
        assert np.isclose(Y[500], 1.0)
        assert Y.shape == X.shape

    def test_lorentzian_njit_empty_peaks(self):
        X = np.linspace(0, 100, 100)
        Y = lorentzian_njit(X, np.array([]), np.array([]), 10.0)
        assert np.all(Y == 0.0)

    def test_diversity_function_njit(self):
        a = np.array([0.0, 1.0, 0.0])
        b = np.array([0.0, 0.0, 0.0])
        w = np.array([1.0, 1.0, 1.0])
        result = diversity_function_njit(a, b, w, 1.0)
        assert np.isclose(result, np.sqrt(1.0/3.0))


class TestBaseGraph:

    @pytest.fixture
    def base(self, mock_protocol, mock_logger):
        g = BaseGraph.__new__(BaseGraph)
        g.confs = []
        g.protocol = mock_protocol
        g.graph_type = "IR"
        g.log = mock_logger
        g.invert = False
        g.fwhm_user = None
        g.shift_user = None
        g.read_population = None
        g.definition = 3
        g.interested_area = None
        BaseGraph.__post_init__(g)
        return g

    def test_post_init(self, base):
        assert base.X is not None
        assert len(base.X) == 1000

    def test_check_conf_inactive(self, base, mock_conformer):
        mock_conformer.active = False
        assert base.check_conf(mock_conformer, base.protocol) is False

    def test_check_conf_no_graph_data(self, base, mock_conformer):
        mock_conformer.active = True
        mock_conformer.graphs_data.__contains__.return_value = False
        assert base.check_conf(mock_conformer, base.protocol) is False

    def test_retrieve_data_empty(self, base, mock_protocol):
        base.confs = []
        base.retrieve_data(mock_protocol)
        assert len(base.energies) == 0

    def test_retrieve_data_invert(self, base, mock_protocol, mock_conformer):
        base.invert = True
        base.confs = [mock_conformer]
        mock_conformer.graphs_data.__contains__.return_value = True
        mock_conformer.graphs_data.has_graph_type.return_value = True
        mock_conformer.graphs_data.__getitem__.return_value = MagicMock()
        mock_conformer.graphs_data.__getitem__.return_value.X = np.array([1000.0, 2000.0])
        mock_conformer.graphs_data.__getitem__.return_value.Y = np.array([1.0, 2.0])
        mock_conformer.energies.__getitem__.return_value.Pop = 0.5
        base.retrieve_data(mock_protocol)
        assert len(base.energies) == 2
        assert np.all(base.impulse <= 0)

    def test_normalize(self, base):
        y = np.array([0.0, 3.0, -1.0])
        norm = base.normalize(y)
        assert np.allclose(norm, [0.0, 1.0, -1.0/3.0])

    def test_normalize_with_bounds(self, base):
        y = np.array([0.0, 3.0, -1.0])
        norm = base.normalize(y, idx_min=0, idx_max=2)
        assert np.allclose(norm, [0.0, 1.0, -1.0/3.0])

    def test_set_boundaries_default(self, base):
        base.set_boundaries()
        assert base.shift_bounds == [0.85, 1.05]
        assert base.fwhm_bounds == [4, 20]

    def test_set_boundaries_user_float(self, base):
        base.shift_user = 0.95
        base.fwhm_user = 12.0
        base.set_boundaries()
        assert base.shift_bounds == [0.95, 0.95]
        assert base.fwhm_bounds == [12.0, 12.0]

    def test_set_boundaries_user_list(self, base):
        base.shift_user = [0.8, 1.2]
        base.fwhm_user = [5, 15]
        base.set_boundaries()
        assert base.shift_bounds == [0.8, 1.2]
        assert base.fwhm_bounds == [5, 15]

    def test_compute_spectrum_empty(self, base):
        base.energies = np.array([])
        base.compute_spectrum()
        base.log.spectra_skip.assert_called()

    def test_compute_spectrum_no_ref(self, base, mock_conformer):
        mock_conformer.graphs_data.__contains__.return_value = True
        mock_conformer.graphs_data.has_graph_type.return_value = True
        mock_conformer.graphs_data.__getitem__.return_value = MagicMock()
        mock_conformer.graphs_data.__getitem__.return_value.X = np.array([1000.0, 2000.0])
        mock_conformer.graphs_data.__getitem__.return_value.Y = np.array([1.0, 2.0])
        mock_conformer.energies.__getitem__.return_value.Pop = 1.0
        base.confs = [mock_conformer]
        base.ref = None
        base.convolute = MagicMock(return_value=np.ones(1000))
        base.compute_spectrum()
        assert base.SHIFT == base.defaults.shift

    def test_diversity_function(self, base):
        base.ref = MagicMock()
        base.ref.weight = np.ones(10)
        a = np.ones(10)
        b = np.zeros(10)
        div = base.diversity_function(a, b)
        assert div > 0

    def test_dump_xy_data(self, base, tmp_path):
        import os
        fname = str(tmp_path / "test.xy")
        base.dump_XY_data(np.array([1.0, 2.0]), np.array([0.5, 1.0]), fname)
        assert os.path.exists(fname)
