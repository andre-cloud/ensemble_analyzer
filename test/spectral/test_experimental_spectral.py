import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from ensemble_analyzer._spectral.experimental import ExperimentalGraph


class TestExperimentalGraph:

    @pytest.fixture
    def exp_graph(self, mock_protocol, mock_logger):
        g = ExperimentalGraph(
            confs=[], protocol=mock_protocol, graph_type="IR",
            log=mock_logger, definition=3
        )
        return g

    def test_post_init_basics(self, exp_graph):
        assert exp_graph.graph_type == "IR"
        assert exp_graph.X is not None
        assert len(exp_graph.X) == 1000
        assert exp_graph.X[0] <= exp_graph.X[-1]

    @patch("ensemble_analyzer._spectral.experimental.np.loadtxt")
    def test_load_file_experimental(self, mock_loadtxt, exp_graph):
        mock_loadtxt.return_value = np.column_stack([
            np.array([500.0, 1500.0, 1000.0]),
            np.array([0.5, 1.0, 0.8])
        ])
        exp_graph.load_file_experimental()
        assert exp_graph.x is not None
        assert exp_graph.y is not None
        assert len(exp_graph.x) == 3

    def test_interpolate(self, exp_graph):
        exp_graph.x = np.array([500.0, 1000.0, 1500.0])
        exp_graph.y = np.array([0.0, 1.0, 0.0])
        Y = exp_graph.interpolate()
        assert len(Y) == len(exp_graph.X)

    def test_process_ir(self, exp_graph):
        exp_graph.x = np.array([500.0, 1000.0, 1500.0])
        exp_graph.y = np.array([0.0, 1.0, 0.0])
        exp_graph.process()
        assert hasattr(exp_graph, "x_min")
        assert hasattr(exp_graph, "x_max")
        assert hasattr(exp_graph, "x_min_idx")
        assert hasattr(exp_graph, "x_max_idx")

    def test_process_uv_with_conversion(self, mock_protocol, mock_logger):
        g = ExperimentalGraph(
            confs=[], protocol=mock_protocol, graph_type="UV",
            log=mock_logger, definition=3
        )
        g.x = np.array([400.0, 500.0, 600.0])
        g.y = np.array([0.5, 1.0, 0.8])
        g.process()
        assert np.any(g.x < 100)

    def test_calc_weighting_no_interest(self, exp_graph):
        exp_graph.X = np.linspace(0, 100, 1001)
        exp_graph.x_min_idx = 100
        exp_graph.x_max_idx = 900
        exp_graph.interested_area = None
        exp_graph.calc_weighting_function()
        assert np.sum(exp_graph.weight > 0) > 0

    def test_calc_weighting_with_interest_area_true(self, exp_graph):
        exp_graph.X = np.linspace(0, 200, 2001)
        exp_graph.x_min_idx = 200
        exp_graph.x_max_idx = 1800
        exp_graph.interested_area = True
        exp_graph.calc_weighting_function()
        assert np.sum(exp_graph.weight > 0) > 0

    def test_calc_weighting_with_list_interest(self, exp_graph):
        exp_graph.X = np.linspace(0, 200, 2001)
        exp_graph.x_min_idx = 200
        exp_graph.x_max_idx = 1800
        exp_graph.interested_area = [500, 1500]
        exp_graph.calc_weighting_function()
        assert np.sum(exp_graph.weight > 0) > 0

    def test_calc_weighting_with_scalar_interest(self, exp_graph):
        exp_graph.X = np.linspace(0, 200, 2001)
        exp_graph.x_min_idx = 200
        exp_graph.x_max_idx = 1800
        exp_graph.interested_area = 1000
        exp_graph.calc_weighting_function()
        assert np.sum(exp_graph.weight > 0) > 0

    def test_calc_weighting_invalid_interest(self, exp_graph):
        exp_graph.X = np.linspace(0, 200, 2001)
        exp_graph.x_min_idx = 200
        exp_graph.x_max_idx = 1800
        exp_graph.interested_area = "invalid"
        with pytest.raises(ValueError, match="Invalid format"):
            exp_graph.calc_weighting_function()

    def test_gau(self, exp_graph):
        exp_graph.X = np.linspace(-10, 10, 1001)
        g = exp_graph.gau(0.0, 1.0)
        assert np.isclose(np.max(g), 1/(1.0*np.sqrt(2*np.pi)), atol=0.01)
