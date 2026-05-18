import os
import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from ensemble_analyzer.graph import eV_to_nm, main_spectra, plot_comparative_graphs


class TestGraph:

    def test_eV_to_nm(self):
        result = eV_to_nm(np.array([1.0, 2.0]))
        assert len(result) == 2
        assert result[0] > result[1]

    @patch("ensemble_analyzer.graph.os.path.exists")
    @patch("ensemble_analyzer.graph.ExperimentalGraph")
    def test_main_spectra_no_ref(self, MockExp, mock_exists, mock_protocol, mock_logger):
        mock_exists.return_value = False
        mock_cls = MagicMock()
        fake_class = {"IR": mock_cls, "VCD": MagicMock(), "UV": MagicMock(), "ECD": MagicMock()}
        with patch.dict("ensemble_analyzer.graph.class_", fake_class, clear=True):
            main_spectra([], mock_protocol, mock_logger, invert=False, interested_area={},
                         shift={"vibro": None, "electro": None}, fwhm={"vibro": None, "electro": None})
            mock_cls.assert_called()
            mock_cls.return_value.compute_spectrum.assert_called()

    @patch("ensemble_analyzer.graph.os.path.exists")
    @patch("ensemble_analyzer.graph.ExperimentalGraph")
    def test_main_spectra_with_ref(self, MockExp, mock_exists, mock_protocol, mock_logger):
        def exists_side(path):
            return "ref.dat" in path
        mock_exists.side_effect = exists_side
        mock_cls = MagicMock()
        fake_class = {"IR": mock_cls}
        with patch.dict("ensemble_analyzer.graph.class_", fake_class, clear=True):
            main_spectra([], mock_protocol, mock_logger, invert=False,
                         interested_area={"vibro": None, "electro": None},
                         shift={"vibro": None, "electro": None}, fwhm={"vibro": None, "electro": None})
            MockExp.assert_called()
            mock_cls.assert_called()

    @patch("ensemble_analyzer.graph.ComparedGraph")
    @patch("ensemble_analyzer.graph.os.path.exists")
    def test_plot_comparative_graphs(self, mock_exists, MockCompared, mock_logger):
        mock_exists.return_value = False
        comp_instance = MagicMock()
        comp_instance.data = [1]
        MockCompared.return_value = comp_instance
        plot_comparative_graphs(mock_logger, show=True, nm=True)
        comp_instance.plot.assert_called_with(show=True, show_ref_weight=False)

    @patch("ensemble_analyzer.graph.ComparedGraph")
    @patch("ensemble_analyzer.graph.os.path.exists")
    def test_plot_comparative_graphs_empty(self, mock_exists, MockCompared, mock_logger):
        mock_exists.return_value = True
        comp_instance = MagicMock()
        comp_instance.data = []
        MockCompared.return_value = comp_instance
        plot_comparative_graphs(mock_logger)
        comp_instance.plot.assert_not_called()
