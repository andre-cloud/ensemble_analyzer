import numpy as np
from unittest.mock import MagicMock, patch
from ensemble_analyzer.constants import eV_to_nm
from ensemble_analyzer.graph import main_spectra, plot_comparative_graphs


class TestGraph:

    def test_eV_to_nm(self):
        result = eV_to_nm(np.array([1.0, 2.0]))
        assert len(result) == 2
        assert result[0] > result[1]

    @patch("ensemble_analyzer.graph.os.path.exists")
    @patch("ensemble_analyzer.graph.ExperimentalGraph")
    @patch("ensemble_analyzer.graph.BaseGraph")
    def test_main_spectra_no_ref(self, mock_cls, MockExp, mock_exists, mock_protocol, mock_logger):
        mock_exists.return_value = False
        main_spectra([], mock_protocol, mock_logger, invert=False, interested_area={},
                     shift={"vibro": None, "electro": None}, fwhm={"vibro": None, "electro": None})
        assert [c.kwargs["graph_type"] for c in mock_cls.call_args_list] == ["IR", "VCD", "UV", "ECD"]
        mock_cls.return_value.compute_spectrum.assert_called()
        MockExp.assert_not_called()

    @patch("ensemble_analyzer.graph.os.path.exists")
    @patch("ensemble_analyzer.graph.ExperimentalGraph")
    @patch("ensemble_analyzer.graph.BaseGraph")
    def test_main_spectra_with_ref(self, mock_cls, MockExp, mock_exists, mock_protocol, mock_logger):
        mock_exists.side_effect = lambda path: path.endswith("ir_ref.dat")
        main_spectra([], mock_protocol, mock_logger, invert=False,
                     interested_area={"vibro": None, "electro": None},
                     shift={"vibro": None, "electro": None}, fwhm={"vibro": None, "electro": None})
        MockExp.assert_called_once()
        assert mock_cls.call_args_list[0].kwargs["ref"] is MockExp.return_value
        assert mock_cls.call_args_list[2].kwargs["ref"] is None

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
