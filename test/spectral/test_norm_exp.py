"""--norm-exp: spectra scaled to the experimental |max| instead of [-1, 1]."""
import numpy as np
import pytest

from ensemble_analyzer._spectral.experimental import ExperimentalGraph
from ensemble_analyzer._spectral.compare import ComparedGraph


@pytest.mark.parametrize("norm_exp, expected", [(False, 1.0), (True, 40.0)])
def test_experimental_and_computed_scale(tmp_path, monkeypatch, mock_protocol, mock_logger,
                                         norm_exp, expected):
    monkeypatch.chdir(tmp_path)
    x = np.linspace(800, 1800, 200)
    y = 40.0 * np.exp(-((x - 1300) / 60) ** 2) - 10.0 * np.exp(-((x - 1000) / 40) ** 2)
    np.savetxt("vcd_ref.dat", np.column_stack([x, y]))

    ref = ExperimentalGraph(confs=[], protocol=mock_protocol, graph_type="VCD",
                            log=mock_logger, definition=3, norm_exp=norm_exp)
    ref.read()
    _, y_ref = np.loadtxt("VCD_ref_norm.xy", unpack=True)
    assert np.max(np.abs(y_ref)) == pytest.approx(expected, rel=1e-2)

    # A computed spectrum on [-1, 1] is plotted at the experimental scale
    np.savetxt("VCD_p1_comp.xy", np.column_stack([ref.X, -y_ref / np.max(np.abs(y_ref))]))
    comp = ComparedGraph(graph_type="VCD", experimental_file="VCD_ref_norm.xy", log=mock_logger)
    i, j = int(comp.bounders[0]), int(comp.bounders[1])
    assert np.max(np.abs(comp.data["1"][1][i:j])) == pytest.approx(np.max(np.abs(y_ref[i:j])))


@pytest.mark.parametrize("graph_type, label", [
    ("UV", "Abs [a.u.]"), ("ECD", "mdeg [°]"), ("IR", "Abs [a.u.]"), ("VCD", "∆A [a.u.]"),
])
def test_norm_exp_y_labels(tmp_path, monkeypatch, mock_logger, graph_type, label):
    from unittest.mock import MagicMock
    from ensemble_analyzer._spectral.graph_default import GraphDefault
    monkeypatch.chdir(tmp_path)
    np.savetxt(f"{graph_type}_ref_norm.xy", np.column_stack([np.arange(10.0), np.ones(10)]))
    np.savetxt(f"{graph_type}_index_lim.xy", np.array([0, 9]))
    np.savetxt(f"{graph_type}_weighted.xy", np.column_stack([np.arange(10.0), np.ones(10)]))

    for norm_exp, expected in [(True, label), (False, GraphDefault(graph_type).axis_label["y"])]:
        comp = ComparedGraph(graph_type=graph_type, experimental_file=f"{graph_type}_ref_norm.xy",
                             log=mock_logger, norm_exp=norm_exp)
        ax = MagicMock()
        comp._configure_axes(ax, in_nm=False)
        ax.set_ylabel.assert_called_with(expected)
