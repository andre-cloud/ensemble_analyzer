from unittest.mock import MagicMock

import numpy as np
import pytest

from ensemble_analyzer._nmr.plot import nmr_spectrum_files, plot_nmr_comparisons


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    x = np.linspace(0, 10, 500)
    for p, c in ((1, 4.0), (3, 4.2)):
        np.savetxt(tmp_path / f"NMR_1H_p{p}_comp.xy", np.column_stack((x, np.exp(-((x - c) / 0.05) ** 2))))
    np.savetxt(tmp_path / "NMR_13C_p3_comp.xy", np.column_stack((x * 20, np.exp(-((x - 5) / 0.05) ** 2))))
    np.savetxt(tmp_path / "IR_p1_comp.xy", np.column_stack((x, x)))
    return tmp_path


def test_files_grouped_by_nucleus(workdir):
    files = nmr_spectrum_files()
    assert sorted(files) == ["13C", "1H"] and sorted(files["1H"]) == [1, 3]
    assert sorted(nmr_spectrum_files(idxs=[3])["1H"]) == [3]


def test_comparison_plots(workdir):
    np.savetxt(workdir / "nmr_1h_ref.dat", np.column_stack((np.linspace(0, 10, 50), np.ones(50))))
    plot_nmr_comparisons(MagicMock())
    for nuc in ("1H", "13C"):
        assert (workdir / f"NMR_{nuc}_comparison.png").exists()
        assert (workdir / f"NMR_{nuc}_comparison.pickle").exists()
