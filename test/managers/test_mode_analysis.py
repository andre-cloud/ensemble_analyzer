import pytest
import numpy as np
from ensemble_analyzer.mode_analysis import NormalModeAnalyzer


class TestIsNullMode:

    @pytest.fixture
    def analyzer(self):
        n_atoms, n_modes = 6, 8
        rng = np.random.default_rng(42)
        modes = rng.normal(0, 0.1, (n_modes, n_atoms, 3))
        return NormalModeAnalyzer(
            modes, np.zeros((n_atoms, 3)), ("C",) * n_atoms
        )

    def test_returns_true_for_zero_displacement(self, analyzer):
        modes = analyzer.normal_modes.copy()
        modes[3] = np.zeros((analyzer.n_atoms, 3))
        analyzer.normal_modes = modes
        assert analyzer.is_null_mode(3)

    def test_returns_false_for_non_zero_displacement(self, analyzer):
        assert not analyzer.is_null_mode(0)

    def test_returns_true_when_all_modes_zero(self, analyzer):
        modes = np.zeros((analyzer.n_modes, analyzer.n_atoms, 3))
        analyzer.normal_modes = modes
        for i in range(analyzer.n_modes):
            assert analyzer.is_null_mode(i)
