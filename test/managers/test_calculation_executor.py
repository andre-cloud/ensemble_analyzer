import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from ensemble_analyzer._managers.calculation_executor import CalculationExecutor
from ensemble_analyzer.constants import regex_parsing

class TestCalculationExecutor:

    @pytest.fixture
    def executor(self, mock_logger):
        config = MagicMock()
        config.temperature = 298.15
        config.cpu = 4
        return CalculationExecutor(config, mock_logger)

    @patch("ensemble_analyzer._managers.calculation_executor.compute_rotational_constants")
    @patch("ensemble_analyzer._managers.calculation_executor.get_conf_parameters")
    def test_execute_success(self, mock_get_params, mock_rot_const, executor, mock_conformer, mock_protocol):
        mock_get_params.return_value = True

        # Use an ML calculator to test the ML execution path
        mock_protocol.calculator = "tblite"
        mock_calc = MagicMock()
        mock_protocol.get_calculator.return_value = (mock_calc, "label")

        mock_atoms = MagicMock()
        mock_conformer.get_ase_atoms.return_value = mock_atoms

        success = executor.execute(1, mock_conformer, mock_protocol)

        assert success is True
        mock_atoms.get_potential_energy.assert_called_once()
        executor.logger.calculation_success.assert_called_once()

    @patch("ensemble_analyzer._managers.calculation_executor.get_conf_parameters")
    def test_execute_failure_parsing(self, mock_get_params, executor, mock_conformer, mock_protocol):
        mock_get_params.return_value = False
        
        # Ensure correct return type for unpacking: calc, label = ...
        mock_protocol.get_calculator.return_value = (MagicMock(), "label")
        
        mock_conformer.get_ase_atoms.return_value = MagicMock()
        
        success = executor.execute(1, mock_conformer, mock_protocol)
        
        assert success is False
        executor.logger.calculation_success.assert_not_called()

    # ------------------------------------------------------------------
    # Null normal mode skipping
    # ------------------------------------------------------------------

    @pytest.fixture
    def null_mode_data(self):
        n_atoms = 6
        freqs = np.array([0.0, 0.0, 0.0, 0.0, 0.0, -108.67, -87.21, 103.76])
        rng = np.random.default_rng(42)
        modes = rng.normal(0, 0.1, (8, n_atoms, 3))
        modes[:5] = 0.0
        modes[5] = 0.0
        modes[6] = 0.0
        modes[6, :3] = [0.5, 0.5, 0.5]
        return freqs, modes

    @pytest.fixture
    def ts_protocol(self):
        proto = MagicMock()
        proto.number = 1
        proto.ts = True
        proto.neg_freq_threshold = 20.0
        proto.ts_target = {"reactive": [0, 1, 2]}
        proto.min_overlap = 50.0
        proto.displace_scale = 0.3
        return proto

    @pytest.fixture
    def opt_protocol(self):
        proto = MagicMock()
        proto.number = 1
        proto.ts = False
        proto.auto_displace = True
        proto.neg_freq_threshold = 20.0
        proto.displace_scale = 0.3
        return proto

    @pytest.fixture
    def conf_with_data(self, null_mode_data):
        freqs, modes = null_mode_data
        n_atoms = modes.shape[1]
        data = MagicMock()
        data.Freq = freqs
        data.NormalModes = modes
        conf = MagicMock()
        conf.number = 1
        conf.active = True
        conf.last_geometry = np.zeros((n_atoms, 3))
        conf.atoms = tuple(f"A{i}" for i in range(n_atoms))
        conf.energies = {1: data}
        return conf

    @pytest.mark.parametrize("calculator", ["orca", "gaussian", "nwchem", "tblite", "aimnet"])
    def test_ts_skips_null_modes(
        self, executor, conf_with_data, ts_protocol, calculator
    ):
        ts_protocol.calculator = calculator
        result = executor._check_imaginary_and_displace(conf_with_data, ts_protocol)
        assert result == (False, None)
        log_messages = [str(c) for c in executor.logger.info.call_args_list]
        assert any("null normal mode(s)" in msg for msg in log_messages)

    def test_ts_null_mode_logs_index_and_freq(
        self, executor, conf_with_data, ts_protocol
    ):
        executor._check_imaginary_and_displace(conf_with_data, ts_protocol)
        log_messages = [str(c) for c in executor.logger.info.call_args_list]
        null_log = [m for m in log_messages if "null normal mode(s)" in m][0]
        assert "5" in null_log
        assert "-108.67" in null_log

    def test_ts_all_modes_null_returns_false(
        self, executor, conf_with_data, ts_protocol
    ):
        conf_with_data.energies[1].NormalModes[:] = 0.0
        result = executor._check_imaginary_and_displace(conf_with_data, ts_protocol)
        assert result == (False, None)

    def test_opt_skips_null_modes(
        self, executor, conf_with_data, opt_protocol
    ):
        result = executor._check_imaginary_and_displace(conf_with_data, opt_protocol)
        success, new_geom = result
        assert success is True
        assert new_geom is not None

    def test_no_null_modes_no_skip_log(
        self, executor, ts_protocol
    ):
        n_atoms = 6
        freqs = np.array([0.0, 0.0, -50.0, -30.0, 100.0, 200.0])
        modes = np.zeros((6, n_atoms, 3))
        modes[2] = np.random.default_rng(0).normal(0, 0.1, (n_atoms, 3))
        modes[3] = np.random.default_rng(1).normal(0, 0.1, (n_atoms, 3))
        data = MagicMock()
        data.Freq = freqs
        data.NormalModes = modes
        conf = MagicMock()
        conf.number = 1
        conf.active = True
        conf.last_geometry = np.zeros((n_atoms, 3))
        conf.atoms = ("C",) * n_atoms
        conf.energies = {1: data}
        executor._check_imaginary_and_displace(conf, ts_protocol)
        log_messages = [str(c) for c in executor.logger.info.call_args_list]
        null_logs = [m for m in log_messages if "null normal mode(s)" in m]
        assert len(null_logs) == 0