import pytest
import numpy as np
from unittest.mock import MagicMock, patch, call
from ensemble_analyzer._managers.calculation_executor import CalculationExecutor
from ensemble_analyzer.constants import EV_TO_EH, regex_parsing

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

    @patch("ensemble_analyzer._managers.calculation_executor.compute_rotational_constants")
    @patch("ensemble_analyzer._managers.calculation_executor.get_conf_parameters")
    def test_ml_energy_unit_conversion(self, mock_get_params, mock_rot_const, executor, mock_conformer, mock_protocol):
        mock_get_params.return_value = True
        mock_protocol.calculator = "tblite"
        mock_protocol.freq = False
        mock_protocol.opt = False
        mock_calc = MagicMock()
        mock_protocol.get_calculator.return_value = (mock_calc, "label")

        mock_atoms = MagicMock()
        mock_atoms.get_potential_energy.return_value = -27.211386245981
        mock_atoms.get_dipole_moment.return_value = np.array([0.0, 0.0, 0.0])
        mock_conformer.get_ase_atoms.return_value = mock_atoms

        conf_energies_mock = MagicMock()
        conf_energies_mock.__contains__.return_value = False
        mock_conformer.energies = conf_energies_mock

        success = executor.execute(1, mock_conformer, mock_protocol)

        assert success is True
        conf_energies_mock.add.assert_called_once()
        args, _ = conf_energies_mock.add.call_args
        record = args[1]
        assert record.E == pytest.approx(-1.0, abs=1e-10), (
            f"Expected -1.0 Eh (from -27.2114 eV), got {record.E}"
        )

    @patch("ensemble_analyzer._managers.calculation_executor.get_conf_parameters")
    def test_execute_failure_parsing(self, mock_get_params, executor, mock_conformer, mock_protocol):
        mock_get_params.return_value = False
        
        # Ensure correct return type for unpacking: calc, label = ...
        mock_protocol.get_calculator.return_value = (MagicMock(), "label")
        
        mock_conformer.get_ase_atoms.return_value = MagicMock()
        
        success = executor.execute(1, mock_conformer, mock_protocol)
        
        assert success is False
        executor.logger.calculation_success.assert_not_called()