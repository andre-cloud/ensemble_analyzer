import pytest
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

    @patch("ensemble_analyzer._managers.calculation_executor.move_files")
    @patch("ensemble_analyzer._managers.calculation_executor.get_conf_parameters")
    def test_execute_success(self, mock_get_params, mock_move, executor, mock_conformer, mock_protocol):
        mock_get_params.return_value = True
        
        # Ensure protocol.calculator is a valid key (e.g., 'orca' or 'gaussian')
        # Note: mock_protocol from fixture already sets calculator="orca"
        # but we ensure get_calculator returns a tuple
        mock_calc = MagicMock()
        mock_protocol.get_calculator.return_value = (mock_calc, "label")
        
        mock_atoms = MagicMock()
        mock_conformer.get_ase_atoms.return_value = mock_atoms
        
        success = executor.execute(1, mock_conformer, mock_protocol)
        
        assert success is True
        mock_atoms.get_potential_energy.assert_called_once()
        mock_move.assert_called_once()
        executor.logger.calculation_success.assert_called_once()

    @patch("ensemble_analyzer._managers.calculation_executor.move_files")
    @patch("ensemble_analyzer._managers.calculation_executor.get_conf_parameters")
    def test_execute_failure_parsing(self, mock_get_params, mock_move, executor, mock_conformer, mock_protocol):
        mock_get_params.return_value = False
        
        # Ensure correct return type for unpacking: calc, label = ...
        mock_protocol.get_calculator.return_value = (MagicMock(), "label")
        
        mock_conformer.get_ase_atoms.return_value = MagicMock()
        
        success = executor.execute(1, mock_conformer, mock_protocol)
        
        assert success is False
        executor.logger.calculation_success.assert_not_called()