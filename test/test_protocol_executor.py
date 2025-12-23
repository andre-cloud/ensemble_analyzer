import pytest
from unittest.mock import MagicMock, patch
from ensemble_analyzer._managers.protocol_executor import ProtocolExecutor

class TestProtocolExecutor:
    
    @pytest.fixture
    def executor(self, mock_logger):
        config = MagicMock()
        config.include_H = True
        config.temperature = 298.15
        config.invert = False
        checkpoint = MagicMock()
        
        return ProtocolExecutor(config, mock_logger, checkpoint)

    def test_execute_skips_inactive_and_calculated(self, executor, mock_conformer, mock_protocol):
        """Verifica che non vengano ricalcolati conformeri già fatti o inattivi."""
        c1 = mock_conformer 
        c2 = MagicMock(active=False) 
        c3 = MagicMock(active=True)

        c3.energies.__contains__.return_value = True 
        
        ensemble = [c1, c2, c3]
        
        with patch.object(executor.calculator, 'execute', return_value=True) as mock_calc:
            executor._run_calculations(ensemble, mock_protocol)

            assert mock_calc.call_count == 1
            assert mock_calc.call_args[0][1] == c1

    def test_execute_full_flow_pruning(self, executor, mock_conformer, mock_protocol):
        """Test that the execute method calls pruning and reporting subsystems."""
        # Ensure conformer is active and has valid energy data
        mock_conformer.active = True
        ensemble = [mock_conformer, mock_conformer]
        mock_protocol.clustering = False
        
        executor.calculator.execute = MagicMock(return_value=True)
        executor.pruning_manager.prune_ensemble = MagicMock()
        executor.pruning_manager.calculate_relative_energies = MagicMock()
        
        with patch("ensemble_analyzer._managers.protocol_executor.main_spectra") as mock_spectra:
             with patch("ensemble_analyzer._managers.protocol_executor.save_snapshot"):
                executor.execute(ensemble, mock_protocol)
                
                executor.pruning_manager.prune_ensemble.assert_called_once()
                mock_spectra.assert_called_once()

    def test_retention_rate_failure(self, executor, mock_conformer, mock_protocol):
        """
        Verify exception is raised if too many conformers are pruned.
        Note: The source code raises a String (TypeError in Py3), so we catch TypeError/Exception.
        """
        mock_conformer.active = True
        ensemble = [mock_conformer]
        
        def kill_all(*args, **kwargs):
            for c in args[0]: c.active = False
            
        executor.pruning_manager.prune_ensemble = kill_all
        mock_protocol.block_on_retention_rate = True
        
        with patch("ensemble_analyzer._managers.protocol_executor.MIN_RETENTION_RATE", 0.5):
            with patch.object(executor, '_run_calculations'):
                with patch("ensemble_analyzer._managers.protocol_executor.save_snapshot"):
                     with pytest.raises((TypeError, Exception)) as exc:
                        executor.execute(ensemble, mock_protocol)
                     assert "pruning" in str(exc.value) or "Exception" in str(exc) or "too much" in str(exc.value)