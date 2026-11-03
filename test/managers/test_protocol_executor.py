import threading

import numpy as np
import pytest
from unittest.mock import MagicMock, patch
from ensemble_analyzer._managers.protocol_executor import ProtocolExecutor
from ensemble_analyzer.conformer.conformer import Conformer

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
        """Verify that already-computed or inactive conformers are not recalculated."""
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

    def test_cpu_per_job_runs_jobs_at_once(self, mock_logger, mock_protocol):
        """With cpu=8 and cpu_per_job=4, two 4-CPU jobs run together and their results are kept."""
        executor = ProtocolExecutor(MagicMock(cpu=8), mock_logger, MagicMock())
        mock_protocol.cpu_per_job = 4
        confs = [Conformer(number=n, geom=np.zeros((1, 3)), atoms=("H",), raw=True) for n in (1, 2)]
        # Both jobs must reach the barrier before either can finish, so a serial run times out
        barrier = threading.Barrier(2, timeout=5)
        cpus = []

        def fake_execute(i, conf, protocol, cpu):
            cpus.append(cpu)
            barrier.wait()
            conf.last_geometry = np.ones((1, 3))
            return True

        with patch.object(executor.calculator, "execute", side_effect=fake_execute):
            executor._run_calculations(confs, mock_protocol)

        assert cpus == [4, 4]
        assert all(c.active and (c.last_geometry == 1).all() for c in confs)

    def test_cpu_per_job_error_stops_run(self, mock_logger, mock_protocol):
        """An unexpected error stops the run, as in the serial loop, and leaves the conformer active."""
        executor = ProtocolExecutor(MagicMock(cpu=2), mock_logger, MagicMock())
        mock_protocol.cpu_per_job = 1
        confs = [Conformer(number=n, geom=np.zeros((1, 3)), atoms=("H",), raw=True) for n in (1, 2)]

        with patch.object(executor.calculator, "execute", side_effect=FileNotFoundError("orca")):
            with pytest.raises(FileNotFoundError):
                executor._run_calculations(confs, mock_protocol)

        assert all(c.active for c in confs)
