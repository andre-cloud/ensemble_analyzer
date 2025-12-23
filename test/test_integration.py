"""
Integration Tests.
Simulates the full protocol execution workflow without external QM calls.
"""

import pytest
from unittest.mock import MagicMock, patch
import numpy as np

from ensemble_analyzer._managers.protocol_executor import ProtocolExecutor
from ensemble_analyzer._conformer.conformer import Conformer
from ensemble_analyzer._protocol.protocol import Protocol

@patch("ensemble_analyzer._conformer.conformer.mkdir")
def test_full_workflow_simulation(mock_mkdir, mock_logger):
    """
    Simulates a complete cycle: Input Ensemble -> Protocol (Calc + Pruning) -> Output.
    Verifies that the executor coordinates managers correctly.
    """

    conf1 = Conformer(1, np.array([[0,0,0]]), ("C",))
    conf2 = Conformer(2, np.array([[1,0,0]]), ("C",))
    
    # === MOCKING ENERGY STORE ===
    # Must mock log_info to return a tuple of 7 values
    
    # Conf 1: Stable
    conf1.energies = MagicMock()
    conf1.energies.get_energy.return_value = -100.0
    conf1.energies.log_info.return_value = (-100.0, 0.0, -100.0, 10.0, 0.0, 50.0, 1.0)
    conf1.energies.__getitem__.return_value = MagicMock(E=-100.0, G=-100.0, zpve=0.1, H=-99.9)

    # Conf 2: High Energy (Unstable)
    conf2.energies = MagicMock()
    conf2.energies.get_energy.return_value = -99.968
    conf2.energies.log_info.return_value = (-99.968, 0.0, -99.968, 10.0, 20.0, 50.0, 1.0)
    conf2.energies.__getitem__.return_value = MagicMock(E=-99.968, G=-99.968, zpve=0.1, H=-99.8)

    # Replace the actual method on instances to avoid StopIteration on side_effects
    # This ensures that no matter how many times get_energy is called, it returns the fixed value
    conf1.get_energy = MagicMock(return_value=-100.0)
    conf2.get_energy = MagicMock(return_value=-99.968)

    ensemble = [conf1, conf2]
    
    config = MagicMock(include_H=True, temperature=298.15, invert=False)
    
    # Setup Protocol with all required attributes for logging
    protocol = MagicMock(spec=Protocol)
    protocol.number = 1
    protocol.thrGMAX = 10.0
    protocol.thrG = 0.5
    protocol.thrB = 10.0
    protocol.clustering = False
    protocol.no_prune = False
    protocol.graph = False
    protocol.read_population = False
    protocol.block_on_retention_rate = False
    protocol.functional = "B3LYP"
    protocol.basis = "6-31G*"
    protocol.calculation_level = "DFT"
    protocol.verbal_internals.return_value = []
    protocol.calculator = "orca"
    protocol.get_calculator.return_value = (MagicMock(), "opt")
    protocol.monitor_internals = []
    
    ckpt_mgr = MagicMock()
    executor = ProtocolExecutor(config, mock_logger, ckpt_mgr)
    
    with patch.object(executor.calculator, 'execute', return_value=True) as mock_calc:
        with patch("ensemble_analyzer._managers.protocol_executor.main_spectra"):
            with patch("ensemble_analyzer._managers.protocol_executor.save_snapshot"):
                
                # EXECUTE
                executor.execute(ensemble, protocol)
                
                # ASSERTIONS
                assert mock_calc.call_count == 2
                
                # Conf2 should be deactivated (Energy > thrGMAX)
                assert conf1.active is True
                assert conf2.active is False
                
                assert mock_logger.pruning_summary.called
                ckpt_mgr.save.assert_called()