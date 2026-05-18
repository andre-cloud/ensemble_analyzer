import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from ensemble_analyzer._managers.pruning_manager import PruningManager

class TestPruningManager:
    
    @pytest.fixture
    def manager(self, mock_logger):
        return PruningManager(logger=mock_logger, include_H=True)

    def test_filter_by_energy_window(self, manager, mock_conformer):
        c1 = mock_conformer
        c2 = MagicMock()
        c2.active = True
        c2.number = 2
        
        # Explicitly return float to allow numpy array creation
        c2.get_energy.return_value = -99.9 
        c2.energies.get_energy.return_value = -99.9 
        
        ensemble = [c1, c2]
        
        manager._filter_by_energy_window(ensemble, protocol_number=1, threshold=5.0)
        
        assert c1.active is True
        assert c2.active is False
        manager.logger.table.assert_called()

    def test_remove_duplicates_logic(self, manager, mock_protocol):
        c1 = MagicMock()
        c1.number = 1
        c1.active = True
        c1.energies.__getitem__().B = 0 
        c1.energies.get_energy.return_value = -100.0
        c1.get_energy.return_value = -100.0
        c1.rotatory = 10.0 # Scalar
        c1.moment = 1.0

        c1.distance_matrix.return_value = np.zeros((3,3))

        c2 = MagicMock()
        c2.number = 2
        c2.active = True
        c2.energies.get_energy.return_value = -100.0
        c2.get_energy.return_value = -100.0
        c2.rotatory = 10.01 # Scalar
        c2.moment = 1.0
        c2.distance_matrix.return_value = np.zeros((3,3))
        
        ensemble = [c1, c2]
        
        mock_protocol.thrG = 1.0 
        mock_protocol.thrB = 1.0 
        mock_protocol.thrGMAX = 1000 
        
        with patch.object(manager, '_calculate_rmsd', return_value=0.0):
             manager._remove_duplicates(ensemble, mock_protocol)
        
        assert c2.active is False
        assert c2.diactivated_by == c1.number

    def test_boltzmann_distribution(self, manager, mock_conformer):
        c1 = mock_conformer 
        c2 = MagicMock(active=True)
        # Return float for numpy math
        c2.get_energy.return_value = -100.00159
        c2.energies.get_energy.return_value = -100.00159

        ensemble = [c1, c2]
        proto = MagicMock(number=1)
        
        manager.calculate_relative_energies(ensemble, temperature=298.15, protocol=proto)
        
        assert c2.energies.last().Erel == 0.0
        assert c2.energies.last().Pop > 50.0