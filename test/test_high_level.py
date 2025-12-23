import pytest
from unittest.mock import MagicMock, patch
from ensemble_analyzer.clustering import get_ensemble
from ensemble_analyzer.graph import main_spectra

def test_get_ensemble_wrapper(mock_logger, mock_conformer):
    ensemble = [mock_conformer]
    with patch("ensemble_analyzer.clustering.ClusteringManager") as MockManager:
        MockManager.return_value.reduce_by_clusters.return_value = ensemble
        res = get_ensemble(ensemble, mock_logger, sort=True)
        assert len(res) == 1
