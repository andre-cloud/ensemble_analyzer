"""
Tests for Clustering and PCA module.
"""
import pytest
import numpy as np
from unittest.mock import patch, MagicMock
from ensemble_analyzer.clustering import validate_possible_PCA, execute_PCA
from ensemble_analyzer._clustering.cluster_manager import ClusteringManager, ClusteringConfig, PCAResult

class TestClustering:
    
    @pytest.fixture
    def manager(self, mock_logger):
        config = ClusteringConfig(include_H=True, set_cluster_attribute=True)
        return ClusteringManager(logger=mock_logger, config=config)

    @pytest.mark.parametrize("n_confs, n_clusters, expected", [
        (5, 2, True), (1, None, False), (3, 4, False)
    ])
    def test_validate_possible_pca(self, n_confs, n_clusters, expected, mock_logger, mock_conformer):
        ensemble = [mock_conformer] * n_confs
        with patch("ensemble_analyzer.clustering.MIN_CONFORMERS_FOR_PCA", 2):
            assert validate_possible_PCA(ensemble, mock_logger, n_clusters) is expected

    @patch("ensemble_analyzer.clustering.ClusteringManager")
    def test_execute_pca_config_passing(self, MockManager, mock_logger, mock_conformer):
        with patch("ensemble_analyzer.clustering.validate_possible_PCA", return_value=True):
            execute_PCA([mock_conformer], 3, "out.png", "Title", mock_logger, False, False)
            args, kwargs = MockManager.call_args
            assert kwargs['config'].n_clusters == 3

    def test_reduce_by_clusters(self, manager, mock_conformer):
        c1 = MagicMock(active=True, cluster=0)
        c1.__lt__ = lambda s, o: True 
        c2 = MagicMock(active=True, cluster=0)
        c3 = MagicMock(active=True, cluster=1)
        ensemble = [c1, c2, c3]
        reduced = manager.reduce_by_clusters(ensemble)
        active = [c for c in reduced if c.active]
        assert len(active) == 2 

    def test_perform_pca_full_pipeline(self, manager, mock_conformer):
        """Test full PCA pipeline execution including plotting hooks."""
        mock_conformer.last_geometry = np.zeros((3, 3))
        mock_conformer.atoms = ["C", "H", "H"]
        # Ensure get_energy returns float for numpy
        mock_conformer.energies.get_energy.return_value = -100.0
        mock_conformer.color = "#ffffff"
        mock_conformer.number = 1
        
        ensemble = [mock_conformer, mock_conformer]
        
        # Use context managers for clear scope and order
        with patch("ensemble_analyzer._clustering.cluster_manager.PCA") as MockPCA, \
             patch("ensemble_analyzer._clustering.cluster_manager.KMeans") as MockKMeans, \
             patch.object(manager, '_create_visualization') as MockViz:
            
            MockPCA.return_value.fit_transform.return_value = np.random.rand(2, 2)
            MockPCA.return_value.explained_variance_ratio_ = np.array([0.8, 0.1])
            MockKMeans.return_value.fit_predict.return_value = [0, 1]
            
            result = manager.perform_pca(ensemble, n_clusters=2)
            
            assert isinstance(result, PCAResult)
            assert result.n_clusters == 2
            MockViz.assert_called()

    def test_calculate_distance_matrix_eigenvalues(self, manager):
        geom = np.array([[[0,0,0], [1,0,0]], [[0,0,0], [2,0,0]]])
        atoms = np.array(["C", "C"])
        feats = manager._calculate_distance_matrix_eigenvalues(geom, atoms, include_H=True)
        assert feats.shape == (2, 2)
        assert np.allclose(sorted(feats[0]), [-1.0, 1.0])