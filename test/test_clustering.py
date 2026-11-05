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

# ---------------------------------------------------------------------------
# Automatic k search (issue #3)
# ---------------------------------------------------------------------------

from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import cdist

from ensemble_analyzer._clustering.cluster_manager import _candidate_ks


def _blobs(centers, n=300, seed=0):
    X, _ = make_blobs(n_samples=n, centers=centers, n_features=8, cluster_std=1.0, random_state=seed)
    return X


class TestOptimalClusters:

    @pytest.mark.parametrize("centers", [5, 20])
    def test_finds_true_number_of_clusters(self, mock_logger, centers):
        mgr = ClusteringManager(mock_logger, ClusteringConfig())
        k, k_range, scores = mgr.find_optimal_clusters(_blobs(centers))
        assert k == centers
        assert len(k_range) == len(scores)

    @pytest.mark.parametrize("cfg, n, bounds", [
        (dict(), 300, (2, 240)),
        (dict(min_k=5), 300, (5, 240)),
        (dict(max_k=30), 300, (2, 30)),
        (dict(max_k_fraction=0.5), 300, (2, 150)),
        (dict(max_k_fraction=1.0), 40, (2, 39)),  # capped at N-1
        (dict(min_k=0), 30, (2, 24)),
    ])
    def test_k_bounds(self, mock_logger, cfg, n, bounds):
        assert ClusteringManager(mock_logger, ClusteringConfig(**cfg))._k_bounds(n) == bounds

    def test_search_stays_in_bounds(self, mock_logger):
        mgr = ClusteringManager(mock_logger, ClusteringConfig(min_k=3, max_k=12))
        _, k_range, _ = mgr.find_optimal_clusters(_blobs(5, n=100))
        assert list(k_range) == list(range(3, 13))

    def test_parallel_matches_serial(self, mock_logger, monkeypatch):
        monkeypatch.setattr("ensemble_analyzer._clustering.cluster_manager._PARALLEL_MIN_SAMPLES", 0)
        X = _blobs(7, n=200)
        serial = ClusteringManager(mock_logger, ClusteringConfig(n_jobs=1)).find_optimal_clusters(X)
        parallel = ClusteringManager(mock_logger, ClusteringConfig(n_jobs=2)).find_optimal_clusters(X)
        assert serial[0] == parallel[0]
        assert list(serial[1]) == list(parallel[1])
        np.testing.assert_allclose(serial[2], parallel[2], atol=1e-8)

    def test_precomputed_silhouette_matches_default(self):
        X = _blobs(4, n=120)
        labels = KMeans(n_clusters=4, n_init="auto", random_state=42).fit_predict(X)
        D = cdist(X, X)
        np.fill_diagonal(D, 0.0)
        assert silhouette_score(D, labels, metric="precomputed") == pytest.approx(silhouette_score(X, labels))

    def test_scores_match_reference_loop(self, mock_logger):
        X = _blobs(5, n=80)
        _, k_range, scores = ClusteringManager(mock_logger, ClusteringConfig()).find_optimal_clusters(X)
        for k, s in zip(k_range[:10], scores[:10]):
            labels = KMeans(n_clusters=k, n_init="auto", random_state=42).fit_predict(X)
            assert s == pytest.approx(silhouette_score(X, labels))

    def test_adaptive_grid(self):
        ks, step = _candidate_ks(2, 100)
        assert ks == list(range(2, 101)) and step == 1
        ks, step = _candidate_ks(2, 2400)
        assert ks[:49] == list(range(2, 51))
        assert ks[-1] == 2400
        assert step > 1 and len(ks) < 200

    def test_wide_range_refines_around_best(self, mock_logger, monkeypatch):
        monkeypatch.setattr("ensemble_analyzer._clustering.cluster_manager._PARALLEL_MIN_SAMPLES", 0)
        X = _blobs(120, n=1000, seed=1)
        mgr = ClusteringManager(mock_logger, ClusteringConfig(n_jobs=2))
        k, k_range, scores = mgr.find_optimal_clusters(X)
        assert len(k_range) < 300
        assert scores[list(k_range).index(k)] == max(scores)
        # the refinement around the coarse optimum reaches the true k (120) closely
        assert abs(k - 120) <= 2
        assert {k - 1, k, k + 1} <= set(k_range)

    def test_too_few_samples_for_range(self, mock_logger):
        mgr = ClusteringManager(mock_logger, ClusteringConfig(min_k=5))
        k, k_range, scores = mgr.find_optimal_clusters(_blobs(2, n=4))
        assert k == 5 and list(k_range) == [5] and scores == [0.0]


class TestPruningFloor:

    def _confs(self, X):
        confs = []
        for i, x in enumerate(X):
            c = MagicMock(active=True, number=i + 1, color="#000000")
            # 3 atoms: the PCA features are their distance-matrix eigenvalues
            c.last_geometry = np.array([[0, 0, 0], [x[0], 0, 0], [0, x[1], 0]], float)
            c.atoms = np.array(["C", "C", "C"])
            c.energies.get_energy.return_value = float(i)
            confs.append(c)
        return confs

    def _run(self, mock_logger, frac):
        X, _ = make_blobs(n_samples=100, centers=[[2, 2], [6, 6]], n_features=2,
                          cluster_std=0.05, random_state=0)
        confs = self._confs(X)
        mgr = ClusteringManager(mock_logger, ClusteringConfig(min_keep_fraction=frac))
        result = mgr._execute_pca_pipeline(confs, n_clusters=None, include_H=True)
        return result, confs

    def test_no_floor_uses_optimal_k(self, mock_logger):
        result, confs = self._run(mock_logger, 0.0)
        assert result.optimal_k == 2
        assert result.n_clusters == 2
        assert len({c.cluster for c in confs}) == 2

    def test_floor_raises_k_for_pruning(self, mock_logger):
        result, confs = self._run(mock_logger, 0.1)
        assert result.optimal_k == 2          # the real optimum is still reported
        assert result.n_clusters == 10        # ceil(0.1 * 100)
        assert len({c.cluster for c in confs}) == 10

    def test_execute_pca_applies_floor_only_for_pruning(self, mock_logger, mock_conformer):
        with patch("ensemble_analyzer.clustering.ClusteringManager") as M, \
             patch("ensemble_analyzer.clustering.validate_possible_PCA", return_value=True):
            execute_PCA([mock_conformer], None, "o.png", "t", mock_logger, set_=True, n_jobs=4)
            cfg = M.call_args.kwargs["config"]
            assert cfg.min_keep_fraction == pytest.approx(0.1) and cfg.n_jobs == 4
            execute_PCA([mock_conformer], None, "o.png", "t", mock_logger, set_=False)
            assert M.call_args.kwargs["config"].min_keep_fraction == 0.0

    def test_execute_pca_reports_failure(self, mock_logger, mock_conformer):
        with patch("ensemble_analyzer.clustering.ClusteringManager") as M, \
             patch("ensemble_analyzer.clustering.validate_possible_PCA", return_value=True):
            M.return_value.perform_pca.return_value = None
            assert execute_PCA([mock_conformer], None, "o.png", "t", mock_logger) is False


def test_vectorised_features_match_loop(mock_logger):
    rng = np.random.default_rng(0)
    geoms = rng.normal(size=(37, 9, 3))
    atoms = np.array(["C", "H", "O", "H", "C", "N", "H", "C", "H"])
    mgr = ClusteringManager(mock_logger, ClusteringConfig())
    for include_H in (True, False):
        got = mgr._calculate_distance_matrix_eigenvalues(geoms, atoms, include_H)
        mask = slice(None) if include_H else atoms != "H"
        ref = np.array([np.linalg.eigvalsh(cdist(g[mask], g[mask])) for g in geoms])
        np.testing.assert_allclose(got, ref, atol=1e-10)
