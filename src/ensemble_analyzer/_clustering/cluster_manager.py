
from ensemble_analyzer._logger.logger import Logger
from ensemble_analyzer.conformer.conformer import Conformer

from pathlib import Path

from ensemble_analyzer.constants import EH_TO_KCAL, DEFAULT_RESOLUTION, MARKERS

from .cluster_config import ClusteringConfig, PCAResult


from typing import Optional, List, Dict

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.gridspec as gridspec

from scipy.interpolate import griddata
from scipy.spatial.distance import cdist  # scipy.spatial.distance_matrix is deprecated since SciPy 1.18
from joblib import Parallel, delayed
from threadpoolctl import threadpool_limits


# Silhouette search grid: every k up to _DENSE_K, then ~_COARSE_POINTS coarse values
_DENSE_K = 50
_MAX_DENSE_RANGE = 150
_COARSE_POINTS = 100
# Below this ensemble size starting the worker processes costs more than it saves
_PARALLEL_MIN_SAMPLES = 1000


def _candidate_ks(lo: int, hi: int) -> tuple[list[int], int]:
    """k values for the first pass of the silhouette search, and the coarse step.

    Narrow ranges are scanned completely (step 1). Wide ranges are scanned
    densely up to _DENSE_K and then with a step chosen to give about
    _COARSE_POINTS values; the caller refines around the best coarse k.
    """
    if hi - lo + 1 <= _MAX_DENSE_RANGE:
        return list(range(lo, hi + 1)), 1
    dense_hi = max(lo, min(_DENSE_K, hi))
    step = max(2, int(np.ceil((hi - dense_hi) / _COARSE_POINTS)))
    ks = list(range(lo, dense_hi + 1)) + list(range(dense_hi + step, hi + 1, step))
    if ks[-1] != hi:
        ks.append(hi)
    return ks, step


def _score_k(features, dist, k, random_state, single_thread):
    """KMeans with k clusters and its silhouette score on a precomputed distance matrix."""
    def run():
        labels = KMeans(n_clusters=k, n_init='auto', random_state=random_state).fit_predict(features)
        if len(np.unique(labels)) < 2:
            return k, -1.0, labels
        return k, float(silhouette_score(dist, labels, metric="precomputed")), labels

    if single_thread:
        # parallel workers: avoid oversubscription by OpenMP/BLAS threads
        with threadpool_limits(1):
            return run()
    return run()


class ClusteringManager:
    """
    Manager for PCA analysis and K-Means clustering of conformer ensembles.
    
    Handles the calculation of invariant features (EDM eigenvalues), 
    dimensionality reduction, and unsupervised clustering.
    """
    
    def __init__(
        self,
        logger: Logger,
        config: Optional[ClusteringConfig] = None
    ):
        """
        Initialize the manager.

        Args:
            logger (Logger): Logger instance.
            config (Optional[ClusteringConfig]): Clustering settings.
        """

        self.logger = logger
        self.config = config or ClusteringConfig()
        
        self.logger.debug(
            f"ClusteringManager initialized: "
            f"include_H={self.config.include_H}, "
        )
    
    # ====
    # Public API
    # ====
    
    def perform_pca(
        self,
        conformers: List[Conformer],
        n_clusters: Optional[int] = None,
        output_file: str = "pca_analysis.png",
        title: str = "PCA Analysis",
        include_legend: bool = True,
    ) -> Optional[PCAResult]:
        """
        Execute PCA pipeline and generate visualization.

        Args:
            conformers (List[Conformer]): List of conformers to analyze.
            n_clusters (Optional[int]): Target number of clusters (None for auto-detection).
            output_file (str): Output path for the plot.
            title (str): Title of the plot.
            include_legend (bool): Whether to draw the legend.

        Returns:
            Optional[PCAResult]: Result object if successful, None on failure.
        """
        
        
        # Filter active conformers
        active_confs = [c for c in conformers if c.active]
        
        self.logger.info(
            f"Starting PCA analysis on {len(active_confs)} conformers"
        )
        
        self.logger.info(
            f"  Hydrogen inclusion: {self.config.include_H}\n"
            f"  Requested clusters: {n_clusters or 'auto-detect'}\n"
            f"  Set cluster attribute: {self.config.set_cluster_attribute}"
        )
        
        try:
            # Execute PCA pipeline
            result = self._execute_pca_pipeline(
                conformers=active_confs,
                n_clusters=n_clusters,
                include_H=self.config.include_H
            )
            
            # Generate visualization
            self._create_visualization(
                result=result,
                output_file=Path(output_file),
                title=title,
                include_legend=include_legend
            )
            
            self.logger.info(
                f"✓ PCA completed: {result.n_clusters} clusters identified, "
                f"variance explained: {result.explained_variance[:2].sum()*100:.1f}% (PC1+PC2)"
            )
            
            return result
            
        except Exception as e:
            self.logger.error(f"✗ PCA failed: {e}")
            return None
    
    def reduce_by_clusters(
        self,
        conformers: List[Conformer],
        sort_by_energy: bool = False
    ) -> List[Conformer]:
        """
        Reduce the ensemble by keeping only one representative per cluster.

        Retains the lowest energy conformer for each identified cluster.

        Args:
            conformers (List[Conformer]): The clustered ensemble.
            sort_by_energy (bool): Whether to sort the result by energy.

        Returns:
            List[Conformer]: The reduced ensemble (modified in-place).
        """
        
        # Sort if requested
        ensemble = sorted(conformers) if sort_by_energy else conformers[:]
        
        # Track which clusters we've seen
        seen_clusters = []
        deactivated_count = 0
        
        for conf in ensemble:
            if not conf.active:
                continue
            
            # Check if we've seen this cluster
            if conf.cluster in seen_clusters:
                conf.active = False
                deactivated_count += 1
            else:
                seen_clusters.append(conf.cluster)
        
        active_count = len([c for c in ensemble if c.active])
        
        self.logger.info(
            f"Cluster-based reduction: "
            f"{active_count} conformers retained from {len(seen_clusters)} clusters "
            f"({deactivated_count} deactivated)"
        )
        
        return ensemble
    
    # ===
    # Private Methods - PCA Pipeline
    # ===
    
    def _execute_pca_pipeline(
        self,
        conformers: List[Conformer],
        n_clusters: Optional[int],
        include_H: bool
    ) -> PCAResult:
        """
        Execute complete PCA pipeline.
        
        Steps:
        1. Calculate distance matrices
        2. Compute eigenvalues (rotation/translation invariant features)
        3. Perform PCA transformation
        4. Determine optimal cluster number (if not specified)
        5. Perform K-Means clustering
        6. Optionally set cluster attribute on conformers
        
        Args:
            conformers: Active conformers
            n_clusters: Number of clusters (None = auto-detect)
            include_H: Include hydrogen in distance matrix
            
        Returns:
            PCAResult with all analysis data
        """
        
        # Step 1: Extract data
        self.logger.debug("Extracting conformer data...")
        geometries = np.array([c.last_geometry for c in conformers])
        colors = [c.color for c in conformers]
        numbers = [c.number for c in conformers]
        energies = np.array([c.energies.get_energy() for c in conformers], dtype=float)
        if np.isnan(energies).any():
            energies = None  # Energies missing: no energy surface
        else:
            energies = (energies - energies.min()) * EH_TO_KCAL  # Relative energies
        
        # Step 2: Calculate distance matrix eigenvalues
        self.logger.debug("Calculating distance matrices...")
        eigenvalues = self._calculate_distance_matrix_eigenvalues(
            geometries=geometries,
            atoms=conformers[0].atoms,
            include_H=include_H
        )
        
        # Step 3: PCA transformation
        self.logger.debug("Performing PCA transformation...")
        n_components = min(eigenvalues.shape[0], eigenvalues.shape[1])
        pca = PCA(n_components=n_components, random_state=self.config.random_state)
        pca_scores = pca.fit_transform(eigenvalues)
        
        # Step 4: Determine cluster number
        optimal_k = k_range = scores = None
        cluster_labels = None
        if n_clusters is None:
            self.logger.debug("Auto-detecting optimal cluster number...")
            optimal_k, k_range, scores, cluster_labels = self._silhouette_search(pca_scores)
            self.logger.info(f"  Optimal clusters detected: {optimal_k}")
            n_clusters = optimal_k

            floor = self._pruning_floor(len(pca_scores))
            if floor > optimal_k:
                n_clusters = floor
                cluster_labels = None
                self.logger.info(
                    f"  Pruning floor ({self.config.min_keep_fraction:.0%} of "
                    f"{len(pca_scores)} conformers) → using {n_clusters} clusters"
                )

        # Step 5: K-Means clustering (the search already gives the labels of its best k)
        if cluster_labels is None:
            self.logger.debug(f"Performing K-Means with {n_clusters} clusters...")
            kmeans = KMeans(
                n_clusters=n_clusters,
                n_init='auto',
                random_state=self.config.random_state
            )
            cluster_labels = kmeans.fit_predict(pca_scores)

        # Step 6: Set cluster attribute
        if self.config.set_cluster_attribute:
            for conf, cluster_id in zip(conformers, cluster_labels):
                conf.cluster = int(cluster_id)
        
        return PCAResult(
            scores=pca_scores,
            clusters=cluster_labels,
            colors=colors,
            numbers=numbers,
            energies=energies,
            explained_variance=pca.explained_variance_ratio_,
            n_clusters=n_clusters,
            original_features=eigenvalues,
            components=pca.components_,
            optimal_k=optimal_k,
            k_range=list(k_range) if k_range is not None else None,
            silhouette_scores=list(scores) if scores is not None else None,
        )
    
    def _calculate_distance_matrix_eigenvalues(
        self,
        geometries: np.ndarray,
        atoms: np.ndarray,
        include_H: bool
    ) -> np.ndarray:
        """
        Calculate eigenvalues of distance matrices.
        
        This provides rotation/translation invariant features for PCA.
        Vectorised over conformers, in chunks to bound memory.
        
        Args:
            geometries: (n_conformers, n_atoms, 3) array
            atoms: Array of atom symbols
            include_H: Include hydrogen atoms
            
        Returns:
            (n_conformers, n_eigenvalues) array
        """
        coords = np.asarray(geometries, dtype=float)
        if not include_H:
            coords = coords[:, np.asarray(atoms) != "H"]

        n_conf, n_at = coords.shape[:2]
        # ~64 MB of (chunk, n, n, 3) differences per chunk
        chunk = max(1, int(8e6 // max(1, n_at * n_at * 3)))
        out = np.empty((n_conf, n_at))
        for i in range(0, n_conf, chunk):
            c = coords[i:i + chunk]
            dist = np.linalg.norm(c[:, :, None, :] - c[:, None, :, :], axis=-1)
            # Real-symmetric matrix: eigvalsh gives real, ascending-sorted
            # eigenvalues, so feature columns correspond across conformers.
            # (np.linalg.eig returns them in arbitrary order.)
            out[i:i + chunk] = np.linalg.eigvalsh(dist)
        return out

    def _k_bounds(self, n_samples: int) -> tuple[int, int]:
        """Silhouette search bounds: [max(min_k, 2), min(frac*N, max_k, N-1)]."""
        lo = max(int(self.config.min_k), 2)
        hi = int(n_samples * self.config.max_k_fraction)
        if self.config.max_k is not None:
            hi = min(hi, int(self.config.max_k))
        hi = min(hi, n_samples - 1)
        return lo, hi

    def _pruning_floor(self, n_samples: int) -> int:
        """Minimum number of clusters in auto mode (0 if no floor)."""
        frac = self.config.min_keep_fraction or 0.0
        if frac <= 0:
            return 0
        return min(int(np.ceil(frac * n_samples)), n_samples)

    def _silhouette_search(self, features: np.ndarray) -> tuple:
        """
        Silhouette search over k, returning also the labels of the best k.

        The pairwise distance matrix is computed once and shared by every
        silhouette evaluation; the k values are evaluated in parallel when
        ``config.n_jobs != 1`` and the ensemble has at least
        ``_PARALLEL_MIN_SAMPLES`` conformers. Wide ranges are scanned on an adaptive grid
        (every k up to ``_DENSE_K``, then a coarse step refined around the
        best coarse k).

        Returns:
            tuple: (optimal_k, k_range, scores, labels)
        """
        n = len(features)
        lo, hi = self._k_bounds(n)

        if hi < lo:
            self.logger.warning(
                f"Cannot optimize clusters: max_k ({hi}) < min_k ({lo}). "
                f"Using k={lo}"
            )
            labels = KMeans(
                n_clusters=min(lo, n), n_init='auto', random_state=self.config.random_state
            ).fit_predict(features)
            return lo, [lo], [0.0], labels

        dist = cdist(features, features)
        np.fill_diagonal(dist, 0.0)

        results: Dict[int, tuple] = {}
        coarse, step = _candidate_ks(lo, hi)
        self._evaluate_ks(features, dist, coarse, results)
        if step > 1:
            best = max(results, key=lambda k: results[k][0])
            fine = [k for k in range(max(lo, best - step + 1), min(hi, best + step - 1) + 1)
                    if k not in results]
            self._evaluate_ks(features, dist, fine, results)

        k_range = sorted(results)
        scores = [results[k][0] for k in k_range]
        optimal_k = k_range[int(np.argmax(scores))]

        self.logger.debug(
            f"  Silhouette search: {len(k_range)} k values in [{lo}, {hi}], "
            f"best k={optimal_k} (score={results[optimal_k][0]:.3f})"
        )
        return optimal_k, k_range, scores, results[optimal_k][1]

    def _evaluate_ks(self, features, dist, ks, results) -> None:
        if not ks:
            return
        n_jobs = self.config.n_jobs or 1
        seed = self.config.random_state
        if n_jobs == 1 or len(features) < _PARALLEL_MIN_SAMPLES:
            out = [_score_k(features, dist, k, seed, False) for k in ks]
        else:
            out = Parallel(n_jobs=n_jobs, backend="loky")(
                delayed(_score_k)(features, dist, k, seed, True) for k in ks
            )
        for k, score, labels in out:
            results[k] = (score, labels)

    def find_optimal_clusters(self, features: np.ndarray) -> tuple:
        """
        Find optimal number of clusters using silhouette score.
        
        Args:
            features: Feature matrix (n_samples, n_features)
            
        Returns:
            tuple: (optimal_k, k_range, silhouette_scores)
        """
        optimal_k, k_range, scores, _ = self._silhouette_search(features)
        return optimal_k, k_range, scores
    
    # ===
    # Private Methods - Visualization
    # ===
    
    def _create_visualization(
        self,
        result: PCAResult,
        output_file: Path,
        title: str,
        include_legend: bool
    ) -> None:
        """
        Create PCA visualization with energy contours.
        
        Args:
            result: PCA analysis results
            output_file: Path for output image
            title: Plot title
            include_legend: Include conformer legend
        """
        
        self.logger.debug(f"Generating visualization: {output_file}")
        
        # Setup figure
        fig = plt.figure(figsize=(12, 9))
        
        # if include_legend:
        #     plt.subplots_adjust(bottom=0.3, right=0.65, left=0.10)
        
        rcParams.update({"figure.autolayout": True})
        
        # Extract data
        x = result.scores[:, 0]
        y = result.scores[:, 1]
        z = result.energies

        if z is None:
            # No energies: scatter only
            ax_main = fig.add_subplot(111)
        else:
            # Create grid: main plot + colorbar
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 0.1], hspace=0.3)
            ax_main = fig.add_subplot(gs[0])
            ax_colorbar = fig.add_subplot(gs[1])

            # Create interpolated grid for contours
            resolution = DEFAULT_RESOLUTION
            xi = np.linspace(x.min(), x.max(), resolution)
            yi = np.linspace(y.min(), y.max(), resolution)
            xi, yi = np.meshgrid(xi, yi)
            zi = griddata((x, y), z, (xi, yi), method='linear')

            # Plot energy surface
            im = ax_main.pcolormesh(
                xi, yi, zi,
                shading='auto',
                cmap='coolwarm',
                alpha=0.75
            )

            # Add contour lines
            ax_main.contour(
                xi, yi, zi,
                levels=6,
                colors='grey',
                linestyles='--',
                linewidths=0.5,
                alpha=0.6
            )

            # Colorbar
            cbar = plt.colorbar(im, cax=ax_colorbar, orientation='horizontal')
            cbar.set_label('Relative Energy [kcal/mol]', fontsize=10)
        
        # Plot conformers with cluster-specific markers
        for idx in range(len(x)):
            marker = self._get_cluster_marker(result.clusters[idx])
            ax_main.scatter(
                x[idx], y[idx],
                c=result.colors[idx],
                marker=marker,
                s=100,
                linewidths=0.5,
                label=f"Conf {result.numbers[idx]} (C{result.clusters[idx]})",
                zorder=10
            )
        
        # Axes and grid
        ax_main.axhline(0, color='#353535', linestyle='--', alpha=0.3, linewidth=1)
        ax_main.axvline(0, color='#353535', linestyle='--', alpha=0.3, linewidth=1)
        ax_main.grid(alpha=0.2, linestyle=':')
        
        # Labels
        variance_pc1 = result.explained_variance[0] * 100
        variance_pc2 = result.explained_variance[1] * 100
        
        ax_main.set_xlabel(
            f'Principal Component 1 ({variance_pc1:.1f}% variance)',
            fontsize=11
        )
        ax_main.set_ylabel(
            f'Principal Component 2 ({variance_pc2:.1f}% variance)',
            fontsize=11
        )
        ax_main.set_title(
            f'{title}\n({len(result.numbers)} conformers, '
            f'{result.n_clusters} clusters)',
            fontsize=13,
            fontweight='bold'
        )
        
        # Legend
        if include_legend:
            n_items = len(result.numbers)
            items_per_col = 30
            n_cols = max(1, (n_items // items_per_col) + (1 if n_items % items_per_col else 0))
            font_size = 9 if n_cols < 4 else 7

            ax_main.legend(
                loc='upper left',
                bbox_to_anchor=(1.02, 1.0),
                fancybox=True,
                shadow=True,
                ncol=n_cols,
                fontsize=font_size,
                markerscale=0.8,
                borderaxespad=0.
            )
        else:
            plt.tight_layout()
        
        # Save
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        self.logger.debug(f"✓ Visualization saved: {output_file}")
    
    @staticmethod
    def _get_cluster_marker(cluster_id: int) -> str:
        """Get marker style for cluster ID"""
        return MARKERS[cluster_id % len(MARKERS)]

