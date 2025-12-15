
from ensemble_analyzer._logger.logger import Logger
from ensemble_analyzer._conformer.conformer import Conformer

from pathlib import Path

from ensemble_analyzer.constants import *

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
from scipy.spatial import distance_matrix

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
            self.logger.error(f"{self.logger.FAIL} PCA failed: {e}")
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
        seen_clusters = set()
        deactivated_count = 0
        
        for conf in ensemble:
            if not conf.active:
                continue
            
            # Check if we've seen this cluster
            if conf.cluster in seen_clusters:
                conf.active = False
                deactivated_count += 1
            else:
                seen_clusters.add(conf.cluster)
        
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
        energies = np.array([c.energies.get_energy() for c in conformers])
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
        if n_clusters is None:
            self.logger.debug("Auto-detecting optimal cluster number...")
            n_clusters = self._find_optimal_clusters(pca_scores)
            self.logger.info(f"  Optimal clusters detected: {n_clusters}")
        
        # Step 5: K-Means clustering
        self.logger.debug(f"Performing K-Means with {n_clusters} clusters...")
        kmeans = KMeans(
            n_clusters=n_clusters,
            n_init='auto',
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
            n_clusters=n_clusters
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
        
        Args:
            geometries: (n_conformers, n_atoms, 3) array
            atoms: Array of atom symbols
            include_H: Include hydrogen atoms
            
        Returns:
            (n_conformers, n_eigenvalues) array
        """
        
        eigenvalues_list = []
        
        for geom in geometries:
            # Filter atoms if needed
            if include_H:
                coords = geom
            else:
                mask = atoms != "H"
                coords = geom[mask]
            
            # Calculate distance matrix
            dist_mat = distance_matrix(coords, coords)
            
            # Compute eigenvalues
            eigenvals, _ = np.linalg.eig(dist_mat)
            eigenvalues_list.append(eigenvals.real)  # Take real part
        
        return np.array(eigenvalues_list)
    
    def _find_optimal_clusters(self, features: np.ndarray) -> int:
        """
        Find optimal number of clusters using silhouette score.
        
        Args:
            features: Feature matrix (n_samples, n_features)
            
        Returns:
            Optimal number of clusters
        """
        
        min_k = max(int(len(features)*.1), 2) # Set as 10% of the ensemble length or 2
        max_k = int(len(features)*.8) # Set as the 80% of the ensemble length
        
        if max_k < min_k:
            self.logger.warning(
                f"Cannot optimize clusters: max_k ({max_k}) < min_k ({min_k}). "
                f"Using min_k={min_k}"
            )
            return min_k
        
        k_range = range(min_k, max_k + 1)
        silhouette_scores = []
        
        for k in k_range:
            kmeans = KMeans(
                n_clusters=k,
                n_init='auto',
            )
            labels = kmeans.fit_predict(features)
            score = silhouette_score(features, labels)
            silhouette_scores.append(score)
        
        optimal_k = k_range[np.argmax(silhouette_scores)]
        best_score = max(silhouette_scores)
        
        self.logger.debug(
            f"  Silhouette scores: k={optimal_k} (score={best_score:.3f})"
        )
        
        return optimal_k
    
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
        
        # Create grid: main plot + colorbar
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 0.1], hspace=0.3)
        ax_main = fig.add_subplot(gs[0])
        ax_colorbar = fig.add_subplot(gs[1])
        
        # Extract data
        x = result.scores[:, 0]
        y = result.scores[:, 1]
        z = result.energies
        
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

