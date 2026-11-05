from dataclasses import dataclass

from typing import List, Optional
import numpy as np


@dataclass
class PCAResult:
    """
    Container for PCA and Clustering analysis results.
    """
    scores: np.ndarray              # PCA coordinates (n_samples, n_components)
    clusters: np.ndarray            # Cluster labels (n_samples,)
    colors: List[str]               # Hex colors for each point
    numbers: List[int]              # Conformer IDs
    energies: Optional[np.ndarray]  # Relative energies array (None if unavailable)
    explained_variance: np.ndarray  # Variance ratio per component
    n_clusters: Optional[int] = None # Number of clusters used for the labels
    original_features: Optional[np.ndarray] = None  # Added for before/after plot
    components: Optional[np.ndarray] = None         # Added for loading plot
    optimal_k: Optional[int] = None                 # k picked by the silhouette search (auto mode)
    k_range: Optional[List[int]] = None             # k values tried by the search
    silhouette_scores: Optional[List[float]] = None # silhouette score for each k in k_range


@dataclass
class ClusteringConfig:
    """
    Configuration parameters for clustering operations.
    """
    n_clusters: Optional[int] = None    # Target clusters (None = auto)
    include_H: bool = False             # Include Hydrogens in distance matrix
    set_cluster_attribute: bool = True  # Write cluster ID to conformer objects
    min_k: int = 2                      # Min clusters for silhouette search
    max_k_fraction: float = 0.8         # Max clusters for silhouette search, as a fraction of N
    max_k: Optional[int] = None         # Optional absolute cap on the silhouette search
    min_keep_fraction: float = 0.0      # Auto mode only: use at least ceil(fraction*N) clusters
                                        # (pruning floor, so that pruning never keeps too few)
    n_jobs: int = 1                     # Parallel workers for the k search
    random_state: int = 42              # Seed for reproducibility