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
    energies: np.ndarray            # Relative energies array
    explained_variance: np.ndarray  # Variance ratio per component
    n_clusters: Optional[int] = None # Number of clusters found/used


@dataclass
class ClusteringConfig:
    """
    Configuration parameters for clustering operations.
    """
    n_clusters: Optional[int] = None    # Target clusters (None = auto)
    include_H: bool = True              # Include Hydrogens in distance matrix
    set_cluster_attribute: bool = True  # Write cluster ID to conformer objects
    min_k: int = 2                      # Min clusters for silhouette search
    max_k: int = 30                     # Max clusters for silhouette search
    random_state: int = 42              # Seed for reproducibility