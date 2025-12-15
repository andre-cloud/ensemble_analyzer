from dataclasses import dataclass

from typing import List, Optional
import numpy as np


@dataclass
class PCAResult:
    """Results from PCA analysis"""
    scores              : np.ndarray                # PCA-transformed coordinates (n_conformers, n_components)
    clusters            : np.ndarray                # Cluster assignments (n_conformers,)
    colors              : List[str]                 # Conformer colors
    numbers             : List[int]                 # Conformer IDs
    energies            : np.ndarray                # Relative energies [Eh]
    explained_variance  : np.ndarray                # Variance explained by each component
    n_clusters          : Optional[int]     = None  # Number of clusters used


@dataclass
class ClusteringConfig:
    """Configuration for clustering operations"""
    n_clusters              : Optional[int]     = None  # If None, auto-detect optimal
    include_H               : bool              = True  # Include hydrogen in distance matrix
    set_cluster_attribute   : bool              = True  # Set cluster ID on Conformer objects
    min_k                   : int               = 2     # Minimum clusters for silhouette search
    max_k                   : int               = 30    # Maximum clusters for silhouette search
    random_state            : int               = 42    # Fix value for reproducibility of the random initiation of cluster points
