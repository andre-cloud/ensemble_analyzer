
from ensemble_analyzer._conformer.conformer import Conformer
from ensemble_analyzer._logger.logger import Logger

from ensemble_analyzer._clustering.cluster_config import ClusteringConfig
from ensemble_analyzer._clustering.cluster_manager import ClusteringManager

from ensemble_analyzer.constants import * 

from typing import List, Optional, Union



def execute_PCA(
    confs: List[Conformer],
    ncluster: Optional[int],
    fname: str,
    title: str,
    log: Logger,
    set_: bool = True,
    include_H: bool = True,
    legend: bool = True
) -> bool:
    """Execute PCA-based clustering on the conformer ensemble.

    Args:
        confs (List[Conformer]): List of conformers.
        ncluster (Optional[int]): Number of clusters (None = auto-detect).
        fname (str): Output filename for the plot.
        title (str): Plot title.
        log (Logger): Logger instance.
        set_ (bool): Set cluster attribute on conformers.
        include_H (bool): Include hydrogen in distance matrix.
        legend (bool): Include legend in plot.

    Returns:
        bool: True if PCA was performed successfully, False otherwise.
    """
    
    config = ClusteringConfig(
        n_clusters=ncluster,
        include_H=include_H,
        set_cluster_attribute=set_
    )
    
    manager = ClusteringManager(logger=log, config=config)
    

    if validate_possible_PCA(ensemble=confs, logger=log, n_clusters=ncluster):
        manager.perform_pca(
            conformers=confs,
            n_clusters=ncluster,
            output_file=fname,
            title=title,
            include_legend=legend,
        )
        return True
    return False

def validate_possible_PCA(ensemble: List[Conformer], logger: Logger, n_clusters: Optional[Union[int, bool]]) -> bool:
    """Check preconditions for PCA execution.

    Verifies that enough active conformers exist and that the number of
    clusters does not exceed the ensemble size.

    Args:
        ensemble (List[Conformer]): List of conformers.
        logger (Logger): Logger instance.
        n_clusters (Optional[Union[int, bool]]): Requested number of clusters.

    Returns:
        bool: True if PCA can proceed, False otherwise.
    """

    ensemble = [conf for conf in ensemble if conf.active]
    if len(ensemble) < MIN_CONFORMERS_FOR_PCA:
        logger.warning(
            f"PCA skipped: only {len(ensemble)} active conformers "
            f"(minimum {MIN_CONFORMERS_FOR_PCA} required)"
        )
        return False
    
    if n_clusters and len(ensemble) < n_clusters:
        logger.warning(
            f"PCA skipped: n_clusters ({n_clusters}) >= "
            f"n_conformers ({len(ensemble)})"
        )
        return False

    return True


def get_ensemble(
    confs: List[Conformer],
    log : Logger,
    sort: bool = False
) -> List[Conformer]:
    """
    Get pruned ensemble
    
    Args:
        confs: Conformer ensemble
        log: Logger instance
        sort: Sort by energy
        
    Returns:
        Reduced ensemble
    """
    
    manager = ClusteringManager(logger=log)
    return manager.reduce_by_clusters(confs, sort_by_energy=sort)


# ===
# CLI for Standalone Usage
# ===

if __name__ == "__main__":
    """
    Standalone CLI for PCA analysis.
    """

