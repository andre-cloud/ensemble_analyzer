
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
    """    
    
    Args:
        confs: List of conformers
        ncluster: Number of clusters (None = auto-detect)
        fname: Output filename
        title: Plot title
        log: Logger instance
        set_: Set cluster attribute on conformers
        include_H: Include hydrogen in distance matrix
        legend: Include legend in plot
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

def validate_possible_PCA(ensemble: List[Conformer], logger: Logger, n_clusters: Optional[Union[int, bool]]):

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

