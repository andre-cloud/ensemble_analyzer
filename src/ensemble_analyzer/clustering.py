
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
) -> None:
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

    import argparse
    from unittest.mock import Mock
    from ensemble_analyzer.ensemble_io import read_ensemble, save_snapshot
    
    parser = argparse.ArgumentParser(
        description='Perform PCA analysis and clustering on conformer ensemble',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic PCA with 5 clusters
  python clustering.py ensemble.xyz -nc 5
  
  # Auto-detect optimal clusters, exclude hydrogen
  python clustering.py ensemble.xyz --no-H
  
  # Save reduced ensemble
  python clustering.py ensemble.xyz -nc 10 --save-reduced
  
  # Custom output and title
  python clustering.py ensemble.xyz -o my_pca.png --title "Drug Conformers"
        """
    )
    
    parser.add_argument(
        'file',
        help='Input ensemble file (XYZ format)'
    )
    parser.add_argument(
        '-nc', '--ncluster',
        type=int,
        default=None,
        help='Number of clusters (default: auto-detect using silhouette score)'
    )
    parser.add_argument(
        '--no-H',
        action='store_false',
        dest='include_H',
        help='Exclude hydrogen atoms from distance matrix calculation'
    )
    parser.add_argument(
        '--no-legend',
        action='store_false',
        dest='legend',
        help='Exclude conformer legend from plot'
    )
    parser.add_argument(
        '--title',
        default='PCA Cluster Analysis',
        help='Plot title (default: "PCA Cluster Analysis")'
    )
    parser.add_argument(
        '-o', '--output',
        default='cluster.png',
        help='Output image filename (default: cluster.png)'
    )
    parser.add_argument(
        '--save-reduced',
        action='store_true',
        help='Save cluster-reduced ensemble to clustered.xyz'
    )
    parser.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='Output image DPI (default: 300)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    args = parser.parse_args()
    
    # Setup mock logger with print statements
    class CLILogger:
        def __init__(self, verbose=False):
            self.verbose = verbose
        
        def info(self, msg):
            print(f"ℹ️ {msg}")
        
        def debug(self, msg):
            if self.verbose:
                print(f"🔍 {msg}")
        
        def warning(self, msg):
            print(f"⚠️  {msg}")
        
        def error(self, msg):
            print(f"❌ {msg}")
    
    logger = CLILogger(verbose=args.verbose)
    
    # Load ensemble
    print(f"\n{'='*60}")
    print(f"PCA CLUSTERING ANALYSIS")
    print(f"{'='*60}\n")
    
    logger.info(f"Loading ensemble from {args.file}...")
    try:
        ensemble = read_ensemble(args.file, logger, raw=True)
        logger.info(f"✓ Loaded {len(ensemble)} conformers")
    except Exception as e:
        logger.error(f"Failed to load ensemble: {e}")
        exit(1)
    
    # Perform PCA
    config = ClusteringConfig(
        n_clusters=args.ncluster,
        include_H=args.include_H,
        set_cluster_attribute=True
    )
    
    manager = ClusteringManager(logger=logger, config=config)
    
    result = manager.perform_pca(
        conformers=ensemble,
        n_clusters=args.ncluster,
        output_file=args.output,
        title=args.title,
        include_legend=args.legend
    )
    
    if result:
        print(f"\n{'='*60}")
        print(f"RESULTS")
        print(f"{'='*60}")
        print(f"✓ Clusters identified: {result.n_clusters}")
        print(f"✓ Variance explained (PC1): {result.explained_variance[0]*100:.1f}%")
        print(f"✓ Variance explained (PC2): {result.explained_variance[1]*100:.1f}%")
        print(f"✓ Total variance (PC1+PC2): {result.explained_variance[:2].sum()*100:.1f}%")
        print(f"✓ Output saved: {args.output}")
        
        # Save reduced ensemble if requested
        if args.save_reduced:
            logger.info("\nReducing ensemble by clusters...")
            reduced = manager.reduce_by_clusters(ensemble)
            save_snapshot("clustered.xyz", reduced, logger)
            active = len([c for c in reduced if c.active])
            logger.info(f"✓ Reduced ensemble: {active} conformers → clustered.xyz")
        
        print(f"{'='*60}\n")
        
    else:
        logger.error("PCA analysis failed or was skipped")
        exit(1)