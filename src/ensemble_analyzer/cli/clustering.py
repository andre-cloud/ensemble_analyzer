import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from unittest.mock import Mock

from ensemble_analyzer.ensemble_io import read_ensemble, save_snapshot
from ensemble_analyzer.clustering import *
from ensemble_analyzer._clustering.cluster_config import ClusteringConfig
from ensemble_analyzer._clustering.cluster_manager import ClusteringManager

def parse_args():
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

    parser.add_argument('file', help='Input ensemble file (XYZ format)')
    parser.add_argument('-nc', '--ncluster', type=int, default=None,
                        help='Number of clusters (default: auto-detect using silhouette score)')
    parser.add_argument('--no-H', action='store_false', dest='include_H',
                        help='Exclude hydrogen atoms from distance matrix calculation')
    parser.add_argument('--no-legend', action='store_false', dest='legend',
                        help='Exclude conformer legend from plot')
    parser.add_argument('--title', default='PCA Cluster Analysis',
                        help='Plot title (default: "PCA Cluster Analysis")')
    parser.add_argument('-o', '--output', default='cluster.png',
                        help='Output image filename (default: cluster.png)')
    parser.add_argument('--write', action='store_true',
                        help='Save cluster-reduced ensemble to clustered.xyz')
    parser.add_argument('--dpi', type=int, default=300,
                        help='Output image DPI (default: 300)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')

    args = parser.parse_args()
    return args

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

def plot_component_analysis(result, output_base, logger):
    """
    Generate loading plot showing feature contributions to principal components.
    """
    out_dir = os.path.dirname(output_base) or "."
    base_name = os.path.splitext(os.path.basename(output_base))[0]
    
    logger.info("Generating component loading analysis...")
    
    if getattr(result, 'components', None) is None:
        logger.warning("No component loadings available, skipping loading plot")
        return
    
    components = result.components[:2]  # First 2 PCs
    n_features = components.shape[1]
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    for idx, (ax, pc) in enumerate(zip(axes, components)):
        pc_variance = result.explained_variance[idx] * 100
        feature_indices = np.arange(n_features)
        colors = ['#2E86AB' if x >= 0 else '#A23B72' for x in pc]
        
        ax.bar(feature_indices, pc, color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
        ax.set_xlabel('Feature Index', fontsize=10)
        ax.set_ylabel('Loading Value', fontsize=10)
        ax.set_title(f'PC{idx+1} Loadings ({pc_variance:.1f}% variance)', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Highlight top contributors
        top_indices = np.argsort(np.abs(pc))[-5:]
        for top_idx in top_indices:
            ax.text(top_idx, pc[top_idx], f'{top_idx}', 
                   ha='center', va='bottom' if pc[top_idx] > 0 else 'top',
                   fontsize=8, fontweight='bold')
    
    plt.tight_layout()
    loading_file = os.path.join(out_dir, f"{base_name}_loadings.png")
    plt.savefig(loading_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"✓ Component loadings saved: {loading_file}")

def plot_before_after_pca(result, output_base, logger):
    """
    Side-by-side comparison of data before and after PCA transformation,
    with PCA axes drawn in the original space.
    """
    out_dir = os.path.dirname(output_base) or "."
    base_name = os.path.splitext(os.path.basename(output_base))[0]
    
    logger.info("Generating before/after PCA comparison...")
    
    if getattr(result, 'original_features', None) is None:
        logger.warning("Original features not available, skipping before/after plot")
        return
    
    original = result.original_features
    transformed = result.scores
    labels = result.clusters
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    unique_labels = np.unique(labels)
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_labels)))
    
    # Plot 1: Original space (first 2 features)
    ax1 = axes[0]
    for idx, label in enumerate(unique_labels):
        mask = labels == label
        ax1.scatter(original[mask, 0], original[mask, 1], 
                   c=[colors[idx]],
                   alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
    
    # Draw PCA axes in the original space if components are available
    if getattr(result, 'components', None) is not None:
        # Calculate the center of the original data (since PCA is centered)
        mean_0 = np.mean(original[:, 0])
        mean_1 = np.mean(original[:, 1])
        
        # Extract the projection of the first two PCs onto the first two original features
        v1_0, v1_1 = result.components[0, 0], result.components[0, 1]
        v2_0, v2_1 = result.components[1, 0], result.components[1, 1]
        
        # Scale vectors based on the standard deviation of the scores
        scale1 = np.std(transformed[:, 0]) * 2
        scale2 = np.std(transformed[:, 1]) * 2
        
        # Plot PC1 Axis (Red)
        ax1.plot([mean_0 - v1_0*scale1, mean_0 + v1_0*scale1], 
                 [mean_1 - v1_1*scale1, mean_1 + v1_1*scale1], 
                 color='red', linewidth=2, linestyle='--')
        ax1.text(mean_0 + v1_0*scale1, mean_1 + v1_1*scale1, 'PC1', color='red', fontsize=12, fontweight='bold')
        
        # Plot PC2 Axis (Green)
        ax1.plot([mean_0 - v2_0*scale2, mean_0 + v2_0*scale2], 
                 [mean_1 - v2_1*scale2, mean_1 + v2_1*scale2], 
                 color='green', linewidth=2, linestyle='--')
        ax1.text(mean_0 + v2_0*scale2, mean_1 + v2_1*scale2, 'PC2', color='green', fontsize=12, fontweight='bold')

    ax1.set_xlabel('Feature 1 (Original Space)', fontsize=11)
    ax1.set_ylabel('Feature 2 (Original Space)', fontsize=11)
    ax1.set_title('Original Feature Space with PC Axes', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: PCA space (PC1 vs PC2)
    ax2 = axes[1]
    for idx, label in enumerate(unique_labels):
        mask = labels == label
        ax2.scatter(transformed[mask, 0], transformed[mask, 1],
                   c=[colors[idx]],
                   alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
    
    var1, var2 = result.explained_variance[0] * 100, result.explained_variance[1] * 100
    ax2.set_xlabel(f'PC1 ({var1:.1f}% variance)', fontsize=11)
    ax2.set_ylabel(f'PC2 ({var2:.1f}% variance)', fontsize=11)
    ax2.set_title(f'PCA-Transformed Space', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    comparison_file = os.path.join(out_dir, f"{base_name}_before_after.png")
    plt.savefig(comparison_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"✓ Before/After comparison saved: {comparison_file}")

def plot_clustering_metrics(result, output_base, logger):
    """
    Generate evaluation plots: Scree plot with cumulative variance and Silhouette scores.
    """
    out_dir = os.path.dirname(output_base) or "."
    base_name = os.path.splitext(os.path.basename(output_base))[0]
    
    logger.info("Generating clustering metric plots...")
    
    variance = result.explained_variance * 100
    cumulative_variance = np.cumsum(variance)
    n_components = len(variance)
    
    fig, ax1 = plt.subplots(figsize=(10, 6))
    
    x_pos = np.arange(1, n_components + 1)
    bars = ax1.bar(x_pos, variance, alpha=0.7, color='#3498db', 
                   edgecolor='black', linewidth=0.8, label='Individual Variance')
    
    ax1.set_xlabel('Principal Component', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Explained Variance (%)', fontsize=12, fontweight='bold', color='#3498db')
    ax1.tick_params(axis='y', labelcolor='#3498db')
    # ax1.set_xticks(x_pos)
    ax1.grid(True, alpha=0.3, axis='y', linestyle='--')
    
    ax2 = ax1.twinx()
    line = ax2.plot(x_pos, cumulative_variance, color='#e74c3c', marker='o', 
                    linewidth=2.5, markersize=8, label='Cumulative Variance')
    
    ax2.set_ylabel('Cumulative Variance (%)', fontsize=12, fontweight='bold', color='#e74c3c')
    ax2.tick_params(axis='y', labelcolor='#e74c3c')
    ax2.set_ylim([0, 105])
    
    ax2.axhline(y=90, color='gray', linestyle='--', linewidth=1.5, alpha=0.7, label='90% threshold')
    
    for i, (bar, val) in enumerate(zip(bars, variance)):
        if i < 5:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{val:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='lower right', fontsize=10)
    
    plt.title('Scree Plot: Explained Variance per Component', fontsize=13, fontweight='bold', pad=15)
    plt.tight_layout()
    
    scree_file = os.path.join(out_dir, f"{base_name}_scree.png")
    plt.savefig(scree_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"✓ Scree plot saved: {scree_file}")

    features = result.scores
    min_k = max(int(len(features) * 0.1), 2)
    max_k = int(len(features) * 0.8)
    
    if max_k > min_k:
        k_range = range(min_k, max_k + 1)
        scores = []
        for k in k_range:
            kmeans = KMeans(n_clusters=k, n_init='auto', random_state=500)
            labels = kmeans.fit_predict(features)
            scores.append(silhouette_score(features, labels))
        
        plt.figure(figsize=(8, 5))
        plt.plot(k_range, scores, marker='s', color='orange')
        plt.axvline(
            x=result.n_clusters, color='red', linestyle='--', 
            label=f'Selected k={result.n_clusters}'
        )
        plt.title('Silhouette Score vs Number of Clusters')
        plt.xlabel('Number of Clusters (k)')
        plt.ylabel('Silhouette Score')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        sil_file = os.path.join(out_dir, f"{base_name}_silhouette.png")
        plt.savefig(sil_file, dpi=300, bbox_inches='tight')
        plt.close()
        logger.info(f"✓ Silhouette plot saved: {sil_file}")
    else:
        logger.warning("Not enough data points to generate Silhouette plot")


def plot_3d_original_space(result, output_base, logger):
    """
    Generate a 3D scatter plot of the first 3 original features to illustrate
    the complexity of the data before PCA transformation.
    """
    out_dir = os.path.dirname(output_base) or "."
    base_name = os.path.splitext(os.path.basename(output_base))[0]
    
    logger.info("Generating 3D original space plot...")
    
    if getattr(result, 'original_features', None) is None or result.original_features.shape[1] < 3:
        logger.warning("Not enough original features available for 3D plot, skipping")
        return
    
    original = result.original_features
    labels = result.clusters
    
    fig, ax = plt.subplots(figsize=(5, 6))
 
    
    unique_labels = np.unique(labels)
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_labels)))
    
    for idx, label in enumerate(unique_labels):
        mask = labels == label
        ax.scatter(
            original[mask, 0], 
            original[mask, 1], 
            c=[colors[idx]], 
            # label=f'Cluster {label}',
            alpha=0.7, 
            s=40, 
            edgecolors='black', 
            linewidth=0.5
        )
    
    ax.set_xlabel('Feature 1', fontsize=10)
    ax.set_ylabel('Feature 2', fontsize=10)
    # ax.set_zlabel('Feature 3', fontsize=10)
    ax.set_title('Original Features (Pre-PCA)', fontsize=12, fontweight='bold')
    ax.grid()
    # Imposta un angolo di visualizzazione iniziale
    # ax.view_init(elev=20, azim=45)
    
    plt.tight_layout()
    plot_file = os.path.join(out_dir, f"{base_name}_3d_original.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"✓ 3D Original Space plot saved: {plot_file}")

def main():
    args = parse_args()
    logger = CLILogger(verbose=args.verbose)

    print(f"\n{'='*60}")
    print(f"PCA CLUSTERING ANALYSIS")
    print(f"{'='*60}\n")

    logger.info(f"Loading ensemble from {args.file}...")
    try:
        ensemble = read_ensemble(args.file, raw=True)
        logger.info(f"✓ Loaded {len(ensemble)} conformers")
    except TypeError:
        try:
            ensemble = read_ensemble(args.file, logger, raw=True)
            logger.info(f"✓ Loaded {len(ensemble)} conformers")
        except Exception as e:
            logger.error(f"Failed to load ensemble: {e}")
            exit(1)
    except Exception as e:
        logger.error(f"Failed to load ensemble: {e}")
        exit(1)

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
        plot_clustering_metrics(result, args.output, logger)
        plot_component_analysis(result, args.output, logger)
        plot_before_after_pca(result, args.output, logger)

        print(f"\n{'='*60}")
        print(f"RESULTS")
        print(f"{'='*60}")
        print(f"✓ Clusters identified: {result.n_clusters}")
        print(f"✓ Variance explained (PC1): {result.explained_variance[0]*100:.1f}%")
        print(f"✓ Variance explained (PC2): {result.explained_variance[1]*100:.1f}%")
        print(f"✓ Total variance (PC1+PC2): {result.explained_variance[:2].sum()*100:.1f}%")
        print(f"✓ Output saved: {args.output}")
        
        if getattr(result, 'components', None) is not None:
            print(f"\n{'='*60}")
            print(f"COMPONENT ANALYSIS")
            print(f"{'='*60}")
            for i in range(min(3, len(result.components))):
                pc = result.components[i]
                top_features = np.argsort(np.abs(pc))[-5:][::-1]
                print(f"PC{i+1} top features: {list(top_features)} (variance: {result.explained_variance[i]*100:.1f}%)")
            print(f"{'='*60}")

        if args.write:
            logger.info("\nReducing ensemble by clusters...")
            print([i.cluster for i in ensemble])
            reduced = manager.reduce_by_clusters(ensemble, True)
            save_snapshot("clustered.xyz", reduced, logger)
            active = len([c for c in reduced if c.active])
            logger.info(f"✓ Reduced ensemble: {active} conformers → clustered.xyz")
    
        print(f"{'='*60}\n")
    else:
        logger.error("PCA analysis failed or was skipped")
        exit(1)

if __name__ == '__main__':
    main()