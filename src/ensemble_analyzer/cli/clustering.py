import argparse
import logging
import os
from typing import Any

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402  (after the backend is chosen)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Perform PCA analysis and clustering on conformer ensemble',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
# Basic PCA with 5 clusters
python clustering.py ensemble.xyz -nc 5
# Auto-detect optimal clusters, include hydrogen
python clustering.py ensemble.xyz --include-H
# Save reduced ensemble
python clustering.py ensemble.xyz -nc 10 --save-reduced
# Custom output and title
python clustering.py ensemble.xyz -o my_pca.png --title "Drug Conformers"
        """
    )

    parser.add_argument('file', help='Input ensemble file (XYZ format)')
    parser.add_argument('-nc', '--ncluster', type=int, default=None,
                        help='Number of clusters (default: auto-detect using silhouette score)')
    parser.add_argument('--include-H', action='store_true', dest='include_H',
                        help='Include hydrogen atoms in distance matrix calculation (default: excluded)')
    parser.add_argument('--no-legend', action='store_false', dest='legend',
                        help='Exclude conformer legend from plot')
    parser.add_argument('--title', default='PCA Cluster Analysis',
                        help='Plot title (default: "PCA Cluster Analysis")')
    parser.add_argument('-o', '--output', default='cluster.png',
                        help='Output image filename (default: cluster.png)')
    parser.add_argument('--write', action='store_true',
                        help='Save cluster-reduced ensemble to clustered.xyz')
    parser.add_argument('-cpu', '--cpu', type=int, default=1, dest='cpu',
                        help='Parallel workers for the automatic cluster search (default: 1)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')

    args = parser.parse_args()
    return args


def plot_component_analysis(result: Any, out_dir: str, base_name: str, logger: logging.Logger) -> None:

    logger.info("Generating component loading analysis...")

    components = result.components[:2]
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


def plot_before_after_pca(result: Any, out_dir: str, base_name: str, logger: logging.Logger) -> None:

    logger.info("Generating before/after PCA comparison...")

    original = result.original_features
    transformed = result.scores
    labels = result.clusters

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    unique_labels = np.unique(labels)
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_labels)))

    ax1 = axes[0]
    for idx, label in enumerate(unique_labels):
        mask = labels == label
        ax1.scatter(original[mask, 0], original[mask, 1],
                    c=[colors[idx]],
                    alpha=0.6, s=50, edgecolors='black', linewidth=0.5)

    mean_0 = np.mean(original[:, 0])
    mean_1 = np.mean(original[:, 1])

    v1_0, v1_1 = result.components[0, 0], result.components[0, 1]
    v2_0, v2_1 = result.components[1, 0], result.components[1, 1]

    scale1 = np.std(transformed[:, 0]) * 2
    scale2 = np.std(transformed[:, 1]) * 2

    ax1.plot([mean_0 - v1_0 * scale1, mean_0 + v1_0 * scale1],
             [mean_1 - v1_1 * scale1, mean_1 + v1_1 * scale1],
             color='red', linewidth=2, linestyle='--')
    ax1.text(mean_0 + v1_0 * scale1, mean_1 + v1_1 * scale1,
             'PC1', color='red', fontsize=12, fontweight='bold')

    ax1.plot([mean_0 - v2_0 * scale2, mean_0 + v2_0 * scale2],
             [mean_1 - v2_1 * scale2, mean_1 + v2_1 * scale2],
             color='green', linewidth=2, linestyle='--')
    ax1.text(mean_0 + v2_0 * scale2, mean_1 + v2_1 * scale2,
             'PC2', color='green', fontsize=12, fontweight='bold')

    ax1.set_xlabel('Feature 1 (Original Space)', fontsize=11)
    ax1.set_ylabel('Feature 2 (Original Space)', fontsize=11)
    ax1.set_title('Original Feature Space with PC Axes', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    ax2 = axes[1]
    for idx, label in enumerate(unique_labels):
        mask = labels == label
        ax2.scatter(transformed[mask, 0], transformed[mask, 1],
                    c=[colors[idx]],
                    alpha=0.6, s=50, edgecolors='black', linewidth=0.5)

    var1, var2 = result.explained_variance[0] * 100, result.explained_variance[1] * 100
    ax2.set_xlabel(f'PC1 ({var1:.1f}% variance)', fontsize=11)
    ax2.set_ylabel(f'PC2 ({var2:.1f}% variance)', fontsize=11)
    ax2.set_title('PCA-Transformed Space', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    comparison_file = os.path.join(out_dir, f"{base_name}_before_after.png")
    plt.savefig(comparison_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"✓ Before/After comparison saved: {comparison_file}")


def plot_clustering_metrics(result: Any, out_dir: str, base_name: str, logger: logging.Logger, explicit_k: int | None = None) -> None:
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
    ax1.grid(True, alpha=0.3, axis='y', linestyle='--')

    ax2 = ax1.twinx()
    ax2.plot(x_pos, cumulative_variance, color='#e74c3c', marker='o',
                    linewidth=2.5, markersize=8, label='Cumulative Variance')

    ax2.set_ylabel('Cumulative Variance (%)', fontsize=12, fontweight='bold', color='#e74c3c')
    ax2.tick_params(axis='y', labelcolor='#e74c3c')
    ax2.set_ylim([0, 105])

    ax2.axhline(y=90, color='gray', linestyle='--', linewidth=1.5, alpha=0.7, label='90% threshold')

    for i, (bar, val) in enumerate(zip(bars, variance)):
        if i < 5:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width() / 2., height,
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

    if explicit_k is None:
        # reuse the search already run by perform_pca
        k_range, scores = result.k_range, result.silhouette_scores
    else:
        k_range, scores = range(explicit_k, explicit_k + 1), []

    if len(k_range) > 1:
        plt.figure(figsize=(8, 5))
        plt.plot(k_range, scores, marker='s', color='orange')
        selected_k = getattr(result, 'optimal_k', None) or result.n_clusters
        plt.axvline(
            x=selected_k, color='red', linestyle='--',
            label=f'Selected k={selected_k}',
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


def main() -> None:
    args = parse_args()

    from ensemble_analyzer.ensemble_io import read_ensemble, save_snapshot
    from ensemble_analyzer._clustering.cluster_manager import ClusteringConfig, ClusteringManager
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(message)s")
    logger = logging.getLogger("enan_clustering")

    print(f"\n{'=' * 60}")
    print("PCA CLUSTERING ANALYSIS")
    print(f"{'=' * 60}\n")

    logger.info(f"Loading ensemble from {args.file}...")

    try:
        ensemble = read_ensemble(args.file, logger, raw=True)
        logger.info(f"✓ Loaded {len(ensemble)} conformers")
    except Exception as e:
        logger.error(f"Failed to load ensemble: {e}")
        exit(1)

    config = ClusteringConfig(
        include_H=args.include_H,
        set_cluster_attribute=True,
        n_jobs=max(1, args.cpu),
    )

    manager = ClusteringManager(logger=logger, config=config)

    result = manager.perform_pca(
        conformers=ensemble,
        n_clusters=args.ncluster,
        output_file=args.output,
        title=args.title,
        include_legend=args.legend,
    )

    if result:
        out_dir = os.path.dirname(args.output) or "."
        base_name = os.path.splitext(os.path.basename(args.output))[0]
        plot_clustering_metrics(result, out_dir, base_name, logger, explicit_k=args.ncluster)
        plot_component_analysis(result, out_dir, base_name, logger)
        plot_before_after_pca(result, out_dir, base_name, logger)

        print(f"\n{'=' * 60}")
        print("RESULTS")
        print(f"{'=' * 60}")
        print(f"✓ Clusters identified: {result.n_clusters}")
        print(f"✓ Variance explained (PC1): {result.explained_variance[0] * 100:.1f}%")
        print(f"✓ Variance explained (PC2): {result.explained_variance[1] * 100:.1f}%")
        print(f"✓ Total variance (PC1+PC2): {result.explained_variance[:2].sum() * 100:.1f}%")
        print(f"✓ Output saved: {args.output}")

        print(f"\n{'=' * 60}")
        print("COMPONENT ANALYSIS")
        print(f"{'=' * 60}")
        for i in range(min(3, len(result.components))):
            pc = result.components[i]
            top_features = np.argsort(np.abs(pc))[-5:][::-1]
            print(f"PC{i + 1} top features: {list(top_features)} "
                  f"(variance: {result.explained_variance[i] * 100:.1f}%)")
        print(f"{'=' * 60}")

        if args.write:
            logger.info("\nReducing ensemble by clusters...")
            reduced = manager.reduce_by_clusters(ensemble, True)
            save_snapshot("clustered.xyz", reduced, logger)
            active = len([c for c in reduced if c.active])
            logger.info(f"✓ Reduced ensemble: {active} conformers → clustered.xyz")

        print(f"{'=' * 60}\n")
    else:
        logger.error("PCA analysis failed or was skipped")
        exit(1)


if __name__ == '__main__':
    main()
