#!/opt/miniconda3/bin/python


from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.gridspec as gridspec

from typing import Union, List, Tuple


from scipy.interpolate import griddata
from scipy.spatial import distance_matrix

from src._conformer.conformer import Conformer
from src._logger.logger import Logger

import sys
import pickle as pl
import numpy as np

MARKERS = [
    ".", ",", "o", "v", "^", "<",
    ">", "1", "2", "3",
    "4", "8", "s", "p",
    "*", "h", "H", "+",
    "x", "D", "d", "|",
    "_", "P", "X",
]


def calc_distance_matrix(coords:np.ndarray, atoms:List[str], include_H:bool=True) -> Tuple[np.ndarray, np.ndarray]:
    """Calculate the eigenvalues and eigenvector of the distance matrix for each conformer

    Args:
        coords (np.ndarray): Array of geometries
        atoms (List[str]): Atom list
        include_H (bool, optional): Include H in the distance matrix. Defaults to True.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Array of the array of eigenvalues, Array of the array of the eigenvector
    """
    natoms = coords[0].shape if include_H else coords[0][atoms != "H"].shape
    dist = np.zeros((coords.shape[0], natoms[0], natoms[0]))
    evalue_dist, evector_dist = [], []

    for idx, _ in enumerate(coords):
        c = coords[idx] if include_H else coords[idx][atoms != "H"]
        dist[idx] = distance_matrix(c, c)
        eva, eve = np.linalg.eig(dist[idx])
        evalue_dist.append(eva)
        evector_dist.append(eve)

    return np.array(evalue_dist), np.array(evector_dist)


def get_best_ncluster(coords:List[List[int]]) -> int:
    """
    Obtain the best number of cluster based on the maximization of the silhouette

    Args:
        coords (List[List[int]]): List of list of the eigenvalues to perfomre the PCA

    Returns:
        int: Best number of cluster, defined with KMeans
    """

    k_range = range(10, 30)
    silhouette_scores = []

    for k in k_range:
        kmeans = KMeans(n_clusters=k, n_init='auto')
        labels = kmeans.fit_predict(coords)
        score = silhouette_score(coords, labels)
        silhouette_scores.append(score)
    return k_range[np.argmax(silhouette_scores)]


def calc_pca(confs: List[Conformer], cluster: bool = False, ncluster: Union[int, None] = None, set_: bool = True, include_H: bool =True,) -> Tuple[np.ndarray, Union[np.ndarray, List[int]], List[str], List[int], np.ndarray]:
    """
    Function that execute the actual PCA analysis.
    It wants to understand how conformations differ from each other based on their overall Cartesian coordinates

    Args:
        confs (List[Conformer]): Ensemble
        cluster (bool, optional): Execute the cluster action. Defaults to False.
        ncluster (Union[int, None], optional): Number of cluster to form using the KMean analysis. Defaults to None.
        set_ (bool, optional): Set the cluster number to the Conformer class. Defaults to True.
        include_H (bool, optional): Include H in the distance matrix. Defaults to True.

    Returns:
        Tuple[np.ndarray, Union[np.ndarray, List[int]], List[str], List[int], np.ndarray]: PCA scores, Cluster numbers, Colors of the conformer, Number of the conformer, Relative energy of the conformers 
    """
    

    

    # fetch all geometries and reshaping them to create the correct 2D matrix
    data = np.array([conf.last_geometry for conf in confs if conf.active])
    colors = [conf.color for conf in confs if conf.active]
    numbers = [conf.number for conf in confs if conf.active]
    energy = np.array([conf.energies.get_energy() for conf in confs if conf.active])
    energy -= min(energy)

    evalue_dist, _ = calc_distance_matrix(
        data, atoms=confs[0].atoms, include_H=include_H
    )

    # perform PCA analysis with number of components as minimum between number of
    # n. confs and whole geom
    pca = PCA(n_components=min(evalue_dist.shape[0], evalue_dist.shape[1]))
    pca.fit(evalue_dist)
    pca_scores = pca.transform(evalue_dist)

    # getting the best number of clusters or create an array of 1
    if cluster:
        if not ncluster:
            n_c = get_best_ncluster(evalue_dist)
        else:
            n_c = ncluster

        kmeans = KMeans(n_clusters=n_c, n_init="auto")
        clusters = kmeans.fit_predict(pca_scores)

        if set_:
            for idx, conf in enumerate(confs):
                conf.cluster = int(clusters[idx])
        else:
            if not any([conf.cluster is None  for conf in confs]): 
                clusters = [conf.cluster for conf in confs]
            else:
                clusters = [1 for _ in confs]

    else:
        clusters = [1 for _ in confs]

    return pca_scores, clusters, colors, numbers, energy


def obtain_markers_from_cluster(cluster: int):
    """Obtain a different marker from the marker library for different conformers

    Args:
        cluster (int): Cluster Number
    """
    return MARKERS[cluster % (len(MARKERS))]


def save_PCA_snapshot(fname: str, title: str, pca_scores: np.ndarray, clusters: np.ndarray, colors: List[str], numbers: List[int], z: np.ndarray, legend: bool = True,) -> None:  
    """
    Graph and save the image of the PCA analysis

    Args:
        fname (str): Filename to save the graphs
        title (str): Title of the graph
        pca_scores (np.ndarray): PCA transformation
        clusters (np.ndarray): Clustered coordinates
        colors (List[str]): List of colors for each conformer
        numbers (List[int]): List of number for each conformer
        z (np.ndarray): List of energy of the conformers
        legend (bool, optional): Plot the legend of the graph. Defaults to True.
    """
    

    fig = plt.figure(figsize=(10, 8))
    if legend:
        plt.subplots_adjust(bottom=0.3, right=0.6, left=0.115)
    rcParams.update({"figure.autolayout": True})

    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 0.1], hspace=0.3)

    ax = fig.add_subplot(gs[0])
    color_axis = fig.add_subplot(gs[1])

    x_ = pca_scores[:, 0]
    y_ = pca_scores[:, 1]

    resolution = 500  # Grid resolution
    xi = np.linspace(min(x_), max(x_), resolution)
    yi = np.linspace(min(y_), max(y_), resolution)
    xi, yi = np.meshgrid(xi, yi)
    zi = griddata((x_, y_), z, (xi, yi), method="linear")

    im = ax.pcolormesh(xi, yi, zi, shading="auto", cmap="coolwarm", alpha=0.75)
    ax.contour(xi, yi, zi, "--", levels=6, colors="grey", linewidths=0.5, alpha=0.6)

    cbar = plt.colorbar(im, cax=color_axis, orientation="horizontal")
    cbar.set_label("Potential energy [kcal/mol]")

    print(obtain_markers_from_cluster, clusters)
    for x, y, m, c, n in zip(
        pca_scores[:, 0],
        pca_scores[:, 1],
        np.array(list(map(obtain_markers_from_cluster, clusters))),
        colors,
        numbers,
    ):
        ax.scatter(x, y, c=c, marker=m, label=f"CONF {n}")

    ax.set_xlabel("Principal Component 1")
    ax.set_ylabel("Principal Component 2")
    ax.set_title(title)

    y = ax.get_ylim()
    x = ax.get_xlim()

    ax.vlines(0, y[0], y[1], "#353535", "--", alpha=0.2)
    ax.hlines(0, x[0], x[1], "#353535", "--", alpha=0.2)

    ax.set_xlim(x)
    ax.set_ylim(y)
    if legend:
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.05, 1.0),
            fancybox=True,
            shadow=True,
            ncol=min(max(1, len(numbers) // 10), 6),
            borderaxespad=0.0,
            fontsize=6,
        )
    else:
        plt.tight_layout()

    with open(fname.replace(".png", ".pickle"), "wb") as f:
        pl.dump(fig, f)

    plt.savefig(fname, dpi=300)
    return None


def perform_PCA(confs: List[Conformer], ncluster: int, fname: str, title: str, log: Logger, set_ : bool =True, include_H: bool=True, legend: bool =True) -> None:
    """
    Perform a PCA analysis

    Args:
        confs (List[Conformer]): List of all the active conformers
        ncluster (int): Number of cluster to group the ensemble
        fname (str): Filename for the graph
        title (str): Title of the graph
        log (Logger): Logger
        set_ (bool, optional): Set the cluster number to the conformer. Defaults to True.
        include_H (bool, optional): Include the H in the distance matrix. Defaults to True.
        legend (bool, optional): Plot the legend of the graph. Defaults to True.
    """
    log.info("Starting PCA analysis")

    PERFORM = len(confs) > ncluster if nc else True

    if not (len(confs)>10 and PERFORM) : 
        log.warning(f"{log.WARNING} PCA not performed, length of ensemble too small ({len(confs)}). To perform PCA, ensemble must be bigger than 10 conformers and Cluster must be smaller than ensemble length.")
        return None
    
    nc = None
    if ncluster:
        nc = ncluster if len(confs) > ncluster else None

    if nc:
        log.info(f"\tUsing number_of_cluster={nc}.")
    else:
        log.info(f"\tEstimating the best number for clustering.")

    log.info(
        f"\tThe cluster will be set: {set_}.\t\nPCA will include hydrogen atoms: {include_H}"
    )

    pca_scores, clusters, colors, numbers, energy = calc_pca(
        confs, ncluster=nc, cluster=True, set_=set_, include_H=include_H
    )

    save_PCA_snapshot(
        fname, title, pca_scores, clusters, colors, numbers, energy, legend=legend
    )

    return None


def get_ensemble(confs: List[Conformer], sort:bool=False) -> List[Conformer]:
    """Get the new ensemble with the clustered part

    Args:
        confs (List[Conformer]): Ensemble
        sort (bool, optional): Sort the enesemble by energy. Defaults to False.

    Returns:
        List[Conformer]: Ensemble it self
    """
    if sort:
        tmp = sorted(confs)
    else:
        tmp = confs[:]

    clust_pres = []
    for i in tmp:
        if not i.active:
            continue
        if i.cluster not in clust_pres:
            clust_pres.append(i.cluster)
        else:
            i.active = False

    return tmp


if __name__ == "__main__":  # pragma: no cover:
    from src.ensemble_io import read_ensemble
    from src.ensemble_io import save_snapshot
    import sys
    import mock
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Ensemble to be clusterd")
    parser.add_argument(
        "-nc",
        "--ncluster",
        help="Number of families to cluster. Defaults 5",
        default=None,
        type=int,
    )
    parser.add_argument(
        "--no-H", help="Exclude hydrogen atoms in the PCA", action="store_false"
    )
    parser.add_argument(
        "--no-legend", help="Exclude legend from PCA graph", action="store_false"
    )
    parser.add_argument("--title", help="Title for the PCA graph", default="Cluster")
    args = parser.parse_args()

    # Load the XYZ file
    xyz_file = read_ensemble(args.file, mock.MagicMock(), raw=True)

    perform_PCA(
        xyz_file,
        args.ncluster,
        "cluster.png",
        args.title,
        mock.MagicMock(),
        include_H=args.no_H,
        legend=args.no_legend,
    )

    xyz_file_new = get_ensemble(xyz_file)

    save_snapshot("clustered.xyz", xyz_file_new, mock.MagicMock())
