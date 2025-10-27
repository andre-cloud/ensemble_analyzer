#!/opt/miniconda3/bin/python

import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import rcParams
from typing import Union
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata
from scipy.spatial import distance_matrix

import sys
import pickle as pl


plt.set_loglevel("error")

MARKERS = [
    ".",
    ",",
    "o",
    "v",
    "^",
    "<",
    ">",
    "1",
    "2",
    "3",
    "4",
    "8",
    "s",
    "p",
    "*",
    "h",
    "H",
    "+",
    "x",
    "D",
    "d",
    "|",
    "_",
    "P",
    "X",
]


def calc_distance_matrix(coords, atoms, include_H=True):
    dist = np.zeros((coords.shape[0], coords.shape[1], coords.shape[1]))
    evalue_dist, evector_dist = [], []

    for idx, _ in enumerate(coords):
        c = coords[idx] if include_H else coords[idx][atoms!='H']
        dist[idx] = distance_matrix(c, c)
        eva, eve = np.linalg.eig(dist[idx])
        evalue_dist.append(eva)
        evector_dist.append(eve)

    # diff = coords[:, :, np.newaxis, :] - coords[:, np.newaxis, :, :]
    # dist = np.linalg.norm(diff, axis=-1)
    # print(dist)

    return np.array(evalue_dist), np.array(evector_dist)


def get_best_ncluster(coords):
    """
    Obtain the best number of cluster based on the maximization of the silhouette

    :param coords: array of the coordinates of each atom
    :type coords: 2D-array

    :return: best number of clusters
    :rtype: int
    """

    k_range = range(10, 30)
    silhouette_scores = []

    for k in k_range:
        kmeans = KMeans(n_clusters=k, n_init=10)
        labels = kmeans.fit_predict(coords)
        score = silhouette_score(coords, labels)
        silhouette_scores.append(score)
    return k_range[np.argmax(silhouette_scores)]


def calc_pca(confs: list, cluster=False, ncluster: Union[int, None] = None, set = True, include_H = True) -> tuple:
    """
    Function that execute the actual PCA analysis.
    It wants to understand how conformations differ from each other based on their overall Cartesian coordinates

    :param confs: whole list of the confomers
    :type confs: list
    :param cluster: execute the cluster action. Default is False
    :type cluster: bool, optional
    :param ncluster: number of cluster to form using the KMean analysis. Defualt is None
    :type ncluster: int

    :return: PCA transformation (pca_energy), Clustered coordinates (clusters), colors and number of the conformer and energy
    :rtype: tuple
    """

    # fetch all geometries and reshaping them to create the correct 2D matrix
    data = np.array([conf.last_geometry for conf in confs if conf.active])
    colors = [conf.color for conf in confs if conf.active]
    numbers = [conf.number for conf in confs if conf.active]
    energy = np.array([conf.get_energy for conf in confs if conf.active])
    energy -= min(energy)

    evalue_dist, _ = calc_distance_matrix(data, atoms=confs[0].atoms, include_H=include_H)

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

        # Cluster the data
        if not confs[0].cluster:
            kmeans = KMeans(n_clusters=n_c, n_init=10)
            clusters = kmeans.fit_predict(pca_scores)
            if set: 
                for idx, conf in enumerate(confs):
                    conf.cluster = int(clusters[idx])
        else:
            clusters = [conf.cluster for conf in confs]

    else:
        clusters = [1 for _ in confs]

    return pca_scores, clusters, colors, numbers, energy


def obtain_markers_from_cluster(cluster: int):
    """
    Obtain a different marker from the marker library for different conformers

    :param cluster: the cluster number
    :return: marker
    :rtype: matplotlib.lines
    """
    return MARKERS[cluster % (len(MARKERS))]


def save_PCA_snapshot(
    fname: str,
    title: str,
    pca_scores: np.ndarray,
    clusters: np.ndarray,
    colors: list,
    numbers: list,
    z: list,
    legend: bool = True
):
    """
    Graph and save the image of the PCA analysis

    :param fname: filename to save the graphs
    :type fname: str
    :param title: title of the graph
    :type title: str
    :param pca_scores: PCA transformation
    :type pca_scores: np.array
    :param clusters: Clustered coordinates
    :type clusters: np.array
    :param z: list of energy of the conformers
    :type z: np.array

    :rtype: None
    """

    fig = plt.figure(figsize=(10,8))
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


    with open(fname.replace('.png', '.pickle'), 'wb') as f:
            pl.dump(fig, f)

    plt.savefig(fname, dpi=300)
    return None


def perform_PCA(confs: list, ncluster: int, fname: str, title: str, log, set=True, include_H=True, legend=True) -> None:
    """
    Perform a PCA analysis

    :param confs:  list of all the active conformers
    :type confs: list
    :param ncluster:  number of cluster to group the ensemble
    :type ncluster: int
    :param fname:  filename for the graph
    :type fname: str
    :param title:  title of the graph
    :type title: str
    :param log:  logger instance
    :type log: logging

    :rtype: None
    """
    log.info("Starting PCA analysis")
    nc = ncluster if len(confs) > ncluster else len(confs) - 1
    if nc <= 2:
        return None
    pca_scores, clusters, colors, numbers, energy = calc_pca(
        confs, ncluster=nc, cluster=True, set=set, include_H=include_H
    )
    save_PCA_snapshot(fname, title, pca_scores, clusters, colors, numbers, energy, legend=legend)

    return None


def get_ensemble(confs, sort=False):
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
    from ioFile import read_ensemble
    from ioFile import save_snapshot
    import sys
    import mock
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='Ensemble to be clusterd')
    parser.add_argument('-nc', '--ncluster', help='Number of families to cluster. Defaults 5', default=5, type=int)
    parser.add_argument('--no-H', help='Exclude hydrogen atoms in the PCA', action='store_false')
    parser.add_argument('--no-legend', help='Exclude legend from PCA graph', action='store_false')
    parser.add('--title', help='Title for the PCA graph', Default='Cluster')
    args = parser.parse_args()

    # Load the XYZ file
    xyz_file = read_ensemble(args.file, mock.MagicMock(), raw=True)

    perform_PCA(xyz_file, args.ncluster, "cluster.png", args.title, mock.MagicMock(), include_H=args.no_H, legend=args.no_legend)

    xyz_file_new = get_ensemble(xyz_file)

    save_snapshot("clustered.xyz", xyz_file_new, mock.MagicMock())
