import os, re

import numpy as np
from typing import List, Tuple, Optional

from ensemble_analyzer._conformer.conformer import Conformer
from ensemble_analyzer._logger.logger import Logger


def _parse_xyz_str(fl: List[str], raw: bool =False) -> Tuple[np.ndarray, np.ndarray, Optional[float]]:
    """Parse an xyz geom descriptor

    Args:
        fl (List[str]): Lines of the xyzfile of only one geometry
        raw (bool, optional): Convert the energy in the comment line. Defaults to False.

    Returns:
        Tuple[np.ndarray, np.ndarray, float]: Atoms, Geometry and eventually Energy
    """
    e = None
    if raw:
        e = float(re.findall(r'([- ]\d*\.\d*)$', fl[1].strip())[0])
    fl = fl[2:]
    atoms, geom = [], []
    for line in fl:
        a, *g = line.split()
        atoms.append(a)
        geom.append(g)
    return np.array(atoms), np.array(geom, dtype=float), e


def read_ensemble(file: str, log:Logger, raw: bool=False) -> list:
    """
    Read the initial ensemble and return the ensemble list
    Not only XYZ file is supported. OBABEL is required

    :param file: initial ensemble file
    :type file: str
    :param charge: charge of the molecule
    :type charge: int
    :param multiplicity: multiplicity of the molecule
    :type multiplicity: int
    :param log: logger instance
    :type log: logger

    :return: whole ensemble list as Conformer instances
    :rtype: list
    """

    confs = []

    if not file.endswith(".xyz"):
        raise "Ensemble file must be an XYZ (multi)geometry file"

    with open(file) as f:
        fl = f.readlines()

    n_atoms = int(fl[0])
    old_idx = 0
    counter = 1
    for i in range(0, len(fl) + 1, n_atoms + 2):
        if i == old_idx:
            continue
        atoms, geom, e = _parse_xyz_str(fl[old_idx:i], raw=raw)
        confs.append(Conformer(counter, geom=geom, atoms=atoms))
        old_idx = i
        counter += 1

    return confs


def save_snapshot(output: str, confs: List[Conformer], log: Logger):
    """
    Save an XYZ file to store a bunch of geometries

    :param output: output filename
    :type output: str
    :param confs: list of all active conformers
    :type confs: list(Conformer)
    :param log: logger instance
    :type log: logger
    :return: None
    """
    log.debug("Saving snapshot of the ensemble")
    xyzs = []
    for conf in confs:
        xyz_data = conf.write_xyz()
        if xyz_data:
            xyzs.append(xyz_data)
    with open(output, "w") as f:
        f.write("\n".join(xyzs))
    return None