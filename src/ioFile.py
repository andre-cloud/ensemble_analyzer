try:
    from src.IOsystem import _parse_xyz_str
    from src.conformer import Conformer
except ImportError as e:  # pragma: no cover
    print(e)
    from IOsystem import _parse_xyz_str
    from conformer import Conformer


import os


def convert_file(file) -> str:
    """
    Convert the input file into xyz multigeometry XYZ file.
    OPENBABEL is required

    :param file: input filename
    :type file: str

    :return: input converted filename
    :rtype: str
    """
    output = "_".join(file.split(".")[:-1]) + ".xyz"
    os.system(f"obabel {file} -O{output}")
    return output


def read_ensemble(file, log, raw=False) -> list:
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
        file = convert_file(file)

    with open(file) as f:
        fl = f.readlines()

    n_atoms = int(fl[0])
    old_idx = 0
    counter = 1
    for i in range(0, len(fl) + 1, n_atoms + 2):
        if i == old_idx:
            continue
        atoms, geom, e = _parse_xyz_str(fl[old_idx:i], raw=raw)
        confs.append(
            Conformer(counter, geom=geom, atoms=atoms)
        )
        if raw:
            confs[-1].energies = {"0": {"E": e * 627.51, "G": e * 627.51}}
        old_idx = i
        counter += 1

    return confs


def save_snapshot(output, confs, log):
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
