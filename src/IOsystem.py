import numpy as np
import json
import os
import shutil


def _parse_xyz_str(fl: str, raw=False):
    """
    Parse an xyz geom descriptor

    :param fl: string of the file splitted in a single geometry
    :type fl: str

    :return: list of atoms, XYZ position of the atoms
    :rtype: tuple
    """
    e = None
    if raw:
        e = float(fl[1].split()[0].strip())
    fl = fl[2:]
    atoms, geom = [], []
    for line in fl:
        a, *g = line.split()
        atoms.append(a)
        geom.append(g)
    return np.array(atoms), np.array(geom, dtype=float), e


def mkdir(directory: str):
    """
    Create a directory, raising an error if the directory already exists

    :param directory: directory name
    :type directory: str

    :return: None
    """
    # if os.path.exists(directory):
    #     raise IOError(f"Directory {directory} already exists. Going to exit!")
    # shutil.rmtree(directory)
    os.makedirs(directory, exist_ok=True)
    return None


class SerialiseEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj.__dict__


def tail(file_path, num_lines=100):
    """Tail an output file

    :param file_path: output filename path
    :type file_path: str
    :param num_lines: number of lines to print
    :type num_lines: int
    :return: the tail of the file
    :rtype: str
    """
    with open(file_path) as f:
        fl = f.readlines()

    return "".join(fl[-num_lines:])


def move_files(conf, protocol, label):
    files = [str(f) for f in os.listdir(os.getcwd()) if str(f).startswith(label)]
    dest_folder = os.path.join(os.getcwd(), conf.folder, f"protocol_{protocol.number}")
    mkdir(os.path.join(os.getcwd(), conf.folder, f"protocol_{protocol.number}"))
    for file in files:
        src = os.path.join(os.getcwd(), file)
        dst = os.path.join(dest_folder, f"{conf.number}_p{protocol.number}_{file}")
        shutil.copyfile(src, dst)


if __name__ == "__main__":
    print(tail("orca.out", 4))
