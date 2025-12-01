from pathlib import Path
import shutil, json

import numpy as np


def mkdir(directory: str) -> bool:
    """Create a directory, raising an error if the directory already exists

    Args:
        directory (str): Folder to be created
    """

    # os.makedirs(directory, exist_ok=True)
    Path(directory).mkdir(parents=True, exist_ok=True)

    return True


def move_files(conf, protocol, label: str) -> None:
    """Move files from one directory to another

    Args:
        conf (Conformer): Conformer instance
        protocol (Protocol): Protocol instance
        label (str): label of the calculator
    """

    cwd = Path.cwd()
    files = [f for f in cwd.iterdir() if f.name.startswith(label)]
    dest_folder = cwd / conf.folder / f"protocol_{protocol.number}"
    mkdir(str(dest_folder))

    for file in files:
        dst = dest_folder / f"{conf.number}_p{protocol.number}_{file.name}"
        shutil.move(str(file), str(dst))


def tail(file_path: str, num_lines:int=100) -> str:
    """Tail an output file

    Args:
        file_path (str): File to be tailed
        num_lines (int, optional): Number of lines to be shown. Defaults to 100.

    Returns:
        str: Tail
    """
    with Path(file_path).open() as f:
        fl = f.readlines()

    return "".join(fl[-num_lines:])


class SerialiseEncoder(json.JSONEncoder):
    """Serialization Encoder
    """
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj.__dict__