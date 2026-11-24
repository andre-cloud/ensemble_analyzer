"""Starting geometry of an NMR reference molecule.

Priority: the user's ``xyz`` file, then the geometry bundled with EnAn when ``name`` is a
standard reference (``parameters_file/nmr_references``, case-insensitive, with aliases), then
``smiles`` through RDKit (ETKDG + MMFF/UFF). The geometry is only a starting point: EnAn
optimises the reference before computing its shielding.
"""

import json
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np

_DIR = files("ensemble_analyzer").joinpath("parameters_file/nmr_references")
_INDEX: Dict[str, Dict] = json.loads(_DIR.joinpath("index.json").read_text())

RDKIT_HINT = (
    "SMILES references need RDKit: conda install -c conda-forge rdkit "
    "(or pip install \"ensemble-analyzer[nmr]\")"
)


@dataclass
class ReferenceGeometry:
    symbols: List[str]
    coords: np.ndarray
    charge: int
    mult: int
    source: str


def _key(text: str) -> str:
    return "".join(ch for ch in str(text).lower() if ch.isalnum())


def bundled_reference(name: str) -> Optional[Dict]:
    """Entry of a bundled reference matching ``name`` or one of its aliases, else None."""
    wanted = _key(name)
    for canonical, entry in _INDEX.items():
        if wanted in {_key(canonical), *(_key(a) for a in entry.get("aliases", []))}:
            return {"canonical": canonical, **entry}
    return None


def rdkit_available() -> bool:
    try:
        import rdkit  # noqa: F401
    except ImportError:
        return False
    return True


def needs_rdkit(ref: Dict) -> bool:
    """True when the geometry of ``ref`` can only come from its SMILES."""
    return not ref.get("xyz") and bundled_reference(ref.get("name", "")) is None and bool(ref.get("smiles"))


def _read_xyz(text: str) -> tuple:
    lines = text.splitlines()
    n = int(lines[0].split()[0])
    symbols, coords = [], []
    for line in lines[2:2 + n]:
        parts = line.split()
        symbols.append(parts[0].capitalize())
        coords.append([float(x) for x in parts[1:4]])
    return symbols, np.array(coords, dtype=float)


def _from_smiles(smiles: str, seed: int = 7) -> ReferenceGeometry:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
    except ImportError as exc:
        raise ImportError(RDKIT_HINT) from exc
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for NMR reference: {smiles!r}")
    mol = Chem.AddHs(mol)
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=seed))
    if not cids:
        raise ValueError(f"RDKit could not embed the NMR reference {smiles!r}")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
    else:
        res = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=2000)
    best = cids[int(np.argmin([e for _, e in res]))]
    conf = mol.GetConformer(best)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    return ReferenceGeometry(
        symbols, coords, Chem.GetFormalCharge(mol),
        Descriptors.NumRadicalElectrons(mol) + 1, f"smiles:{smiles}",
    )


def resolve_reference_geometry(ref: Dict, base_dir: Union[str, Path, None] = None) -> ReferenceGeometry:
    """Starting geometry, charge and multiplicity of a reference (see module docstring).

    ``charge``/``mult`` in ``ref`` override the bundled or SMILES-derived values.
    """
    if ref.get("xyz"):
        path = Path(ref["xyz"])
        if not path.is_absolute() and base_dir is not None:
            path = Path(base_dir) / path
        symbols, coords = _read_xyz(path.read_text())
        geom = ReferenceGeometry(symbols, coords, 0, 1, f"xyz:{ref['xyz']}")
    elif (entry := bundled_reference(ref.get("name", ""))) is not None:
        symbols, coords = _read_xyz(_DIR.joinpath(entry["file"]).read_text())
        geom = ReferenceGeometry(symbols, coords, int(entry.get("charge", 0)),
                                 int(entry.get("mult", 1)), f"bundled:{entry['canonical']}")
    elif ref.get("smiles"):
        geom = _from_smiles(ref["smiles"])
    else:
        raise ValueError(
            f"NMR reference {ref.get('name')!r}: give 'xyz' or 'smiles', or use a bundled "
            f"reference ({', '.join(_INDEX)})"
        )
    if ref.get("charge") is not None:
        geom.charge = int(ref["charge"])
    if ref.get("mult") is not None:
        geom.mult = int(ref["mult"])
    return geom
