import os
import shutil
from ase.calculators.nwchem import NWChem

from enan_calculators._input import append_to_input

NWCHEM_COMMAND = (
    os.getenv("NWCHEMCOMMAND")
    or shutil.which("nwchem_openmpi")
    or shutil.which("nwchem")
)


def create_nwchem_calc(charge, mult, method, basis, solvent=None, cpu=1, add_input="", label=None, command=None):
    kw = {
        "theory": "dft",
        "xc": method,
        "basis": basis,
        "dft": {"mult": mult},
    }
    if charge != 0:
        kw["charge"] = charge
    if solvent:
        kw["cosmo"] = {"solvent": solvent.lower()}
    memory_mb = cpu * 5000
    kw["memory"] = f"{memory_mb} mb"

    calculator = NWChem(label=label or "nwchem", command=command or NWCHEM_COMMAND, **kw)

    if add_input.strip():
        append_nwchem_block(calculator, add_input)

    return calculator


def append_nwchem_block(calculator, text: str) -> None:
    """Append ``text`` to the NWChem input (chains earlier patches)."""
    append_to_input(calculator, text, "write_input", calculator.input_filename)
