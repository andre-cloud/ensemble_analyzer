import os
import shutil
from pathlib import Path
from ase.calculators.nwchem import NWChem

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

    extra = []
    if add_input.strip():
        extra.append(add_input)

    if extra:
        block = "\n\n".join(extra)
        original = calculator.write_input

        def patched_write_input(atoms, properties=None, system_changes=None) -> None:
            original(atoms, properties, system_changes)
            inp = Path(calculator.directory) / calculator.input_filename()
            with open(inp, "a") as f:
                f.write("\n" + block + "\n")

        calculator.write_input = patched_write_input

    return calculator
