import shutil
import os

from ase.calculators.nwchem import NWChem
from ensemble_analyzer._calculators.base import BaseCalc, register_calculator

from typing import Tuple


NWCHEM_COMMAND = (
    os.getenv("NWCHEMCOMMAND")
    or shutil.which("nwchem_openmpi")
    or shutil.which("nwchem")
)


@register_calculator("nwchem")
class NWChemCalc(BaseCalc):
    """
    Calculator wrapper for NWChem.
    Handles input generation for SP, OPT, and FREQ jobs.
    """

    label = "nwchem"

    def common_str(self) -> dict:
        kw = {
            "theory": "dft",
            "xc": self.protocol.functional,
            "basis": self.protocol.basis,
            "dft": {"mult": self.protocol.mult},
        }
        if self.protocol.charge != 0:
            kw["charge"] = self.protocol.charge

        if self.protocol.solvent:
            solv = self.protocol.solvent.solvent
            if solv:
                kw["cosmo"] = {"solvent": solv.lower()}

        memory_mb = self.cpu * 2000
        kw["memory"] = f"{memory_mb} mb"

        return kw

    def _std_calc(self) -> Tuple[NWChem, str]:
        kw = self.common_str()
        ase_label = f"{self.conf.folder}/protocol_{self.protocol.number}/{self.conf.number}_p{self.protocol.number}_nwchem"

        command = NWCHEM_COMMAND
        if "nwchem_openmpi" in command and not any(
            x in command for x in ("mpirun", "mpiexec")
        ):
            command = f"mpirun -np {self.cpu} {command}"
        command = f"{command} PREFIX.nwi > PREFIX.nwo"

        calculator = NWChem(label=ase_label, command=command, **kw)

        if self.protocol.add_input.strip():
            add_input = self.protocol.add_input
            original = calculator.write_input

            def patched_write_input(atoms, properties=None, system_changes=None):
                original(atoms, properties, system_changes)
                from pathlib import Path
                inp = Path(calculator.directory) / calculator.input_filename()
                with open(inp, "a") as f:
                    f.write("\n" + add_input + "\n")

            calculator.write_input = patched_write_input

        return calculator, "nwchem"

    def single_point(self) -> Tuple[NWChem, str]:
        return self._std_calc()

    def optimisation(self) -> Tuple[NWChem, str]:
        calc, label = self._std_calc()
        calc.parameters["task"] = "optimize"
        if self.protocol.freq:
            raise NotImplementedError("Frequency in NWChem are NOT supported for now in EnAn.")
        return calc, label

    def frequency(self) -> Tuple[NWChem, str]:
        raise NotImplementedError("Frequency in NWChem are NOT supported for now in EnAn.")
