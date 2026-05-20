import shutil
import os
from pathlib import Path
import warnings
from typing import Any, Tuple

from ase.calculators.nwchem import NWChem
from ensemble_analyzer._calculators.base import BaseCalc, register_calculator


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

    def common_str(self) -> dict:
        """Build the common NWChem input keyword dictionary.

        Includes theory, functional, basis, memory, and optional solvent,
        charge, and orbital reading directives.

        Returns:
            dict: NWChem input keywords.
        """
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

        memory_mb = self.cpu * 5000
        kw["memory"] = f"{memory_mb} mb"

        current_movecs = os.path.abspath(
            f"{self.conf.folder}/protocol_{self.protocol.number}/{self.conf.number}_p{self.protocol.number}_nwchem/{self.conf.number}_p{self.protocol.number}_nwchem.movecs"
        ).replace('\\', '/')

        if getattr(self.protocol, "read_orbitals", None):
            # If read_orbitals is an int, read from that specific protocol. Otherwise, default to previous.
            prev_p_num = self.protocol.read_orbitals
            old_movecs = os.path.abspath(
                f"{self.conf.folder}/protocol_{prev_p_num}/{self.conf.number}_p{self.protocol.number}_nwchem/{self.conf.number}_p{prev_p_num}_nwchem.movecs"
            ).replace('\\', '/')
            
            kw["dft"]["vectors"] = f'input "{old_movecs}" output "{current_movecs}"'
        else:
            kw["dft"]["vectors"] = f'output "{current_movecs}"'

        return kw

    def _build_constraints(self) -> str:
        """Build the NWChem constraints block from the protocol constraints.

        Supports freeze-cartesian (list of lists of atoms), bond/angle/dihedral
        constraints, and raw constraint strings.

        Returns:
            str: NWChem constraints block, or empty string if no constraints.
        """
        if not self.constrains:
            return ""

        if isinstance(self.constrains, str):
            raw = self.constrains
        elif isinstance(self.constrains, list):
            if all(isinstance(c, list) for c in self.constrains):
                atoms = []
                for constraint in self.constrains:
                    atoms.extend(str(idx + 1) for idx in constraint)
                raw = f"fix atom {' '.join(atoms)}"
            elif all(isinstance(c, str) for c in self.constrains):
                raw = "\n".join(self.constrains)
            else:
                raw = str(self.constrains)
        else:
            raw = str(self.constrains)

        if "constraints" not in raw.lower():
            return f"constraints\n  {raw}\nend"
        return raw

    def _std_calc(self) -> Tuple[NWChem, str]:
        """Build a standard NWChem calculator instance with common keywords.

        Handles MPI command construction, patched write_input for extra
        input blocks and constraints.

        Returns:
            Tuple[NWChem, str]: ASE NWChem calculator and label string.
        """
        kw = self.common_str()
        ase_label = f"{self.conf.folder}/protocol_{self.protocol.number}/{self.conf.number}_p{self.protocol.number}_nwchem"

        command = NWCHEM_COMMAND
        if "nwchem_openmpi" in command and not any(
            x in command for x in ("mpirun", "mpiexec")
        ):
            command = f"mpirun -np {self.cpu} {command}"
        command = f"{command} {self.conf.number}_p{self.protocol.number}_nwchem.nwi > {self.conf.number}_p{self.protocol.number}_nwchem.nwo"

        calculator = NWChem(label=ase_label, command=command, **kw)

        extra = []
        if self.protocol.add_input.strip():
            extra.append(self.protocol.add_input)

        constraints_block = self._build_constraints()
        if constraints_block:
            extra.append(constraints_block)

        if extra:
            block = "\n\n".join(extra)
            original = calculator.write_input

            def patched_write_input(atoms: Any, properties: Any = None, system_changes: Any = None) -> None:
                """Patch write_input to append extra blocks to the NWChem input file."""
                original(atoms, properties, system_changes)
                inp = Path(calculator.directory) / calculator.input_filename()
                with open(inp, "a") as f:
                    f.write("\n" + block + "\n")

            calculator.write_input = patched_write_input

        return calculator, ase_label

    def single_point(self) -> Tuple[NWChem, str]:
        """Configure a single-point energy calculation with NWChem.

        Returns:
            Tuple[NWChem, str]: ASE NWChem calculator and label string.
        """
        calc, label = self._std_calc()
        if "task" not in self.protocol.add_input:
            calc.parameters["task"] = "energy"
        
        return calc, label

    def optimisation(self) -> Tuple[NWChem, str]:
        """Configure a geometry optimisation with NWChem.

        Returns:
            Tuple[NWChem, str]: ASE NWChem calculator and label string.
        """
        calc, label = self._std_calc()
        calc.parameters["task"] = "optimize"
        if self.protocol.freq: 
            calc.parameters["task"] += "\ntask dft freq"
        return calc, label

    def frequency(self) -> Tuple[NWChem, str]:
        """Configure a frequency calculation with NWChem.

        Returns:
            Tuple[NWChem, str]: ASE NWChem calculator and label string.
        """
        calc, label = self._std_calc()
        calc.parameters["task"] = "freq"
        return calc, label