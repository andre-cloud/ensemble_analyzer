import shutil
import os
from pathlib import Path
from typing import Any, Tuple

from ase.calculators.nwchem import NWChem
from ensemble_analyzer.calculators.base import BaseCalc, register_calculator
from ensemble_analyzer.protocol.dispersion import dispersion_keyword


NWCHEM_COMMAND = (
    os.getenv("NWCHEMCOMMAND")
    or shutil.which("nwchem_openmpi")
    or shutil.which("nwchem")
)


@register_calculator("nwchem")
class NWChemCalc(BaseCalc):

    def common_str(self) -> dict:
        nroots = self.protocol.nroots
        is_tddft = isinstance(nroots, int) and nroots > 0
        kw = {
            "theory": "tddft" if is_tddft else "dft",
            "xc": self.protocol.functional,
            "basis": self.protocol.basis,
            "dft": {"mult": self.protocol.mult},
        }
        disp = dispersion_keyword(self.protocol.dispersion, "nwchem")
        if disp:
            kw["dft"]["disp"] = {"vdw": disp}
        if is_tddft:
            kw["nroots"] = nroots
            tda = self.protocol.tda
            kw["tddft"] = {"tda": True} if isinstance(tda, bool) and tda else True
        if self.protocol.charge != 0:
            kw["charge"] = self.protocol.charge

        if self.protocol.solvent:
            solv = self.protocol.solvent.solvent
            if solv:
                kw["cosmo"] = {"solvent": solv.lower()}

        memory_mb = self.cpu * 5000
        kw["memory"] = f"{memory_mb} mb"

        kw["dft"]["vectors"] = self._vectors_directive()

        return kw

    def _movecs_path(self, protocol_number: Any) -> str:
        """MO vectors written by the NWChem job of a protocol step.

        ASE uses ``<conf>/protocol_<n>/<label>`` as ``permanent_dir`` and NWChem
        names the vectors after the ``start`` label: ``<label>.movecs``.
        """
        label = f"{self.conf.number}_p{protocol_number}_nwchem"
        return self._build_path(
            self.conf.folder, f"protocol_{protocol_number}", label, f"{label}.movecs"
        )

    def _vectors_directive(self) -> str:
        """``vectors`` line of the dft block (initial guess from ``read_orbitals``).

        The output file is always given explicitly: without it NWChem overwrites
        the input vectors of the previous step.
        """
        current = self._movecs_path(self.protocol.number)
        previous = getattr(self.protocol, "read_orbitals", None)
        # step "0" / 0 is a valid reference: test for emptiness, not truthiness
        if previous is not None and previous is not False and str(previous).strip():
            return f'input "{self._movecs_path(previous)}" output "{current}"'
        return f'output "{current}"'

    def _build_constraints(self) -> str:
        if not self.constraints:
            return ""

        if isinstance(self.constraints, str):
            raw = self.constraints
        elif isinstance(self.constraints, list):
            if all(isinstance(c, list) for c in self.constraints):
                atoms = []
                for constraint in self.constraints:
                    atoms.extend(str(idx + 1) for idx in constraint)
                raw = f"fix atom {' '.join(atoms)}"
            elif all(isinstance(c, str) for c in self.constraints):
                raw = "\n".join(self.constraints)
            else:
                raw = str(self.constraints)
        else:
            raw = str(self.constraints)

        if "constraints" not in raw.lower():
            return f"constraints\n  {raw}\nend"
        return raw

    def _build_calculator(self) -> Tuple[Any, str]:
        from enan_calculators import get_ase_calculator

        ase_label = f"{self.conf.folder}/protocol_{self.protocol.number}/{self.conf.number}_p{self.protocol.number}_nwchem"
        command = NWCHEM_COMMAND
        if "nwchem_openmpi" in command and not any(
            x in command for x in ("mpirun", "mpiexec")
        ):
            command = f"mpirun -np {self.cpu} {command}"
        command = f"{command} {self.conf.number}_p{self.protocol.number}_nwchem.nwi > {self.conf.number}_p{self.protocol.number}_nwchem.nwo"

        calculator = get_ase_calculator(
            "nwchem",
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            method=self.protocol.functional,
            basis=self.protocol.basis,
            solvent=self.protocol.solvent.solvent if self.protocol.solvent else None,
            cpu=self.cpu,
            add_input=self.protocol.add_input,
            label=ase_label,
            command=command,
        )

        # `disp vdw N` must sit inside the DFT block, not in the appended add_input
        disp = dispersion_keyword(self.protocol.dispersion, "nwchem")
        if disp:
            calculator.parameters.setdefault("dft", {})["disp"] = {"vdw": disp}

        # Initial guess from a previous step (read_orbitals) + explicit output vectors
        calculator.parameters.setdefault("dft", {})["vectors"] = self._vectors_directive()

        constraints_block = self._build_constraints()
        if constraints_block:
            original = calculator.write_input

            def patched_write_input(atoms, properties=None, system_changes=None) -> None:
                original(atoms, properties, system_changes)
                inp = Path(calculator.directory) / calculator.input_filename()
                with open(inp, "a") as f:
                    f.write("\n" + constraints_block + "\n")

            calculator.write_input = patched_write_input

        return calculator, "nwchem"

    def _add_sp_keywords(self, calc: NWChem) -> None:
        if "task" not in self.protocol.add_input:
            calc.parameters["task"] = "energy"

    def _add_opt_keywords(self, calc: NWChem) -> None:
        calc.parameters["task"] = "optimize ts" if self.protocol.ts else "optimize"
        if self.protocol.ts:
            add_input = self.protocol.add_input or ""
            if 'inhess' not in add_input:
                src = self._find_hessian_source_protocol()
                calc.parameters["driver"] = {"inhess": 2 if src is not None else 1}
        if self.protocol.freq:
            nroots = self.protocol.nroots
            freq_task = "tddft" if (isinstance(nroots, int) and nroots > 0) else "dft"
            calc.parameters["task"] += f"\ntask {freq_task} freq"

    def _add_freq_keywords(self, calc: NWChem) -> None:
        calc.parameters["task"] = "freq"

    def _add_tddft_keywords(self, calc: NWChem) -> None:
        """TD-DFT single points (``task tddft energy``) when ``nroots`` is set.

        Optimisation/frequency steps keep their own task strings (see
        ``_add_opt_keywords``).
        """
        nroots = self.protocol.nroots
        if not (isinstance(nroots, int) and nroots > 0):
            return
        if self.protocol.opt or calc.parameters.get("task") != "energy":
            return
        tda = self.protocol.tda
        calc.parameters["theory"] = "tddft"
        calc.parameters["tddft"] = {
            "nroots": nroots,
            ("cis" if isinstance(tda, bool) and tda else "rpa"): None,
        }
