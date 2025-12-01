from ase.calculators.orca import ORCA, OrcaProfile
from ensemble_analyzer._calculators.base import BaseCalc, register_calculator
import shutil
import os

VERSION = None

try:
    ORCA_COMMAND = os.getenv("ORCACOMMAND") or shutil.which("orca")
    orca_profile = OrcaProfile(command=ORCA_COMMAND)
    VERSION = int(os.getenv("ORCAVERSION")[0])
except Exception:
    orca_profile = None


@register_calculator("orca")
class OrcaCalc(BaseCalc):
    label = "orca"
    VERSION = VERSION if VERSION else 0

    def common_str(self):

        if self.protocol.solvent:
            if "xtb" in self.protocol.functional.lower():
                solv = f"ALPB({self.protocol.solvent.solvent})"
            elif self.protocol.solvent.solvent.strip():
                solv = f" {self.protocol.solvent}"
            else:
                solv = f" CPCM"
        else:
            solv = ""

        si = f"{self.protocol.functional} {self.protocol.basis} {solv} nopop"

        ob = (
            f"%pal nprocs {self.cpu} end "
            + self.protocol.add_input.format(CONF=self.conf.folder)
            + (" %maxcore 5000" if "maxcore" not in self.protocol.add_input else "")
        )

        return si, ob

    def _std_calc(self):
        si, ob = self.common_str()

        calculator = ORCA(
            profile=orca_profile,
            label="orca",
            orcasimpleinput=si,
            orcablocks=ob,
            charge=self.protocol.charge,
            mult=self.protocol.mult,
        )

        if self.protocol.read_orbitals:
            calculator.parameters["orcasimpleinput"] += " moread"
            calculator.parameters[
                "orcablocks"
            ] += f'\n%moinp "{self.conf.folder}/protocol_{self.protocol.read_orbitals}/{self.conf.number}_p{self.protocol.read_orbitals}_orca.gbw"\n'

        if "freq" in self.protocol.add_input.lower():
            calculator.parameters["orcablocks"] += "\n%freq vcd true end\n"

        return calculator, "orca"

    def single_point(self):
        return self._std_calc()

    def optimisation(self):
        calc, label = self._std_calc()
        calc.parameters["orcasimpleinput"] += " opt"
        if self.constrains:
            text = "\n%geom Constraints "
            for i in self.constrains:
                text += " {C " + str(i) + " C}"
            text += "end end\n"
            calc.parameters["orcasimpleinput"] += text
            
        if self.protocol.freq:
            calc.parameters["orcasimpleinput"] += " freq"
            if self.VERSION > 5:
                calc.parameters["orcablocks"] += "\n%freq vcd true end\n"

        return calc, label

    def frequency(self):
        calc, label = self._std_calc()
        calc.parameters["orcasimpleinput"] += " freq"
        if self.VERSION > 5:
            calc.parameters["orcablocks"] += "\n%freq vcd true end\n"
        return calc, label