import os
from ase.calculators.gaussian import Gaussian
from ensemble_analyzer._calculators.base import BaseCalc, register_calculator


@register_calculator("gaussian")
class GaussianCalc(BaseCalc):
    """
    ASE-compatible Gaussian calculator wrapper.
    This class follows the same structure as OrcaCalc,
    encapsulating the logic for SP, OPT, and FREQ calculations.
    """

    def common_str(self):
        """
        Build the common Gaussian input string (route section)
        based on the Protocol parameters.
        """

        # Construct solvent model (if any)
        solv = ""
        if self.protocol.solvent:
            if self.protocol.solvent.smd:
                solv = f" SCRF=(SMD,Solvent={self.protocol.solvent.solvent})"
            else:
                solv = f" SCRF=(CPCM,Solvent={self.protocol.solvent.solvent})"

        # Basic route section
        route = f"# {self.protocol.functional}/{self.protocol.basis}{solv}"

        # Add user-specified custom input
        if self.protocol.add_input.strip():
            route += " " + self.protocol.add_input.strip()

        if self.protocol.read_orbitals:
            route += " guess=read"

        # Default options: no population output
        # route += " pop=none"

        return route

    def _std_calc(self):
        route = self.common_str()

        calc = Gaussian(
            label="gaussian",
            output_type='N',
            mem=f"{self.cpu*2}GB",
            extra=route,
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            nprocshared=self.cpu,
        )
        if self.protocol.read_orbitals:
            calc.oldchk = (
                f"{self.conf.folder}/protocol_{self.protocol.read_orbitals}.chk"
            )

        return calc, "gaussian"

    def single_point(self):

        calc, label = self._std_calc()
        return calc, label

    def optimisation(self):

        calc, label = self._std_calc()
        if not self.protocol.constrains:
            calc.parameters["extra"] += " opt"
        else:
            calc.parameters["extra"] += " opt=(modredudant)"
            redundant = "\n".join([f"X {i+1} F" for i in self.protocol.constrains])
            # Counting in gaussian starts at 1
            
            if calc.parameters.get("addsec"):
                calc.parameters["addsec"] += redundant
            else: 
                calc.parameters["addsec"] = redundant

        if self.protocol.freq: 
            calc.parameters["extra"] += " freq=(HPModes,vcd)"

        return calc, label

    def frequency(self):

        calc, label = self._std_calc()
        calc.parameters["extra"] += " freq=(HPModes,vcd)"
        return calc, label
