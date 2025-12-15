import os
from ase.calculators.gaussian import Gaussian
from ensemble_analyzer._calculators.base import BaseCalc, register_calculator

from typing import Tuple

@register_calculator("gaussian")
class GaussianCalc(BaseCalc):
    """
    ASE-compatible Gaussian calculator wrapper.
    Encapsulates logic for SP, OPT, and FREQ input generation.
    """

    def common_str(self)-> str:
        """
        Build the common Gaussian route section.

        Returns:
            str: Route string (e.g. '# B3LYP/6-31G* SCRF=...').
        """

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

        return route

    def _std_calc(self) -> Tuple[Gaussian, str]:
        """
        Create standard Gaussian calculator.

        Returns:
            Tuple[Gaussian, str]: Initialized calculator and label.
        """

        route = self.common_str()

        calc = Gaussian(
            label="gaussian",
            output_type='N',
            mem=f"{self.cpu*2}GB",
            chk='gaussian.chk',
            extra=route,
            charge=self.protocol.charge,
            mult=self.protocol.mult,
            nprocshared=self.cpu,
        )
        if self.protocol.read_orbitals:
            calc.oldchk = (
                f"{self.conf.folder}/protocol_{self.protocol.read_orbitals}/{self.conf.number}_{self.protocol.read_orbitals}_gaussian.chk"
            )

        return calc, "gaussian"

    def single_point(self) -> Tuple[Gaussian, str]:
        """Configure Single Point calculation."""

        calc, label = self._std_calc()
        return calc, label

    def optimisation(self) -> Tuple[Gaussian, str]:
        """
        Configure Geometry Optimization.
        Handles modredundant constraints.
        """

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

    def frequency(self) -> Tuple[Gaussian, str]:
        """Configure Frequency calculation."""

        calc, label = self._std_calc()
        calc.parameters["extra"] += " freq=(HPModes,vcd)"
        return calc, label
