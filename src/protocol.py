import json
import os
import sys
from ase.calculators.orca import ORCA
from ase.calculators.orca import OrcaProfile


DEBUG = os.getenv("DEBUG")


orca_profile = OrcaProfile(command='/opt/orca/6.0.1/orca')


def load_protocol(file: str):  # pragma: no cover
    default = "ensemble_analyser/parameters_file/default_protocol.json"
    return json.load(open(default if not file else file))


LEVEL_DEFINITION = {
    0: "SP".lower(),  # mere energy calculation
    1: "OPT".lower(),  # optimisation step
    2: "FREQ".lower(),  # single point and frequency analysis
    3: "OPT+FREQ".lower(),  # optimisation and frequency analysis
}


class Solvent:
    """
    Solvent class
    """

    def __init__(self, solv: dict):
        self.solvent = solv["solvent"]
        self.smd = solv["smd"]

    def __str__(self):  # pragma: no cover
        if self.smd:
            return f"SMD({self.solvent})"
        elif self.solvent:
            return f"CPCM({self.solvent})"
        else:
            return "CPCM"

    def __repr__(self):  # pragma: no cover
        if self.smd:
            return f"SMD({self.solvent})"
        elif self.solvent:
            return f"CPCM({self.solvent})"
        else:
            return "CPCM"

    def orca_input_smd(self):
        """Get the ORCA input string for the SMD

        :return: SMD input for ORCA
        :rtype: str
        """
        if self.smd:
            return f'%cpcm smd true smdsolvent "{self.solvent}" end'
        return ""


class Protocol:
    def __init__(
        self,
        number: int,
        functional: str,
        basis: str = "def2-svp",
        solvent: dict = {},
        opt: bool = False,
        freq: bool = False,
        add_input: str = "",
        freq_fact: float = 1,
        graph: bool = False,
        calculator: str = "orca",
        thrG: float = None,
        thrB: float = None,
        thrGMAX: float = None,
        constrains: list = [],
        maxstep: float = 0.2,
        fmax: float = 0.05,
        cluster: bool = False,
    ):
        self.number = number
        self.functional = functional.upper()
        self.basis = basis.upper()
        self.solvent = Solvent(solvent) if solvent else None
        self.opt = opt
        self.freq = freq
        self.add_input = add_input.replace("'", '"')
        self.thrG = thrG
        self.thrB = thrB
        self.thrGMAX = thrGMAX
        self.get_thrs(self.load_threshold())
        self.calculator = calculator
        self.constrains = constrains
        self.maxstep = maxstep
        self.cluster = cluster

        if fmax != 0.05:
            self.fmax = fmax
        elif self.constrains:
            self.fmax = 0.1
        else:  # if an opt freeze, the max force convergence will be lifted.
            self.fmax = fmax

        self.freq_fact = freq_fact
        self.graph = graph

    @property
    def calculation_level(self):
        return LEVEL_DEFINITION[self.number_level].upper()

    @property
    def level(self):
        return f"{self.functional}/{self.basis}" + (
            str(self.solvent) if self.solvent else ""
        )

    @property
    def thr(self):
        return f"\tthrG    : {self.thrG} kcal/mol\n\tthrB    : {self.thrB} cm-1\n\tthrGMAX : {self.thrGMAX} kcal/mol\n"

    @property
    def number_level(self):
        c = 0
        if self.opt:
            c += 1
        if self.freq:
            c += 2
        return c

    def load_threshold(self) -> dict:
        """
        Load default thresholds

        :return: thresholds
        :rtype: dict
        """

        default = os.path.join(
            os.path.abspath(os.path.dirname(__file__)),
            "parameters_file",
            "default_threshold.json",
        )
        return json.load(open(default))

    def get_calculator(self, cpu, charge: int, mult: int, mode: str):
        """
        Get the calculator from the user selector

        :param cpu: allocated CPU
        :type cpu: int
        :param charge: charge of the molecule
        :type charge: int
        :param mult: multiplicity of the molecule
        :type mult: int
        :param mode: type of calculation required. Choose between: opt, freq, energy
        :type mode: str
        """

        calc = {
            "orca": {
                "opt": self.orca_opt,
                "freq": self.orca_freq,
                "energy": self.calc_orca_std,
            },
        }

        return calc[self.calculator][mode](cpu, charge, mult)

    def get_thrs(self, thr_json):
        """
        Get default thrs if not defined by user

        :param thr_json: JSON default thresholds
        :type thr_json: dict
        """
        c = LEVEL_DEFINITION[self.number_level]
        if not self.thrG:
            self.thrG = thr_json[c]["thrG"]
        if not self.thrB:
            self.thrB = thr_json[c]["thrB"]
        if not self.thrGMAX:
            self.thrGMAX = thr_json[c]["thrGMAX"]

    def __str__(self):  # pragma: no cover
        if self.solvent:
            return f"{self.functional}/{self.basis} - {self.solvent}"
        return f"{self.functional}/{self.basis}"

    def __repr__(self):  # pragma: no cover
        if self.solvent:
            return f"{self.functional}/{self.basis} - {self.solvent}"
        return f"{self.functional}/{self.basis}"

    # ORCA CALCULATOR

    def orca_common_str(self, cpu):
        """ORCA common string for the input

        :param cpu: number of CPU
        :type cpu: int
        :return: input string
        :rtype: str
        """

        if self.solvent and not (self.solvent.smd and self.freq):
            if "xtb" in self.functional.lower():
                solv = f"ALPB({self.solvent.solvent})"
            elif self.solvent.solvent.strip():
                solv = f" CPCM({self.solvent.solvent.strip()})"
            else:
                solv = f" CPCM"
        else:
            solv = ""

        si = f"{self.functional} {self.basis} {solv} nopop"

        smd = ""
        if self.solvent and ("xtb" not in self.functional.lower()) and (not self.freq):
            smd = self.solvent.orca_input_smd()

        ob = (
            f"%pal nprocs {cpu} end "
            + smd
            + self.add_input
            + (" %maxcore 6000" if "maxcore" not in self.add_input else "")
        )

        return si, ob

    def calc_orca_std(self, cpu: int, charge: int, mult: int):
        """Standard calculator

        :param cpu: CPU number
        :type cpu: int
        :param charge: charge of the molecule
        :type charge: int
        :param mult: multiplicity of the molecule
        :type mult: int
        :return: calculator and label
        :rtype: tuple
        """

        simple_input, ob = self.orca_common_str(cpu)
        label = "ORCA"
        calculator = ORCA(
            profile=orca_profile,
            label=label,
            orcasimpleinput=simple_input,
            orcablocks=ob,
            charge=charge,
            mult=mult,
        )

        return calculator, label

    def orca_opt(self, cpu: int, charge: int, mult: int):
        """Optimization calculator

        :param cpu: CPU number
        :type cpu: int
        :param charge: charge of the molecule
        :type charge: int
        :param mult: multiplicity of the molecule
        :type mult: int
        :return: calculator and label
        :rtype: tuple
        """
        calculator, label = self.calc_orca_std(cpu, charge, mult)
        calculator.parameters["orcasimpleinput"] += " engrad"

        return calculator, label

    def orca_freq(self, cpu: int, charge: int, mult: int):
        """Frequency calculator

        :param cpu: CPU number
        :type cpu: int
        :param charge: charge of the molecule
        :type charge: int
        :param mult: multiplicity of the molecule
        :type mult: int
        :return: calculator and label
        :rtype: tuple
        """
        calculator, label = self.calc_orca_std(cpu, charge, mult)
        calculator.parameters["orcasimpleinput"] += " freq"

        return calculator, label

    @staticmethod
    def load_raw(json):
        return Protocol(**json)
