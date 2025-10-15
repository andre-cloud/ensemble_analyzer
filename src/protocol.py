import json
import os
import sys
from ase.calculators.orca import ORCA

try:
    from ase.calculators.orca import OrcaProfile
except ImportError:
    pass


DEBUG = os.getenv("DEBUG")

try:
    orca_profile = OrcaProfile(command="/opt/orca/6.1.0/orca")
except TypeError:
    orca_profile = None
except NameError:
    orca_profile = None


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
        mult: int = 1,
        charge: int = 0,
        graph: bool = False,
        calculator: str = "orca",
        thrG: float = None,
        thrB: float = None,
        thrGMAX: float = None,
        constrains: list = [],
        maxstep: float = 0.2,
        fmax: float = 0.05,
        cluster: bool|int = False,
        no_prune: bool = False,
        comment: str = "",
        read_orbitals = "",
        read_population: str|None = None,
    ): 
        self.number = number
        self.functional = functional.upper()
        self.basis = basis.upper() if 'xtb' not in functional.lower() else ""
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
        self.no_prune = no_prune
        self.mult = mult
        self.charge = charge
        self.comment = comment
        self.read_orbitals = read_orbitals # number of protocol to read orbitals from
        self.read_population = read_population
        
        assert self.mult > 0, "Multiplicity must be greater than 0"


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
        return LEVEL_DEFINITION[self.number_level].upper() + (f'-> {self.comment}' if self.comment else "")

    @property
    def level(self):
        return f"{self.functional}/{self.basis}" + (
            ("["+str(self.solvent))+"]" if self.solvent else ""
        ) + (f'-> {self.comment}' if self.comment else "")

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

    def get_calculator(self, cpu, mode: str, conf=None):
        """
        Get the calculator from the user selector

        :param cpu: allocated CPU
        :type cpu: int
        :param mode: type of calculation required. Choose between: opt, freq, energy
        :type mode: str
        :param conf: Conformer instance, if needed
        :type conf: Conformer
        """

        calc = {
            "orca": {
                "opt": self.orca_opt,
                "freq": self.orca_freq,
                "energy": self.calc_orca_std,
            },
        }

        return calc[self.calculator][mode](cpu, conf)

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

        if self.solvent:
            if "xtb" in self.functional.lower():
                solv = f"ALPB({self.solvent.solvent})"
            elif self.solvent.solvent.strip():
                solv = f" {self.solvent}"
            else:
                solv = f" CPCM"
        else:
            solv = ""

        si = f"{self.functional} {self.basis} {solv} nopop"


        ob = (
            f"%pal nprocs {cpu} end "
            + self.add_input
            + (" %maxcore 5000" if "maxcore" not in self.add_input else "")
        )

        return si, ob

    def calc_orca_std(self, cpu: int, conf=None):
        """Standard calculator

        :param cpu: CPU number
        :type cpu: int
        :return: calculator and label
        :rtype: tuple
        """

        simple_input, ob = self.orca_common_str(cpu)
        label = "orca"
        calculator = ORCA(
            profile=orca_profile,
            label=label,
            orcasimpleinput=simple_input,
            orcablocks=ob,
            charge=self.charge,
            mult=self.mult,
        )
        if self.read_orbitals:
            calculator.parameters['orcasimpleinput'] += " moread"
            calculator.parameters["orcablocks"] += f'\n%moinp "{conf.folder}/protocol_{self.read_orbitals}.gbw"\n'

        if 'freq' in self.add_input.lower(): 
            calculator.parameters["orcablocks"] += "\n%freq vcd true end\n"

        return calculator, label

    def orca_opt(self, cpu: int, conf=None):
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
        calculator, label = self.calc_orca_std(cpu, conf)
        calculator.parameters["orcasimpleinput"] += " engrad"

        return calculator, label

    def orca_freq(self, cpu: int, conf=None):
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
        calculator, label = self.calc_orca_std(cpu, conf)
        calculator.parameters["orcasimpleinput"] += " freq"
        calculator.parameters["orcablocks"] += "\n%freq vcd true end\n"

        return calculator, label

    @staticmethod
    def load_raw(json):
        return Protocol(**json)
