import json
import os

from src._calculators.base import CALCULATOR_REGISTRY
from typing import Union, List, Optional

DEBUG = os.getenv("DEBUG")
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
    INTERNALS = {2:'B', 3:'A', 4:'D'}

    def __init__(
        self,
        number: int,
        functional: str,
        basis: Optional[str] = "",
        solvent: Optional[dict] = {},
        opt: Optional[bool] = False,
        freq: Optional[bool] = False,
        add_input: Optional[str] = "",
        freq_fact: Optional[float] = 1,
        mult: int = 1,
        charge: int = 0,
        graph: Optional[bool] = False,
        calculator: str = "orca",
        thrG: Optional[float] = None,
        thrB: Optional[float] = None,
        thrGMAX: Optional[float] = None,
        # thrRMSD_enantio : float = None,
        constrains: Optional[list] = [],
        cluster: Optional[bool|int] = False,
        no_prune: Optional[bool] = False,
        comment: Optional[str] = "",
        read_orbitals: Optional[str] = "",
        read_population: Optional[str|None] = None,
        monitor_internals: Optional[List[List[int]]] = [],
        skip_opt_fail: Optional[bool] = False
    ): 
        self.number = number
        self.functional = functional.upper()
        self.basis = basis.upper() if 'xtb' not in functional.lower() else ""
        self.solvent = Solvent(solvent) if solvent.get('solvent', None) else None
        self.opt = opt
        self.freq = freq
        self.add_input = add_input.replace("'", '"')
        self.thrG = thrG
        self.thrB = thrB
        self.thrGMAX = thrGMAX
        # self.thrRMSD_enantio = thrRMSD_enantio
        self.get_thrs(self.load_threshold())
        self.calculator = calculator
        self.constrains = constrains
        self.cluster = cluster
        self.no_prune = no_prune
        self.mult = mult
        self.charge = charge
        self.comment = comment
        self.read_orbitals = read_orbitals # number of protocol to read orbitals from
        self.read_population = read_population
        self.monitor_internals = monitor_internals
        self.skip_opt_fail = skip_opt_fail
        
        assert self.mult > 0, "Multiplicity must be greater than 0"

        self.freq_fact = freq_fact
        self.graph = graph

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

    def get_calculator(self, cpu, conf=None):
        """
        Get the calculator from the user selector

        :param cpu: allocated CPU
        :type cpu: int
        :param mode: type of calculation required. Choose between: opt, freq, energy
        :type mode: str
        :param conf: Conformer instance, if needed
        :type conf: Conformer
        """

        calc_name = self.calculator.lower()
        if calc_name not in CALCULATOR_REGISTRY:
            raise ValueError(f"Calculator '{calc_name}' not yet registered. "
                            f"Availables: {list(CALCULATOR_REGISTRY.keys())}")

        calc_class = CALCULATOR_REGISTRY[calc_name]
        calc_instance = calc_class(self, cpu, conf)

        mode_map = {
            "opt": calc_instance.optimisation,
            "freq": calc_instance.frequency,
            "energy": calc_instance.single_point,
        }
     
        if self.opt: 
            return mode_map["opt"]()
        if self.freq: 
            return mode_map["freq"]()
        return mode_map["energy"]()


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
        # if not self.thrRMSD_enantio:
        #     self.thrRMSD_enantio = thr_json[c]["thrRMSD_enantio"]

    def verbal_internals(self): 
        internals = []
        for internal in self.monitor_internals: 
            internals.append(f'{self.INTERNALS[len(internal)]} {"-".join([str(i) for i in internal])}')
        return internals
    
    ## Properties

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


    ## Repr methods

    def __str__(self):  # pragma: no cover
        if self.solvent:
            return f"{self.functional}/{self.basis} - {self.solvent}"
        return f"{self.functional}/{self.basis}"

    def __repr__(self):  # pragma: no cover
        if self.solvent:
            return f"{self.functional}/{self.basis} - {self.solvent}"
        return f"{self.functional}/{self.basis}"

    ## Statics

    @staticmethod
    def load_raw(json):
        return Protocol(**json)
