
from dataclasses import dataclass, field
from typing import Union, Optional, List, Dict, Literal
import json

from importlib.resources import files

from ensemble_analyzer._conformer.conformer import Conformer
from ensemble_analyzer._calculators import CALCULATOR_REGISTRY

from ensemble_analyzer._protocol.solvent import Solvent

from pathlib import Path



LEVEL_DEFINITION = {
    0: "SP".lower(),        # Single point calculation
    1: "OPT".lower(),       # Optimization Step
    2: "FREQ".lower(),      # Frequency analysis
    3: "OPT+FREQ".lower(),  # Optimization and Frequency
}

COMPOSITE_METHODS = {
    "xTB"       .lower() : "",
    "HF-3c"     .lower() : "MINIX",
    "B97-3c"    .lower() : "def2-mTZVP",
    "PBEh-3c"   .lower() : "def2-mSVP",
    "r2SCAN-3c" .lower() : "def2-mTZVPP", 
    "wB97X-3c"  .lower() : "vDZP"
}

INTERNALS = {2: 'B', 3: 'A', 4: 'D'}

@dataclass
class Protocol:
    """
    Protocol generation
    """

    number              : int
    
    # Method definition
    functional          : str
    basis               : Optional[str]                 = ""
    solvent             : Optional[Dict]                = None
    calculator          : str                           = "orca"
    
    # Calculation settings
    mult                : int                           = 1
    charge              : int                           = 0
    opt                 : Optional[bool]                = False
    freq                : Optional[bool]                = False
    freq_fact           : Optional[float]               = 1
    constrains          : Optional[list]                = field(default_factory=list)
    read_orbitals       : Optional[str]                 = ""
    add_input           : Optional[str]                 = ""
    
    # Pruning & Clustering
    graph               : Optional[bool]                = False
    no_prune            : Optional[bool]                = False
    cluster             : Optional[Union[bool,int]]     = False

    # Thresholds
    thrG                : Optional[float]               = None
    thrB                : Optional[float]               = None
    thrGMAX             : Optional[float]               = None

    # Logging
    monitor_internals   : Optional[List[List[int]]]     = field(default_factory=list)
    comment             : Optional[str]                 = ""

    # Options
    read_population     : Optional[str|None]            = None
    skip_opt_fail       : Optional[bool]                = False
    skip_retention_rate : Optional[bool]                = False
    

    # ===
    # Thresholds
    # === 

    def load_threshold(self) -> dict:
        default = files("ensemble_analyzer").joinpath("parameters_file/default_threshold.json")

        with open(default, "r") as f:
            return json.load(f)

    def get_thrs(self, thr_json: dict):
        c = LEVEL_DEFINITION[self.number_level]
        if self.thrG is None:
            self.thrG = thr_json[c]["thrG"]
        if self.thrB is None:
            self.thrB = thr_json[c]["thrB"]
        if self.thrGMAX is None:
            self.thrGMAX = thr_json[c]["thrGMAX"]

    # === 
    # Properties
    # ===

    @property
    def number_level(self):
        c = 0
        if self.opt:
            c += 1
        if self.freq:
            c += 2
        return c

    @property
    def calculation_level(self):
        return LEVEL_DEFINITION[self.number_level].upper()

    @property
    def thr(self):
        return (
            f"\tthrG    : {self.thrG} kcal/mol\n"
            f"\tthrB    : {self.thrB} cm-1\n"
            f"\tthrGMAX : {self.thrGMAX} kcal/mol\n"
        )
    
    @property
    def clustering(self):
        if isinstance(self.cluster, bool):
            return self.cluster
        
        if isinstance(self.cluster, int):
            return self.cluster > 1
        
        return False

    # ===
    # Functions
    # === 

    def verbal_internals(self):
        internals = []
        for internal in self.monitor_internals:
            internals.append(f"{INTERNALS[len(internal)]} {'-'.join(str(i) for i in internal)}")
        return internals

    def get_calculator(self, cpu: int, conf:Conformer):

        calc_name = self.calculator.lower()
        if calc_name not in CALCULATOR_REGISTRY:
            raise ValueError(
                f"Calculator '{calc_name}' not yet registered. "
                f"Availables: {list(CALCULATOR_REGISTRY.keys())}"
            )

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
    
    # ===
    # Static Functions
    # ===

    @staticmethod
    def load_raw(json):
        return Protocol(**json)
    
    def __repr__(self): 
        if self.solvent:
            return f"{self.functional}/{self.basis} [{self.solvent}]"
        return f"{self.functional}/{self.basis}"


    # === 
    # Initialization
    # ===

    def __post_init__(self):

        assert (self.mult > 0 and isinstance(self.mult, int)), \
            f"Multiplicity must be greater than 0, given {self.mult}"

        if self.functional.lower() in COMPOSITE_METHODS:
            self.basis = COMPOSITE_METHODS[self.functional.lower()]
        if 'xtb' in self.functional.lower():
            self.basis = ""

        # Load solvent
        if isinstance(self.solvent, dict) and self.solvent.get("solvent", None): 
            self.solvent = Solvent(**self.solvent)
        else: 
            self.solvent = None

        # Clean additional input
        self.add_input = self.add_input.replace("'", "\"")

        # Load eventual more Thresholds
        self.get_thrs(self.load_threshold())



def load_protocol(file: Optional[str]) -> Dict: 
    """Load protocol from JSON files

    Args:
        file (Optional[str]): Protocol Filename

    Returns:
        Dict: Dictionary with all the settings defined by the user
    """
    
    default = files("ensemble_analyzer").joinpath("parameters_file/default_protocol.json")
    return json.load(open(default if not file else file))