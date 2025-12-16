
from dataclasses import dataclass, field
from typing import Union, Optional, List, Dict, Literal, Any
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
    Configuration for a single computational step.

    Stores all parameters required to setup a QM calculation (method, basis,
    thresholds) and post-processing options (pruning, clustering).
    """

    number              : int
    
    # Method definition
    functional                  : str
    basis                       : Optional[str]                 = ""
    solvent                     : Optional[Dict]                = None
    calculator                  : str                           = "orca"
    
    # Calculation settings
    mult                        : int                           = 1
    charge                      : int                           = 0
    opt                         : Optional[bool]                = False
    freq                        : Optional[bool]                = False
    freq_fact                   : Optional[float]               = 1
    constrains                  : Optional[list]                = field(default_factory=list)
    read_orbitals               : Optional[str]                 = ""
    add_input                   : Optional[str]                 = ""
    
    # Pruning & Clustering
    graph                       : Optional[bool]                = False
    no_prune                    : Optional[bool]                = False
    cluster                     : Optional[Union[bool,int]]     = False

    # Thresholds
    thrG                        : Optional[float]               = None
    thrB                        : Optional[float]               = None
    thrGMAX                     : Optional[float]               = None

    # Logging
    monitor_internals           : Optional[List[List[int]]]     = field(default_factory=list)
    comment                     : Optional[str]                 = ""

    # Options
    read_population             : Optional[str|None]            = None
    skip_opt_fail               : Optional[bool]                = False
    block_on_retention_rate     : Optional[bool]                = False
    

    # ===
    # Thresholds
    # === 

    def load_threshold(self) -> dict:
        """
        Load default pruning thresholds from the configuration file.

        Returns:
            dict: Dictionary of thresholds keyed by calculation level (opt, freq, sp).
        """

        default = files("ensemble_analyzer").joinpath("parameters_file/default_threshold.json")

        with open(default, "r") as f:
            return json.load(f)

    def get_thrs(self, thr_json: dict) -> None:
        """
        Update instance thresholds based on the calculation level defaults.

        If specific thresholds (thrG, thrB) are not set by the user, applies
        defaults from the provided dictionary.

        Args:
            thr_json (dict): Dictionary of default thresholds.

        Returns:
            None
        """

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
    def number_level(self) -> int:
        """
        Calculate the complexity level of the protocol.

        Returns:
            int: Level identifier (0=SP, 1=OPT, 2=FREQ, 3=OPT+FREQ).
        """

        c = 0
        if self.opt:
            c += 1
        if self.freq:
            c += 2
        return c

    @property
    def calculation_level(self) -> str:
        """
        Get the string representation of the calculation type.

        Returns:
            str: e.g., "OPT", "FREQ", "SP".
        """

        return LEVEL_DEFINITION[self.number_level].upper()

    @property
    def thr(self):
        return (
            f"\tthrG    : {self.thrG} kcal/mol\n"
            f"\tthrB    : {self.thrB} cm-1\n"
            f"\tthrGMAX : {self.thrGMAX} kcal/mol\n"
        )
    
    @property
    def clustering(self) -> bool:
        """
        Check if clustering is enabled for this step.

        Returns:
            bool: True if clustering is requested.
        """

        if isinstance(self.cluster, bool):
            return self.cluster
        
        if isinstance(self.cluster, int):
            return self.cluster > 1
        
        return False

    # ===
    # Functions
    # === 

    def verbal_internals(self) -> List[str]:
        """
        Generate human-readable strings for monitored internal coordinates.

        Returns:
            List[str]: Descriptions of monitored internals (e.g. "B 0-1", "A 0-1-2").
        """

        internals = []
        for internal in self.monitor_internals:
            internals.append(f"{INTERNALS[len(internal)]} {'-'.join(str(i) for i in internal)}")
        return internals

    def get_calculator(self, cpu: int, conf:Conformer) -> Any:
        """
        Instantiate the appropriate ASE calculator wrapper.

        Args:
            cpu (int): Number of CPUs to assign to the calculator.
            conf (Conformer): The conformer to be calculated.

        Returns:
            Any: Configured calculator instance (e.g. OrcaCalc).

        Raises:
            ValueError: If the calculator type is not registered.
        """

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
    def load_raw(json) -> 'Protocol':
        """
        Factory method to create a Protocol instance from a dictionary.

        Args:
            json_data (dict): Dictionary containing protocol parameters.

        Returns:
            Protocol: Initialized protocol object.
        """
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
    """
    Load the full protocol sequence from a JSON file.

    Args:
        file (Optional[str]): Path to the protocol JSON file. If None, loads default.

    Returns:
        Dict: Dictionary mapping step numbers to protocol configurations.
    """
    
    default = files("ensemble_analyzer").joinpath("parameters_file/default_protocol.json")
    return json.load(open(default if not file else file))