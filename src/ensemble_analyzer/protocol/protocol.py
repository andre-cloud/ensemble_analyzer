
from dataclasses import dataclass, field
from typing import Union, Optional, List, Dict, Literal, Any
import json
from importlib.resources import files
from ensemble_analyzer.conformer.conformer import Conformer
from ensemble_analyzer.protocol.solvent import Solvent
from ensemble_analyzer._managers.prism_backend import ALL_PRUNE_METHODS
from pathlib import Path


def sort_protocols(protocols: List['Protocol']) -> List['Protocol']:
    return sorted(protocols, key=lambda x: x.number)



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

# ML calculators whose models have a single level of theory: `task` is an error here
TASKLESS_CALCULATORS = {"aimnet"}

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
    task                        : Optional[str]                 = None

    # Calculation settings
    mult                        : int                           = 1
    charge                      : int                           = 0
    opt                         : Optional[bool]                = False
    ts                          : Optional[bool]                = False
    freq                        : Optional[bool]                = False
    freq_fact                   : Optional[float]               = 1
    fmax                        : Optional[float]               = 0.01
    maxstep                     : Optional[float]               = 0.2
    maxiter                     : int                           = 100000000
    constraints                  : Optional[list]                = field(default_factory=list)
    read_orbitals               : Optional[str]                 = ""
    add_input                   : Optional[str]                 = ""

    # TD-DFT settings
    nroots                      : Optional[int]                 = None
    tda                         : Optional[bool]                = False
    
    # Pruning & Clustering
    no_prune                    : Optional[bool]                = False
    cluster                     : Optional[Union[bool,int]]     = False
    prune_method                : str                           = "rotconst"   # 'rotconst' | 'prism' | 'prism_rotcorr'

    # PRISM Pruner settings (used only when prune_method starts with 'prism')
    prism_moi                   : bool                          = False   # moment-of-inertia stage (opt-in: 1 % merges distinct conformers of large molecules)
    prism_moi_dev               : float                         = 1e-2    # relative MOI deviation (0.01 = 1 %)
    prism_max_rmsd              : float                         = 0.25    # [Å]
    prism_max_dev               : Optional[float]               = None    # [Å], default 2*prism_max_rmsd
    prism_timeout               : int                           = 60      # [s] per stage

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

    # Line validators [[pattern, expected, threshold], ...]
    validators                  : Optional[list]                = field(default_factory=list)

    # TS mode analysis
    loc_freq                    : Optional[list[list[int]]]     = None
    min_localization            : Optional[float]               = 40.0
    auto_displace               : Optional[bool]                = False
    displace_scale              : Optional[float]               = 0.3
    neg_freq_threshold          : Optional[float]               = 20.0
    

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
    def thr(self) -> str:
        """Formatted string of pruning thresholds.

        Returns:
            str: Multi-line string with thrG, thrB, and thrGMAX values.
        """
        out = (
            f"\tthrG    : {self.thrG} kcal/mol\n"
            f"\tthrB    : {self.thrB} cm-1\n"
            f"\tthrGMAX : {self.thrGMAX} kcal/mol\n"
        )
        if self.prune_method != "rotconst":
            max_dev = self.prism_max_dev if self.prism_max_dev is not None else 2 * self.prism_max_rmsd
            out += (
                f"\tmethod  : {self.prune_method} "
                f"(MOI {'on, ' + str(self.prism_moi_dev) if self.prism_moi else 'off'}; "
                f"RMSD < {self.prism_max_rmsd} Å, max dev < {max_dev} Å; thrB unused)\n"
            )
        return out
    
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

    def get_calculator(self, cpu: int, conf: Conformer, **kwargs) -> Any:
        """
        Instantiate the appropriate ASE calculator wrapper.

        Args:
            cpu (int): Number of CPUs to assign to the calculator.
            conf (Conformer): The conformer to be calculated.
            **kwargs: Extra keyword arguments forwarded to the calculator constructor.

        Returns:
            Any: Configured calculator instance (e.g. OrcaCalc).

        Raises:
            ValueError: If the calculator type is not registered.
        """

        from ensemble_analyzer.calculators import CALCULATOR_REGISTRY

        calc_name = self.calculator.lower()
        if calc_name not in CALCULATOR_REGISTRY:
            raise ValueError(
                f"Calculator '{calc_name}' not yet registered. "
                f"Availables: {list(CALCULATOR_REGISTRY.keys())}"
            )

        calc_class = CALCULATOR_REGISTRY[calc_name]
        calc_instance = calc_class(self, cpu, conf, **kwargs)
        
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
    
    def __repr__(self) -> str:
        """String representation of the protocol step.

        Returns:
            str: e.g. "wB97X-D4rev/def2-QZVPPD" or "r2SCAN-3c/def2-mTZVPP [water]".
        """
        if self.solvent:
            return f"{self.functional}/{self.basis} [{self.solvent}]"
        return f"{self.functional}/{self.basis}"


    # === 
    # Initialization
    # ===

    def __post_init__(self) -> None:
        """Post-initialization: validate multiplicity, resolve composite methods, load solvent and thresholds."""
        assert (self.mult > 0 and isinstance(self.mult, int)), \
            f"Multiplicity must be greater than 0, given {self.mult}"

        if self.ts and not self.opt:
            self.opt = True

        if self.functional.lower() in COMPOSITE_METHODS:
            self.basis = COMPOSITE_METHODS[self.functional.lower()]
        if 'xtb' in self.functional.lower():
            self.basis = ""

        # Load solvent
        if isinstance(self.solvent, dict) and self.solvent.get("solvent", None): 
            self.solvent = Solvent(**self.solvent)
        else: 
            self.solvent = None

        # ML task/head: empty string means "engine default"
        if isinstance(self.task, str):
            self.task = self.task.strip() or None
        if self.task and self.calculator.lower() in TASKLESS_CALCULATORS:
            raise ValueError(
                f"Protocol {self.number}: calculator '{self.calculator}' has no task/head, "
                f"but task={self.task!r} was given. Remove 'task' from this step."
            )

        # Similarity-pruning backend
        self.prune_method = (self.prune_method or "rotconst").strip().lower()
        if self.prune_method not in ALL_PRUNE_METHODS:
            raise ValueError(
                f"Protocol {self.number}: prune_method={self.prune_method!r} not valid. "
                f"Choose one of {ALL_PRUNE_METHODS}."
            )
        if self.prism_max_rmsd <= 0 or self.prism_moi_dev <= 0 or self.prism_timeout <= 0:
            raise ValueError(
                f"Protocol {self.number}: prism_max_rmsd, prism_moi_dev and prism_timeout must be > 0."
            )
        if self.prism_max_dev is not None and self.prism_max_dev <= 0:
            raise ValueError(f"Protocol {self.number}: prism_max_dev must be > 0 (or null).")

        # Clean additional input: quote compat + old-format double-brace unescape
        self.add_input = self.add_input.replace("'", "\"").replace("{{", "{").replace("}}", "}")

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
    with open(default if not file else file) as f:
        return json.load(f)