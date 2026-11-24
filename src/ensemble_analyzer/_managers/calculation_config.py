from dataclasses import dataclass
from typing import Any, Optional, Dict, List
from pathlib import Path
import json

from ensemble_analyzer.protocol.protocol import Protocol

@dataclass
class CalculationConfig:
    """
    Configuration container for global calculation settings.
    
    Stores runtime parameters that are constant across all protocols
    (e.g., resources, physical constants, global flags).
    """

    cpu: int = 1 
    temperature: float = 298.15
    start_from_protocol: int = 0
    include_H: bool = False
    restart: bool = False

    # Graph settings
    definition: int = 4
    fwhm: Optional[Dict[str, Optional[float]]] = None
    shift: Optional[Dict[str, Optional[float]]] = None
    interested: Optional[Dict[str, Optional[float]]] = None
    invert: bool = False

    linear: bool = False
    cut_off: float = 100
    alpha: float = 4
    P: float = 101.325
        
    def __post_init__(self) -> None:
        """
        Initialize default values for optional graph settings.

        Ensures that fwhm, shift, and interested dictionaries are always
        populated with their default 'vibro'/'electro' keys.
        """
        if self.fwhm is None:
            self.fwhm = {'vibro': None, 'electro': None}
        if self.shift is None:
            self.shift = {'vibro': None, 'electro': None}
        if self.interested is None:
            self.interested = {'vibro': None, 'electro': None}


    @classmethod
    def from_args(cls, args, start_from_protocol: int = 0) -> 'CalculationConfig':
        """
        Create a config from command-line arguments.

        If `settings.json` exists its values win (keys missing from it fall back
        to the arguments); otherwise it is created from the arguments.

        Args:
            args (argparse.Namespace): Parsed command-line arguments.
            start_from_protocol (int): Protocol ID to resume from.

        Returns:
            CalculationConfig: Initialized configuration object.
        """
        settings = cls._args_to_dict(args)
        settings_file = Path("settings.json")
        if settings_file.exists():
            with open(settings_file) as f:
                settings.update(json.load(f))
        else:
            with open(settings_file, 'w') as f:
                json.dump(settings, f, indent=4)
        return cls._from_settings(settings, start_from_protocol)

    @classmethod
    def load(cls, filepath: Path = Path("settings.json")) -> 'CalculationConfig':
        """
        Load configuration from a settings JSON file.

        Raises:
            FileNotFoundError: If the settings file is missing.
        """
        with open(filepath) as f:
            return cls._from_settings(json.load(f))

    @classmethod
    def _from_settings(cls, s: dict, start_from_protocol: int = 0) -> 'CalculationConfig':
        """Build a config from a settings dict (keys as in settings.json)."""
        def pair(key):
            return {'vibro': s.get(f"{key}_vibro"), 'electro': s.get(f"{key}_electro")}
        simple = ("cpu", "temperature", "include_H", "definition", "invert",
                  "restart", "linear", "cut_off", "alpha", "P")
        return cls(
            start_from_protocol=start_from_protocol,
            fwhm=pair("fwhm"), shift=pair("shift"), interested=pair("interested"),
            **{k: s[k] for k in simple if k in s},
        )

    @staticmethod
    def _args_to_dict(args: Any) -> dict:
        """Convert the argparse namespace to the settings.json dictionary."""
        return {
            "cpu": args.cpu,
            "temperature": args.temperature,
            "definition": args.definition,
            "fwhm_vibro": args.fwhm_vibro,
            "fwhm_electro": args.fwhm_electro,
            "shift_vibro": args.shift_vibro,
            "shift_electro": args.shift_electro,
            "interested_vibro": args.interest_vibro,
            "interested_electro": args.interest_electro,
            "invert": args.invert,
            "include_H": args.include_H,
            "restart": args.restart,
            "P": args.pressure,
            "linear": args.linear,
            "alpha": args.alpha,
            "cut_off": args.cut_off,
        }

    def create_log(self, protocols: List[Protocol], conformers: int) -> dict:
        """
        Create a summary dictionary of the calculation setup.

        Args:
            protocols (List[Protocol]): The list of protocols to be executed.
            conformers (int): Number of conformers in the ensemble.

        Returns:
            dict: Summary dictionary with setup details.
        """
        return {
            'conformers' : conformers, 
            'protocols' : protocols, 
            'len_protocols' : len(protocols),
            'temperature' : self.temperature, 
            'cpu' : self.cpu, 
            'restart' : self.restart,
        }