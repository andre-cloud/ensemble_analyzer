

from dataclasses import dataclass
from typing import Optional, Dict, List
from pathlib import Path
import json

from ensemble_analyzer._protocol.protocol import Protocol

@dataclass
class CalculationConfig:
    """Configuration for ensemble calculations."""
    cpu: int = 1 
    temperature: float = 298.15
    start_from_protocol: int = 0
    include_H: bool = True
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
        
    def __post_init__(self):
        if self.fwhm is None:
            self.fwhm = {'vibro': None, 'electro': None}
        if self.shift is None:
            self.shift = {'vibro': None, 'electro': None}
        if self.interested is None:
            self.interested = {'vibro': None, 'electro': None}


    @classmethod
    def from_args(cls, args, start_from_protocol: int = 0) -> 'CalculationConfig':
        """
        Create configuration from parsed arguments.
        
        Also handles settings.json persistence:
        - If settings.json exists: load and merge with args
        - If not: create from args and save
        
        Args:
            args: Parsed command-line arguments
            start_from_protocol: Protocol to start from (for restart)
            
        Returns:
            CalculationConfig instance
        """
        settings_file = Path("settings.json")
        
        # Try to load existing settings
        if settings_file.exists():
            with open(settings_file) as f:
                settings = json.load(f)
        else:
            # Create settings from args
            settings = cls._args_to_dict(args)
            
            # Save for future use
            with open(settings_file, 'w') as f:
                json.dump(settings, f, indent=4)
        
        # Create config from settings
        return cls(
            cpu=settings.get("cpu", args.cpu),
            temperature=settings.get("temperature", args.temperature),
            start_from_protocol=start_from_protocol,
            include_H=settings.get("include_H", not args.exclude_H),
            definition=settings.get("definition", args.definition),
            fwhm={
                'vibro': settings.get("fwhm_vibro", args.fwhm_vibro),
                'electro': settings.get("fwhm_electro", args.fwhm_electro)
            },
            shift={
                'vibro': settings.get("shift_vibro", args.shift_vibro),
                'electro': settings.get("shift_electro", args.shift_electro)
            },
            interested={
                'vibro': settings.get("interested_vibro", args.interest_vibro),
                'electro': settings.get("interested_electro", args.interest_electro)
            },
            invert=settings.get("invert", args.invert),
            restart = settings.get("restart", args.restart)
        )
    
    @staticmethod
    def _args_to_dict(args) -> dict:
        """
        Convert argparse namespace to settings dictionary.
        
        Args:
            args: Parsed arguments
            
        Returns:
            Dictionary suitable for settings.json
        """
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
            "include_H": not args.exclude_H,
            "restart": args.restart,
        }
    
    def to_dict(self) -> dict:
        """
        Convert to dictionary for serialization.
        
        Returns:
            Dictionary representation
        """
        d = self._args_to_dict(self)
        
        # Flatten nested dicts for settings.json compatibility
        result = {
            "cpu": d["cpu"],
            "temperature": d["temperature"],
            "definition": d["definition"],
            "include_H": d["include_H"],
            "invert": d["invert"],
        }
        
        # Flatten fwhm
        if d["fwhm"]:
            result["fwhm_vibro"] = d["fwhm"].get("vibro")
            result["fwhm_electro"] = d["fwhm"].get("electro")
        
        # Flatten shift
        if d["shift"]:
            result["shift_vibro"] = d["shift"].get("vibro")
            result["shift_electro"] = d["shift"].get("electro")
        
        # Flatten interested
        if d["interested"]:
            result["interested_vibro"] = d["interested"].get("vibro")
            result["interested_electro"] = d["interested"].get("electro")
        
        return result
    
    def save(self, filepath: Path = Path("settings.json")) -> None:
        """
        Save configuration to file.
        
        Args:
            filepath: Path to settings file
        """
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=4)
    
    @classmethod
    def load(cls, filepath: Path = Path("settings.json")) -> 'CalculationConfig':
        """
        Load configuration from file.
        
        Args:
            filepath: Path to settings file
            
        Returns:
            CalculationConfig instance
            
        Raises:
            FileNotFoundError: If settings file doesn't exist
        """
        if not filepath.exists():
            raise FileNotFoundError(f"Settings file not found: {filepath}")
        
        with open(filepath) as f:
            settings = json.load(f)
        
        return cls(
            cpu=settings["cpu"],
            temperature=settings["temperature"],
            start_from_protocol=settings.get("start_from_protocol", 0),
            include_H=settings.get("include_H", True),
            definition=settings.get("definition", 4),
            fwhm={
                'vibro': settings.get("fwhm_vibro"),
                'electro': settings.get("fwhm_electro")
            },
            shift={
                'vibro': settings.get("shift_vibro"),
                'electro': settings.get("shift_electro")
            },
            interested={
                'vibro': settings.get("interested_vibro"),
                'electro': settings.get("interested_electro")
            },
            invert=settings.get("invert", False),
            restart = settings.get("restart", False)
        )
    
    def validate(self) -> None:
        """
        Validate configuration parameters.
        
        Raises:
            ValueError: If any parameter is invalid
        """
        if self.cpu < 1:
            raise ValueError(f"CPU count must be ≥ 1, got {self.cpu}")
        
        if self.temperature <= 0:
            raise ValueError(f"Temperature must be > 0 K, got {self.temperature}")
        
        if self.definition < 1:
            raise ValueError(f"Definition must be ≥ 1, got {self.definition}")
        
    def create_log(self, protocols: List[Protocol], conformers: int): 
        return {
            'conformers' : conformers, 
            'protocols' : protocols, 
            'len_protocols' : len(protocols),
            'temperature' : self.temperature, 
            'cpu' : self.cpu, 
            'restart' : self.restart,
        }