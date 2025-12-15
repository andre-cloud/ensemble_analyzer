from dataclasses import dataclass, field, asdict
from typing import Optional, Dict, Tuple, Union
import numpy as np



@dataclass
class EnergyRecord:
    """
    Data container for energetic and thermodynamic properties of a conformer.
    """

    E       : float                     = 0.0           # Electronic Energy [Eh]
    G       : float                     = np.nan        # Gibbs Free Energy [Eh]
    H       : float                     = np.nan        # Enthalpy [Eh]
    S       : float                     = np.nan        # Total Entropy [Eh]
    G_E     : float                     = np.nan        # Thermal Correction to Gibbs (G-E)
    zpve    : float                     = np.nan        # Zero Point Vibrational Energy
    B       : Optional[float]           = None          # Rotational Constant Norm [cm-1]
    B_vec   : Optional[np.ndarray]      = None          # Rotational Constants Vector [cm-1]
    m       : Optional[float]           = None          # Dipole Moment Norm [Debye]
    m_vec   : Optional[np.ndarray]      = None          # Dipole Moment Vector
    Pop     : float                     = np.nan        # Boltzmann Population [%]
    time    : Optional[float]           = None          # Calculation elapsed time [s]
    Erel    : float                     = np.nan        # Relative Energy [kcal/mol]
    Freq    : Optional[np.ndarray]      = None          # Vibrational Frequencies [cm-1]

    def as_dict(self) -> dict:
        """Convert record to dictionary, handling numpy arrays."""

        data = asdict(self)
        for key in ['B_vec', 'm_vec', 'Freq']:
            if data[key] is not None:
                data[key] = data[key].tolist()
        return data
    
    @classmethod
    def from_dict(cls, data: dict) -> 'EnergyRecord':
        """Create record from dictionary, restoring numpy arrays."""

        for key in ['B_vec', 'm_vec', 'Freq']:
            if data.get(key) is not None:
                data[key] = np.array(data[key])
        
        return cls(**data)


@dataclass
class EnergyStore:
    """
    Dictionary-like store for EnergyRecords indexed by protocol number.
    """
    
    data: Dict[int, EnergyRecord] = field(default_factory=dict)

    def add(self, protocol_number: int, record: EnergyRecord) -> None:
        """Add a record for a specific protocol step."""

        self.data[int(protocol_number)] = record

    def last(self) -> EnergyRecord:
        """Retrieve the record from the most recent protocol step."""
        
        if not self.data:
            return EnergyRecord()
        last_key = list(self.data.keys())[-1]
        return self.data[last_key]

    def __getitem__(self, protocol_number: int) -> EnergyRecord:
        if self.__contains__(protocol_number=protocol_number):
            return self.data.get(int(protocol_number))
        
        return EnergyRecord()

    def __contains__(self, protocol_number: int) -> bool:
        return int(protocol_number) in self.data

    def as_dict(self):
        """Used for checkpoint serialization"""
        return {k: v.as_dict() for k, v in self.data.items()}
    
    def get_energy(self) -> float: 
        data = self.last()
        if not np.isnan(data.G): 
            return data.G
        return data.E
    
    def set(self, protocol_number: int, property: str, value: Union[float, np.ndarray]):
        if not self.__contains__(protocol_number):
            raise KeyError(f"Protocol {protocol_number} not found in EnergyStore")
        
        if not hasattr(self.data[protocol_number], property):
            raise AttributeError(
                f"EnergyRecord has no attribute '{property}'. "
                f"Valid: E, G, H, S, G_E, zpve, B, B_vec, m, m_vec, Pop, time, Erel, Freq"
            )
        
        setattr(self.data[protocol_number], property, value)
    
    def log_info(self, protocol_number : int) -> Tuple[float]:
        data = self.__getitem__(int(protocol_number))
        erel = f'{data.Erel:.2f}' if not np.isnan(data.Erel) else np.nan
        pop = f'{data.Pop:.2f}' if not np.isnan(data.Pop) else np.nan

        return data.E, data.G_E, data.G, f'{data.B:.5f}', erel, pop, f'{data.time:.2f}'

    def load(self, input_dict):
        self.data = dict()
        for proto_str, vals in input_dict.get('data', {}).items():
            proto = int(proto_str)
                        
            self.data[proto] = EnergyRecord.from_dict(data=vals)

    def get_last_freq(self, protocol_number: int) -> np.ndarray: 
        
        if self.data.__getitem__(int(protocol_number)).get("Freq", None): 
            return self.data.__getitem__(int(protocol_number)).get("Freq")
    
        for i in range(protocol_number-1, -1):   
            if self.data.__getitem__(int(i)).get("Freq", None):
                return self.data.__getitem__(int(i)).get("Freq")

        return np.array([])