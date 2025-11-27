from dataclasses import dataclass, field
from typing import Optional, Dict, Union
import numpy as np



@dataclass
class EnergyRecord:
    E       : float                     = 0.0
    G       : float                     = np.nan
    H       : float                     = np.nan
    G_E     : float                     = np.nan
    zpve    : float                     = np.nan
    B       : Optional[float]           = None
    B_vec   : Optional[np.ndarray]      = None
    m       : Optional[float]           = None
    m_vec   : Optional[np.ndarray]      = None
    Pop     : Optional[float]           = None
    time    : Optional[float]           = None
    Erel    : Optional[float]           = None
    Freq    : Optional[np.ndarray]      = None

    def as_dict(self):
        return {
            "E"     : self.E,
            "G"     : self.G,
            "H"     : self.H,
            "zpve"  : self.zpve,
            "B"     : self.B,
            "B_vec" : self.B_vec,
            "m"     : self.m,
            "m_vec" : self.m_vec,
            "Pop"   : self.Pop,
            "time"  : self.time,
            "Erel"  : self.Erel,
            "G-E"   : self.G_E if not np.isnan(self.G) else np.nan,
            "Freq"  : self.Freq,
        }
    

@dataclass
class EnergyStore:
    data: Dict[str, EnergyRecord] = field(default_factory=dict)

    def add(self, protocol_number: int, record: EnergyRecord):
        self.data[int(protocol_number)] = record

    def last(self) -> EnergyRecord:
        if not self.data:
            return EnergyRecord()
        last_key = list(self.data.keys())[-1]
        return self.data[last_key]

    def __getitem__(self, protocol_number: int) -> EnergyRecord:
        return self.data.get(int(protocol_number), EnergyRecord())

    def __contains__(self, protocol_number: int) -> bool:
        return int(protocol_number) in self.data

    def as_dict(self):
        """Used for checkpoint serialization"""
        return {k: v.as_dict() for k, v in self.data.items()}
    
    def get_energy(self) -> float: 
        data = self.data.last()
        if ~np.isnan(data.G) or data.G is not None: 
            return data.G
        return data.E
