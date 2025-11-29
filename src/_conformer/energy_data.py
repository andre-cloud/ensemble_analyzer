from dataclasses import dataclass, field, asdict
from typing import Optional, Dict, Tuple, Union
import numpy as np



@dataclass
class EnergyRecord:
    E       : float                     = 0.0
    G       : float                     = np.nan
    H       : float                     = np.nan
    S       : float                     = np.nan
    G_E     : float                     = np.nan
    zpve    : float                     = np.nan
    B       : Optional[float]           = None
    B_vec   : Optional[np.ndarray]      = None
    m       : Optional[float]           = None
    m_vec   : Optional[np.ndarray]      = None
    Pop     : float                     = np.nan
    time    : Optional[float]           = None
    Erel    : float                     = np.nan
    Freq    : Optional[np.ndarray]      = None

    def as_dict(self):
        data = asdict(self)
        for key in ['B_vec', 'm_vec', 'Freq']:
            if data[key] is not None:
                data[key] = data[key].tolist()
        return data
    
    @classmethod
    def from_dict(cls, data: dict) -> 'EnergyRecord':
        for key in ['B_vec', 'm_vec', 'Freq']:
            if data.get(key) is not None:
                data[key] = np.array(data[key])
        
        return cls(**data)

    

@dataclass
class EnergyStore:
    data: Dict[int, EnergyRecord] = field(default_factory=dict)

    def add(self, protocol_number: int, record: EnergyRecord):
        self.data[int(protocol_number)] = record

    def last(self) -> EnergyRecord:
        if not self.data:
            return EnergyRecord()
        last_key = list(self.data.keys())[-1]
        return self.data[last_key]

    def __getitem__(self, protocol_number: int) -> EnergyRecord:
        if self.__contains__(protocol_number=protocol_number):
            return self.data.get(int(protocol_number))
        
        return None

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

        return data.E, data.G_E, data.G, data.B, erel, pop, f'{data.time:.2f}'

    def load(self, input_dict):
        self.data = dict()
        for proto_str, vals in input_dict.get('data', {}).items():
            proto = int(proto_str)
                        
            self.data[proto] = EnergyRecord.from_dict(data=vals)