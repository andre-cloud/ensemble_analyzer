

from dataclasses import dataclass
from typing import Optional, Dict



@dataclass
class CalculationConfig:
    """Configuration for ensemble calculations."""
    cpu: int
    temperature: float
    start_from_protocol: int = 0
    include_H: bool = True
    
    # Graph settings
    definition: int = 4
    fwhm: Optional[Dict[str, Optional[float]]] = None
    shift: Optional[Dict[str, Optional[float]]] = None
    interested: Optional[Dict[str, Optional[float]]] = None
    invert: bool = False
    
    def __post_init__(self):
        if self.fwhm is None:
            self.fwhm = {'vibro': None, 'electro': None}
        if self.shift is None:
            self.shift = {'vibro': None, 'electro': None}
        if self.interested is None:
            self.interested = {'vibro': None, 'electro': None}
