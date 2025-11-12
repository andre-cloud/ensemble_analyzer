
from dataclasses import dataclass
from typing import Optional, Union, List
import numpy as np


from src.constants import *

from src._spectral.base import BaseGraph
from src._spectral.experimental import ExperimentalGraph


@dataclass
class ComputedElectronic(BaseGraph):

    ref: Optional[ExperimentalGraph] = None

    def convolute(self, energies: np.ndarray, impulses: np.ndarray, shift: float, fwhm: float):
        
        # POSITIVE SHIFT = BLUE SHIFT
        return self.gaussian(energies + shift, impulses, fwhm)


