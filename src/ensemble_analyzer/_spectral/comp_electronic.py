
from dataclasses import dataclass
from typing import Optional, Union, List
import numpy as np


from ensemble_analyzer.constants import *

from ensemble_analyzer._spectral.base import BaseGraph
from ensemble_analyzer._spectral.experimental import ExperimentalGraph


@dataclass
class ComputedElectronic(BaseGraph):

    ref: Optional[ExperimentalGraph] = None

    def convolute(self, energies: np.ndarray, impulses: np.ndarray, shift: float, fwhm: float):
        
        # POSITIVE SHIFT = BLUE SHIFT
        return self.gaussian(energies + shift, impulses, fwhm)


