
from dataclasses import dataclass
from typing import Optional, Union, List
import numpy as np

from ensemble_analyzer.constants import *

from ensemble_analyzer._spectral.base import BaseGraph
from ensemble_analyzer._spectral.experimental import ExperimentalGraph


@dataclass
class ComputedVibronic(BaseGraph):

    ref: Optional[ExperimentalGraph] = None

    def convolute(self, energies: np.ndarray, impulses: np.ndarray, shift: float, fwhm: float):
        # y_comp = np.zeros_like(self.X)
        # for x, y in zip (energies, impulses): 
        #     y_comp += self.lorentzian(x * shift, y, fwhm)
        return self.lorentzian(energies * shift, impulses, fwhm)