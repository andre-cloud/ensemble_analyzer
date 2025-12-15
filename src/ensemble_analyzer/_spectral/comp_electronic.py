
from dataclasses import dataclass
from typing import Optional, Union, List
import numpy as np


from ensemble_analyzer.constants import *

from ensemble_analyzer._spectral.base import BaseGraph
from ensemble_analyzer._spectral.experimental import ExperimentalGraph


@dataclass
class ComputedElectronic(BaseGraph):

    ref: Optional[ExperimentalGraph] = None

    def convolute(self, energies: np.ndarray, impulses: np.ndarray, shift: float, fwhm: float) -> np.ndarray:
        """
        Apply Gaussian convolution to electronic transitions.

        Used for UV and ECD spectra. The shift is applied as an additive term
        (rigid shift) to the excitation energies.

        Args:
            energies (np.ndarray): Array of excitation energies [eV].
            impulses (np.ndarray): Array of oscillator/rotational strengths.
            shift (float): Additive shift value [eV] (x + shift).
            fwhm (float): Full Width at Half Maximum for Gaussian broadening.

        Returns:
            np.ndarray: Convolved spectral intensity array.
        """
        
        # POSITIVE SHIFT = BLUE SHIFT
        return self.gaussian(energies + shift, impulses, fwhm)


