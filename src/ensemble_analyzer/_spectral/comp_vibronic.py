
from dataclasses import dataclass
from typing import Optional, Union, List
import numpy as np

from ensemble_analyzer.constants import *

from ensemble_analyzer._spectral.base import BaseGraph
from ensemble_analyzer._spectral.experimental import ExperimentalGraph


@dataclass
class ComputedVibronic(BaseGraph):

    ref: Optional[ExperimentalGraph] = None

    def convolute(self, energies: np.ndarray, impulses: np.ndarray, shift: float, fwhm: float) -> np.ndarray:
        """
        Apply Lorentzian convolution to vibronic transitions.

        Used for IR and VCD spectra. The shift is applied as a scaling factor
        (multiplicative) to the frequencies.

        Args:
            energies (np.ndarray): Array of transition frequencies [cm^-1].
            impulses (np.ndarray): Array of transition intensities.
            shift (float): Scaling factor for frequency correction (x * shift).
            fwhm (float): Full Width at Half Maximum for Lorentzian broadening.

        Returns:
            np.ndarray: Convolved spectral intensity array.
        """

        return self.lorentzian(energies * shift, impulses, fwhm)