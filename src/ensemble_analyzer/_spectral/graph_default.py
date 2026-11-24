from dataclasses import dataclass
from typing import Literal

from ensemble_analyzer.constants import VIBRO_OR_ELECTRO


# Defaults per spectrum family: IR/VCD share 'vibro' [cm^-1], UV/ECD share 'electro' [eV]
_DEFAULTS = {
    'vibro': dict(
        start=300, end=7000,              # spectral range
        shift=1, fwhm=10,                 # shift (scaling factor) and FWHM of each peak
        shift_intervals=[0.85, 1.05], fwhm_intervals=[4, 20],
        axis_label={"x": r"Wavenumber $\widetilde \nu$ [$cm^{-1}$]", "y": r"Intensity [$a.u.$]"},
        X_buffer=100, interested_area=[1100, 1400],
    ),
    'electro': dict(
        start=1.25, end=9,
        shift=0, fwhm=0.25,               # additive shift [eV]
        shift_intervals=[-0.5, 0.5], fwhm_intervals=[0.20, 0.50],
        axis_label={"x": r"Energy [$eV$]", "y": r"Intensity [$a.u.$]"},
        X_buffer=0.2, interested_area=[4.96, 6.2],
    ),
}


@dataclass
class GraphDefault:
    """Default spectral parameters (range, shift, FWHM, labels) for each graph type."""

    graph_type: Literal['IR', 'VCD', 'UV', 'ECD']

    def __post_init__(self) -> None:
        for key, value in _DEFAULTS[VIBRO_OR_ELECTRO[self.graph_type]].items():
            setattr(self, key, list(value) if isinstance(value, list) else value)
        self.experimental_fname = f'{self.graph_type.lower()}_ref.dat'
