from dataclasses import dataclass
from typing import Literal



@dataclass
class GraphDefault: 
    graph_type: Literal['IR', 'VCD', 'UV', 'ECD']
    
    def __post_init__(self):
        # units of the graph
        self.units = {
            'IR':   'cm**-1', 
            'VCD':  'cm**-1', 
            'UV':   'eV',
            'ECD':  'eV',
        }[self.graph_type]

        # lower boundaries of the spectra 
        self.start = {
            'IR':   300, 
            'VCD':  300, 
            'UV':   1.25, 
            'ECD':  1.25, 
        }[self.graph_type]

        # upper boundaries of the spectra 
        self.end = {
            'IR':   7000, 
            'VCD':  7000, 
            'UV':   9, 
            'ECD':  9,
        }[self.graph_type]

        # Shift of each peak
        self.shift = {
            'IR':   1, 
            'VCD':  1, 
            'UV':   0, 
            'ECD':  0,
        }[self.graph_type]

        # Full Width Half Maximum of each peak
        self.fwhm = {
            'IR':   10, 
            'VCD':  10, 
            'UV':   0.25, 
            'ECD':  0.25,
        }[self.graph_type]

        # Default interval for the shift
        self.shift_intervals = {
            'IR':   [0.85, 1.05], 
            'VCD':  [0.85, 1.05], 
            'UV':   [-0.5, 0.5], 
            'ECD':  [-0.5, 0.5],
        }[self.graph_type]
        
        # Default interval for the FWHM
        self.fwhm_intervals = {
            'IR':   [   4,   20], 
            'VCD':  [   4,   20], 
            'UV':   [0.20, 0.50], 
            'ECD':  [0.20, 0.50],
        }[self.graph_type]

        self.experimental_fname = f'{self.graph_type.lower()}_ref.dat'

        self.axis_label = {
            "IR": {
                "x": r"Wavenumber $\widetilde \nu$ [$cm^{-1}$]",
                "y": r"Intensity [$a.u.$]",
            },
            "VCD": {
                "x": r"Wavenumber $\widetilde \nu$ [$cm^{-1}$]",
                "y": r"Intensity [$a.u.$]",
            },
            "UV": {"x": r"Energy [$eV$]", "y": r"Intensity [$a.u.$]"},
            "ECD": {"x": r"Energy [$eV$]", "y": r"Intensity [$a.u.$]"},
        }[self.graph_type]

        self.X_buffer = {
            "IR":   100,
            "VCD":  100,
            "UV":   0.2,
            "ECD":  0.2,
        }[self.graph_type]

        self.interested_area = {
            "IR":   [1100, 1400],
            "VCD":  [1100, 1400],
            "UV":   [4.96, 6.2],
            "ECD":  [4.96, 6.2],
        }[self.graph_type]

