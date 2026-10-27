

from ensemble_analyzer.parsers.base import BaseParser, register_parser
from ensemble_analyzer.parsers.vibrations import (
    gaussian_get_freq, gaussian_get_ir_intensities, gaussian_get_normal_modes,
    gaussian_align_normal_modes, to_displacements,
)

import re
import numpy as np
from typing import List, Tuple



@register_parser('gaussian')
class GaussianParser(BaseParser):
    """
    Parser implementation for Gaussian output files.
    """

    REGEX = {
        "B": r"Rotational constants \(GHZ\):\s*(-?\d+.\d*)\s*(-?\d+.\d*)\s*(-?\d+.\d*)",
        'units_B': 'GHz',
        "m": r"X=\s+(-?\d+.\d+)\s+Y=\s+(-?\d+.\d+)\s+Z=\s+(-?\d+.\d+)",
        "E": r"SCF Done.* =\s+(-?\d+.\d+)\s*",
        "break": "\n\n",
        "idx_en_tddft": None,  # index for energy in the UV & ECD table in eV
        "idx_imp_tddft": None,  # index for oscillator strength in the UV table
        "idx_en_ir": None,  # index for energy in the IR table in cm**-1
        "idx_imp_ir": None,  # index for oscillator strength in the IR table
        "idx_imp_vcd": None,  # index for oscillator strength in the VCD table
        "s_freq": "Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering",
        "start_spec": "Excited states from",
        "e_freq": "\n\n\n",
        "idx_freq": 1,  # index for frequency in frequency table
        # frequency section (same keys as frogchem regex_gaussian)
        "freq_section": "Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering",
        "freq_end": "\n\n\n",
        "freq_pattern": r'Frequencies\s*--\s*((?:[+-]?\d+\.\d+\s*)+)',
        "ir_pattern": r'IR Inten\s*--\s*((?:[+-]?\d+\.\d+\s*)+)',
        "rot_pattern": r'Rot\. str\.\s*--\s*((?:[+-]?\d+\.\d+\s*)+)',
        "opt_done": "Optimization completed",
        "geom_start": """Input orientation:                          
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------""",
        "finish": "Normal termination",
        "ext": "log"
    }


    def parse_geom(self) -> np.ndarray:
        """
        Parse final geometry coordinates.

        Returns:
            np.ndarray: Cartesian coordinates array.
        """
        fl = self.get_filtered_text(start=self.regex['geom_start'], end='--')

        pattern = r'(?:\d+)\s+(?:\d+)\s+(?:\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)'

        coords = np.array(re.findall(pattern, fl, flags=re.MULTILINE), dtype=float)
        return coords

    def _get_dipole_text(self) -> str:
        return self.get_filtered_text(start='Dipole moment', end='Quadrupole')
    
    def parse_freq(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Parse vibrational frequencies, IR and VCD intensities (FROGChem parser).

        Returns:
            Tuple: (Frequencies array, IR spectrum array, VCD spectrum array).
        """

        if not self.regex['freq_section'] in self.fl:
            return np.array([]), np.zeros(shape=(1,2)), np.zeros(shape=(1,2))

        frequencies = gaussian_get_freq(self.fl, self.regex)
        ir = gaussian_get_ir_intensities(self.fl, self.regex)

        # VCD: rotational strengths, only printed with freq=VCD
        fl = self.get_filtered_text(start=self.regex['freq_section'], end=self.regex['freq_end'])
        rot_str = np.array([float(x) for m in re.findall(self.regex['rot_pattern'], fl) for x in m.split()])
        if len(rot_str) == len(frequencies) and len(rot_str) > 0:
            vcd = np.column_stack((frequencies, rot_str))
        else:
            vcd = np.zeros(shape=(1, 2))

        return frequencies, ir, vcd

    def parse_normal_modes(self, n_atoms: int) -> np.ndarray:
        """Cartesian displacement vectors of the vibrational modes (FROGChem parser).

        Gaussian prints them in the standard orientation; they are rotated onto the input
        orientation, i.e. the frame of parse_geom.

        Returns:
            np.ndarray: (n_modes, n_atoms, 3), same order as the frequencies of parse_freq.
            Empty (0, n_atoms, 3) array when the frequency section is missing.
        """
        empty = np.empty((0, n_atoms, 3))
        try:
            matrix = gaussian_get_normal_modes(self.fl, self.regex, N=3 * n_atoms)
        except ValueError:
            return empty
        except Exception as err:
            self.log.warning(f"Normal modes could not be parsed: {err}")
            return empty
        if matrix.shape[1] == 0:
            return empty
        matrix = gaussian_align_normal_modes(self.fl, matrix)
        return to_displacements(matrix, n_atoms)

    def parse_tddft(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Parse TD-DFT excited states for UV and ECD spectra.

        Returns:
            Tuple: (UV data array, ECD data array).
        """


        if not self.regex['start_spec'] in self.fl: 
            return np.zeros(shape=(1,2)), np.zeros(shape=(1,2))

        uv = self.get_filtered_text(start='Excitation energies and oscillator strengths', end='***')
        uv_pattern = re.compile(r'Excited State\s+\d+:.*(\d+\.\d+ )eV.*f=(\d+.\d+)')
        impulse = np.array([m for m in uv_pattern.findall(uv)], dtype=np.float64)
        energies, f = impulse[:,0], impulse[:,1]

        ecd = self.get_filtered_text(start='<0|del|b> * <b|rxdel|0> + <0|del|b> * <b|delr+rdel|0>', end='1/2[<0|r|b>*<b|rxdel|0> + (<0|rxdel|b>*<b|r|0>)*]')
        ecd_pattern = re.compile(r'^\s*\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+([+-]?\d+\.\d+)', re.MULTILINE)
        R = np.array([m for m in ecd_pattern.findall(ecd)], dtype=np.float64)

        return np.column_stack((energies,f)), np.column_stack((energies,R))

if __name__ == '__main__':

    import mock
    parser = GaussianParser("files/gaussian.log", mock.MagicMock())
    B,M = parser.parse_B_m()
    print(f'{B,M=}')
    geom = parser.parse_geom()
    print(f'{geom=}')
    energy = parser.parse_energy()
    print(f'{energy=}')
    print(f'{parser.opt_done()=}')
    print(f'{parser.normal_termination()=}')
    frequencies, ir_inten, rot_str = parser.parse_freq()
    print(f'{frequencies, ir_inten, rot_str=}')
    parser = GaussianParser("files/gaussian_tddft.log", mock.MagicMock())
    uv, ecd = parser.parse_tddft()
    print(f'{uv, ecd=}')


    



