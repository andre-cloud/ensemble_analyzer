

from ensemble_analyzer.parsers.base import BaseParser, register_parser

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
        Parse vibrational frequencies, IR and VCD intensities.

        Returns:
            Tuple: (Frequencies array, IR spectrum array, VCD spectrum array).
        """

        if not self.regex['s_freq'] in self.fl: 
            return np.array([]), np.zeros(shape=(1,2)), np.zeros(shape=(1,2))
        
        fl = self.get_filtered_text(start=self.regex['s_freq'], end=self.regex['e_freq'])
        freq_pattern = re.compile(r'Frequencies\s*--\s*((?:[+-]?\d+\.\d+\s*)+)')
        ir_pattern   = re.compile(r'IR Inten\s*--\s*((?:[+-]?\d+\.\d+\s*)+)')
        rot_pattern  = re.compile(r'Rot\. str\.\s*--\s*((?:[+-]?\d+\.\d+\s*)+)')

        # extract and flatten into a single list of floats
        frequencies = np.array([float(x) for m in freq_pattern.findall(fl) for x in m.split()])
        ir_inten    = np.array([float(x) for m in ir_pattern.findall(fl) for x in m.split()])
        rot_str     = np.array([float(x) for m in rot_pattern.findall(fl) for x in m.split()])

        return frequencies, np.column_stack((frequencies,ir_inten)), np.column_stack((frequencies,rot_str))

    def parse_normal_modes(self, n_atoms: int) -> np.ndarray:
        if self.regex['s_freq'] not in self.fl:
            return np.empty((0, n_atoms, 3))
        fl = self.get_filtered_text(start=self.regex['s_freq'], end=self.regex['e_freq'])
        parts = fl.split('Frequencies --')
        if len(parts) < 2:
            return np.empty((0, n_atoms, 3))
        all_modes = []
        for chunk in parts[1:]:
            lines = chunk.split('\n')
            freq_vals = re.findall(r'[-+]?\d+\.\d+(?:[eE][+-]?\d+)?', lines[0])
            n_freq = min(len(freq_vals), 3)
            if n_freq == 0:
                continue
            idx = 1
            while idx < len(lines):
                line = lines[idx].strip()
                if not line:
                    idx += 1
                    continue
                if re.match(r'^\s*\d+\s+\d+', line):
                    break
                idx += 1
            mode_disps = np.zeros((n_freq, n_atoms, 3))
            for i in range(idx, min(idx + n_atoms, len(lines))):
                line = lines[i].strip()
                if not line or not re.match(r'^\d+\s+\d+', line):
                    break
                parts_line = line.split()
                if len(parts_line) < 2 + 3 * n_freq:
                    continue
                atom_idx = int(parts_line[0]) - 1
                vals = [float(x) for x in parts_line[2:]]
                for k in range(n_freq):
                    mode_disps[k, atom_idx, 0] = vals[3 * k]
                    mode_disps[k, atom_idx, 1] = vals[3 * k + 1]
                    mode_disps[k, atom_idx, 2] = vals[3 * k + 2]
            all_modes.append(mode_disps)
        if not all_modes:
            return np.empty((0, n_atoms, 3))
        return np.vstack(all_modes)

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


    



