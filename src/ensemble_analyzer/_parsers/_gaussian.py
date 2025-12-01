

from ensemble_analyzer._parsers.base import BaseParser, register_parser
from ensemble_analyzer.constants import * 

import re
import numpy as np
from typing import List, Tuple



@register_parser('gaussian')
class GaussianParser(BaseParser):

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


    def __init__(self, output_name, log):
        super().__init__(output_name, log)
        
        self.regex = self.REGEX
        self.correct_exiting = self.normal_termination()

        if not self.correct_exiting: 
            self.log.warning(self.skip_message)

    def parse_geom(self):
        fl = self.get_filtered_text(start=self.regex['geom_start'], end='--')

        pattern = r'(?:\d+)\s+(?:\d+)\s+(?:\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)'

        coords = np.array(re.findall(pattern, fl, flags=re.MULTILINE), dtype=float)
        return coords

    def parse_B_m(self):
        match_B = re.findall(self.regex['B'], self.fl)
        if match_B:
            B = np.array(match_B[-1], dtype=float)
            if self.regex['units_B'] != 'cm-1':
                B /= CONVERT_B[self.regex['units_B']]
        else:
            self.log.warning("\tB not found, storing a versor")
            B = np.array([1,0,0])

        fl = self.get_filtered_text(start='Dipole moment', end='Quadrupole')
        match_M = re.findall(self.regex['m'], fl)
        if match_M:
            M = np.array(match_M[-1], dtype=float)
        else: 
            self.log.warning("\tM not found, storing a versor")
            M = np.array([1,0,0])

        return B, M

    def parse_energy(self):
        E = re.findall(self.regex['E'], self.fl)[-1]
        return float(E)
    
    def parse_freq(self):

        if not self.regex['s_freq'] in self.fl: 
            return np.array([]), np.zeros(shape=(1,2)), np.zeros(shape=(1,2))
        
        fl = self.get_filtered_text(start=self.regex['s_freq'], end=self.regex['e_freq'])
        freq_pattern = re.compile(r'Frequencies\s*--\s*((?:[+-]?\d+\.\d+\s*)+)')
        ir_pattern   = re.compile(r'IR Inten\s*--\s*((?:[+-]?\d+\.\d+\s*)+)')
        rot_pattern  = re.compile(r'Rot\. str\.\s*--\s*((?:[+-]?\d+\.\d+\s*)+)')

        # estrai e "flatta" in un'unica lista di float
        frequencies = np.array([float(x) for m in freq_pattern.findall(fl) for x in m.split()])
        ir_inten    = np.array([float(x) for m in ir_pattern.findall(fl) for x in m.split()])
        rot_str     = np.array([float(x) for m in rot_pattern.findall(fl) for x in m.split()])

        return frequencies, np.column_stack((frequencies,ir_inten)), np.column_stack((frequencies,rot_str))

    def parse_tddft(self):

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

    def opt_done(self) -> bool:
        return len(re.findall(self.regex['opt_done'], self.fl)) >= 1

    def normal_termination(self) -> bool:
        return len(re.findall(self.regex['finish'], self.fl)) >= 1


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


    



