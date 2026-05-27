
from ensemble_analyzer.parsers.base import BaseParser, register_parser
from ensemble_analyzer.constants import FACTOR_EV_CM_1

import re
import numpy as np
from typing import List, Tuple



@register_parser('orca')
class OrcaParser(BaseParser):
    """
    Parser implementation for ORCA output files.
    Supports ORCA 5 and ORCA 6 syntax variations.
    """

    REGEX = {
        '6':{
        "B": r"Rotational constants in cm-1:\s*(-?\d+.\d*)\s*(-?\d+.\d*)\s*(-?\d+.\d*)",
        'units_B': 'cm-1',
        "m": r"Total Dipole Moment\s*:\s*([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)",
        "E": r"FINAL SINGLE POINT ENERGY.*?(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
        "start_spec": "SPECTRA",
        "end_spec": "***",
        "s_UV": """ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS    
----------------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength fosc(D2)      D2        DX        DY        DZ   
                      (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)  
----------------------------------------------------------------------------------------------------""",
        "s_ECD": """CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS    
------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength    R        MX        MY        MZ   
                      (eV)      (cm-1)    (nm)   (1e40*cgs)   (au)      (au)      (au)  
------------------------------------------------------------------------------------------""",
        "s_IR": """Mode   freq       eps      Int      T**2         TX        TY        TZ
       cm**-1   L/(mol*cm) km/mol    a.u.
----------------------------------------------------------------------------
""",
        "s_VCD": """Mode   Freq    VCD-Intensity    
       (1/cm) (1E-44*esu^2*cm^2) 
---------------------------------""",
        "break": "\n\n",
        "idx_en_tddft": 3,  # index for energy in the UV & ECD table in eV
        "idx_imp_tddft": 6,  # index for oscillator strength in the UV table
        "idx_en_ir": 1,  # index for energy in the IR table in cm**-1
        "idx_imp_ir": 3,  # index for oscillator strength in the IR table
        "idx_en_vcd": 1,  # index for energy in the VCD table in cm**-1
        "idx_imp_vcd": 2,  # index for oscillator strength in the VCD table
        "s_freq": "VIBRATIONAL FREQUENCIES",
        "e_freq": "------------",
        "idx_freq": 1,  # index for frequency in frequency table
        "opt_done": "THE OPTIMIZATION HAS CONVERGED",
        "geom_start": """CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------""",
        "finish": "ORCA TERMINATED NORMALLY",
        "ext": "out",
        },

        "5": {
        "B": r"Rotational constants in cm-1:\s*(-?\d+.\d*)\s*(-?\d+.\d*)\s*(-?\d+.\d*)",
        'units_B': 'cm-1',
        "m": r"Total Dipole Moment\s*:\s*([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)",
        "E": r"FINAL SINGLE POINT ENERGY\s*(-?\d*.\d*)",
        "start_spec": "SPECTRA",
        "end_spec": "***",
        "s_UV": """ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
-----------------------------------------------------------------------------
State   Energy    Wavelength  fosc         T2        TX        TY        TZ  
        (cm-1)      (nm)                 (au**2)    (au)      (au)      (au) 
-----------------------------------------------------------------------------""",
        "s_ECD": """CD SPECTRUM
-------------------------------------------------------------------
State  Energy     Wavelength     R         MX        MY        MZ   
       (cm-1)       (nm)     (1e40*cgs)   (au)      (au)      (au)  
-------------------------------------------------------------------""",
        "s_IR": """Mode   freq       eps      Int      T**2         TX        TY        TZ
       cm**-1   L/(mol*cm) km/mol    a.u.
----------------------------------------------------------------------------
""",
        "s_VCD": None,
        "break": "\n\n",
        "idx_en_tddft": 1,  # index for energy in the UV & ECD table in eV
        "idx_imp_tddft": 3,  # index for oscillator strength in the UV table
        "idx_en_ir": 1,  # index for energy in the IR table in cm**-1
        "idx_imp_ir": 3,  # index for oscillator strength in the IR table
        "idx_imp_vcd": 2,  # index for oscillator strength in the VCD table
        "s_freq": "VIBRATIONAL FREQUENCIES",
        "e_freq": "------------",
        "idx_freq": 1,  # index for frequency in frequency table
        "opt_done": "THE OPTIMIZATION HAS CONVERGED",
        "geom_start": """CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------""",
        "finish": "ORCA TERMINATED NORMALLY",
        "ext": "out",}
    }
    
    def _init_regex(self) -> None:
        self.version = self.get_version()
        self.regex = self.REGEX[self.version]

    def get_version(self) -> str:
        """
        Detect ORCA major version from output header.

        Returns:
            str: Version string (e.g., '5' or '6').
        """
        find = re.findall(r'Program Version (\d)', self.fl)
        return find[0] if find else "0"

    def parse_geom(self) -> np.ndarray:
        """
        Parse final geometry coordinates.

        Returns:
            np.ndarray: Cartesian coordinates array.
        """
        fl = self.get_filtered_text(start = self.regex['geom_start'], end = self.regex['break'])

        pattern = r'\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'

        coords = np.array(re.findall(pattern, fl, flags=re.MULTILINE), dtype=float)
        return coords
    


    def parse_freq(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Parse vibrational frequencies, IR and VCD intensities.

        Returns:
            Tuple: (Frequencies array, IR spectrum array, VCD spectrum array).
        """
        
        if not self.regex['s_freq'] in self.fl: 
            return np.array([]), np.zeros(shape=(1,2)), np.zeros(shape=(1,2))

        fl = self.get_filtered_text(self.regex['s_freq'], end='\n\n\n')
    
        pattern = r'(?:\d+:)\s*(-?\d+.\d*)'
        # freq
        freq = np.array(re.findall(pattern, fl, flags=re.MULTILINE), dtype=float)
        freq = freq[freq!=0]

        # IR
        ir_text = self.get_filtered_text(start=self.regex['s_IR'], end='\n\n').splitlines()
        ir = np.array(self.parse_table(ir_text, [self.regex['idx_en_ir'], self.regex['idx_imp_ir']]), dtype=np.float64)

        # VCD
        if self.version == '5' or self.regex['s_VCD'] not in self.fl:
            vcd = np.zeros(shape=(1,2))
        else: 
            vcd_text = self.get_filtered_text(start=self.regex['s_VCD'], end='\n\n').splitlines()
            vcd = np.array(self.parse_table(vcd_text, [self.regex['idx_en_vcd'], self.regex['idx_imp_vcd']]), dtype=np.float64)

        return freq, ir, vcd

    def parse_normal_modes(self, n_atoms: int) -> np.ndarray:
        if "NORMAL MODES" not in self.fl:
            return np.empty((0, n_atoms, 3))
        try:
            block = self.fl.split("NORMAL MODES")[-1].split("\n\n\n")[0]
            lines = block.split("\n")
            n_total = 3 * n_atoms
            max_col = -1
            for line in lines:
                stripped = line.strip()
                if not stripped or stripped.startswith("-"):
                    continue
                header_match = re.match(r'^\s*(\d+(?:\s+\d+)*)\s*$', stripped)
                if header_match:
                    cols = [int(x) for x in header_match.group(1).split()]
                    if cols and cols[-1] > max_col:
                        max_col = cols[-1]
            if max_col < 0:
                return np.empty((0, n_atoms, 3))
            matrix = np.zeros((n_total, max_col + 1))
            current_columns = []
            for line in lines:
                stripped = line.strip()
                if not stripped or stripped.startswith("-"):
                    continue
                header_match = re.match(r'^\s*(\d+(?:\s+\d+)*)\s*$', stripped)
                if header_match:
                    current_columns = [int(x) for x in header_match.group(1).split()]
                    continue
                data_match = re.match(
                    r'^\s*(\d+)\s+((?:[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?\s*)+)', stripped
                )
                if data_match and current_columns:
                    row_idx = int(data_match.group(1))
                    coeff_str = data_match.group(2)
                    coeffs = [float(c) for c in re.findall(
                        r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', coeff_str
                    )]
                    for i, c in enumerate(coeffs):
                        if i < len(current_columns) and row_idx < n_total:
                            col = current_columns[i]
                            matrix[row_idx, col] = c
            nonzero = [c for c in range(matrix.shape[1])
                       if np.sum(matrix[:n_total, c] ** 2) > 1e-14]
            if not nonzero:
                return np.empty((0, n_atoms, 3))
            all_modes = np.zeros((len(nonzero), n_atoms, 3))
            for i, col in enumerate(nonzero):
                vec = matrix[:n_total, col]
                all_modes[i] = vec.reshape(n_atoms, 3)
            return all_modes
        except Exception:
            return np.empty((0, n_atoms, 3))

    def parse_tddft(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Parse TD-DFT excited states for UV and ECD spectra.

        Returns:
            Tuple: (UV data array, ECD data array).
        """

        if not self.regex['start_spec'] in self.fl: 
            return np.zeros(shape=(1,2)), np.zeros(shape=(1,2))

        # UV
        uv_text = self.get_filtered_text(start=self.regex['s_UV'], end=self.regex['break']).splitlines()
        uv = np.array(self.parse_table(uv_text, [self.regex['idx_en_tddft'], self.regex['idx_imp_tddft']]), dtype=np.float64)
        if self.version=='5':
            uv[:, 0] = FACTOR_EV_CM_1/uv[:, 0]

        # ECD
        ecd_text = self.get_filtered_text(start=self.regex['s_ECD'], end=self.regex['break']).splitlines()
        ecd = np.array(self.parse_table(ecd_text, [self.regex['idx_en_tddft'], self.regex['idx_imp_tddft']]), dtype=np.float64)
        if self.version=='5':
            ecd[:, 0] = FACTOR_EV_CM_1/ecd[:, 0]

        return uv, ecd

if __name__ == '__main__':

    import mock
    # print('ORCA 6')
    # parser = OrcaParser("files/opt_6.out", mock.MagicMock())
    # B,M = parser.parse_B_m()
    # print(B,M)
    # geom = parser.parse_geom()
    # print(parser.opt_done())
    # print(geom)
    # E = parser.parse_energy()
    # print(E)
    # freq, ir, vcd = parser.parse_freq()
    # print(freq, ir, vcd)
    # parser = OrcaParser("files/tddft_6.out", mock.MagicMock())
    # uv, ecd = parser.parse_tddft()
    # print(uv, ecd)
    # print('='*10)
    # print('ORCA 5')
    # parser = OrcaParser("files/opt_5.out", mock.MagicMock())
    # geom = parser.parse_geom()
    # print(parser.opt_done())
    # print(geom)
    # E = parser.parse_energy()
    # print(E)
    # freq, ir, vcd = parser.parse_freq()
    # print(freq, ir, vcd)
    # parser = OrcaParser("files/tddft_5.out", mock.MagicMock())
    # uv, ecd = parser.parse_tddft()
    # print(uv, ecd)