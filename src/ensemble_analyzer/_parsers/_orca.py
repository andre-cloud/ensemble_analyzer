
from ensemble_analyzer._parsers.base import BaseParser, register_parser
from ensemble_analyzer.constants import * 


import re
import numpy as np
from typing import List, Tuple



@register_parser('orca')
class OrcaParser(BaseParser):

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
    
    def __init__(self, output_name, log):
        super().__init__(output_name, log)
        self.version = self.get_version()
        self.regex = self.REGEX[self.version]

        self.correct_exiting = self.normal_termination()

        if not self.correct_exiting: 
            self.log.warning(self.skip_message)


    def get_version(self):
        find = re.findall(r'Program Version (\d)', self.fl)
        return find[0]

    def parse_geom(self):
        fl = self.get_filtered_text(start = self.regex['geom_start'], end = self.regex['break'])

        pattern = r'\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'

        coords = np.array(re.findall(pattern, fl, flags=re.MULTILINE), dtype=float)
        return coords
    
    def parse_energy(self):
        E = re.findall(self.regex['E'], self.fl)[-1]
        return float(E)

    def parse_B_m(self) -> Tuple:
        match_B = re.findall(self.regex['B'], self.fl)
        if match_B:
            B = np.array(match_B[-1], dtype=float)
            if self.regex['units_B'] != 'cm-1':
                B /= CONVERT_B[self.regex['units_B']]
        else:
            self.log.warning("\tB not found, storing a versor")
            B = np.array([1,0,0])

        match_M = re.findall(self.regex['m'], self.fl)
        if match_M:
            M = np.array(match_M[-1], dtype=float)
        else: 
            self.log.warning("\tM not found, storing a versor")
            M = np.array([1,0,0])

        return B, M

    def parse_freq(self) -> List:
        
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

    def parse_tddft(self):

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

    def opt_done(self) -> bool:
        return len(re.findall(self.regex['opt_done'], self.fl)) == 1
    
    def normal_termination(self) -> bool:
        return len(re.findall(self.regex['finish'], self.fl)) == 1

if __name__ == '__main__':

    import mock
    print('ORCA 6')
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
    parser = OrcaParser("files/tddft_6.out", mock.MagicMock())
    uv, ecd = parser.parse_tddft()
    print(uv, ecd)
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