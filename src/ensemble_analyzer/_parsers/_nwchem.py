from ensemble_analyzer._parsers.base import BaseParser, register_parser
from ensemble_analyzer.constants import *

import re
import numpy as np
from typing import List, Tuple


@register_parser('nwchem')
class NWChemParser(BaseParser):
    """
    Parser implementation for NWChem output files.
    """

    REGEX = {
        "B": r"Rotational Constants[\s\S]*?A=\s*(-?\d+\.\d+)\s+cm-1\s+\([^)]*\)\s*[\s\S]*?B=\s*(-?\d+\.\d+)\s+cm-1\s+\([^)]*\)\s*[\s\S]*?C=\s*(-?\d+\.\d+)\s+cm-1",
        "units_B": "cm-1",
        "m": r"Nuclear Dipole moment[\s\S]*?(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)",
        "E": r"(?:Total DFT energy|SCF energy)\s*=\s*(-?\d+\.\d+)",
        "break": "\n\n",
        "idx_en_tddft": 1,
        "idx_imp_tddft": 2,
        "idx_en_ir": 1,
        "idx_imp_ir": 2,
        "idx_imp_vcd": None,
        "s_freq": "NORMAL MODE EIGENVECTORS IN CARTESIAN COORDINATES",
        "s_IR": "Projected Infra Red Intensities",
        "s_UV": "Excitation energies",
        "s_ECD": "CD Spectrum",
        "geom_start": "Output coordinates in angstroms",
        "finish": "Total times",
        "opt_done": "Optimization converged",
        "ext": "log",
    }

    def __init__(self, output_name: str, log, conf=None) -> None:
        super().__init__(output_name, log, conf)
        self.regex = self.REGEX
        self.correct_exiting = self.normal_termination()
        if not self.correct_exiting:
            self.log.warning(self.skip_message)

    def parse_geom(self) -> np.ndarray:
        fl = self.get_filtered_text(start=self.regex["geom_start"], end="Atomic Mass")

        pattern = r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
        coords = np.array(re.findall(pattern, fl, flags=re.MULTILINE), dtype=float)
        return coords

    def parse_energy(self) -> float:
        match = re.findall(self.regex["E"], self.fl)
        if not match:
            return 0.0
        return float(match[-1])

    def parse_B_m(self) -> Tuple[np.ndarray, np.ndarray]:
        match_B = re.findall(self.regex["B"], self.fl)
        if match_B:
            B = np.array(match_B[-1], dtype=float)
            if self.regex["units_B"] != "cm-1":
                B /= CONVERT_B[self.regex["units_B"]]
        else:
            self.log.warning(f"\t{self.log.WARNING} B not found, calculating with ASE")
            B = self.calculate_B()

        match_M = re.findall(self.regex["m"], self.fl)
        if match_M:
            M = np.array(match_M[-1], dtype=float)
        else:
            self.log.warning(f"\t{self.log.WARNING} M not found, storing a versor")
            M = np.array([1, 0, 0])

        return B, M

    def parse_freq(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        if self.regex["s_IR"] not in self.fl:
            return np.array([]), np.zeros(shape=(1, 2)), np.zeros(shape=(1, 2))

        ir_fl = self.get_filtered_text(start=self.regex["s_IR"], end="\n\n")
        pattern = r"\s+(\d+)\s+([-\d.]+)\s+\|\|\s+[-\d.]+\s+[-\d.]+\s+([-\d.]+)\s+[-\d.]+"
        matches = re.findall(pattern, ir_fl)

        if matches:
            freq = np.array([float(m[1]) for m in matches])
            ir_arr = np.array([float(m[2]) for m in matches])
            mask = np.abs(freq) > 1.0
            freq = freq[mask]
            ir_arr = ir_arr[mask]
            ir = np.column_stack((freq, ir_arr)) if len(freq) > 0 else np.zeros(shape=(1, 2))
        else:
            freq = np.array([])
            ir = np.zeros(shape=(1, 2))

        vcd = np.zeros(shape=(1, 2))
        return freq, ir, vcd

    def parse_tddft(self) -> Tuple[np.ndarray, np.ndarray]:
        uv = np.zeros(shape=(1, 2))
        ecd = np.zeros(shape=(1, 2))

        if self.regex["s_UV"] in self.fl:
            uv_text = self.get_filtered_text(start=self.regex["s_UV"], end="\n\n").splitlines()
            uv_data = self.parse_table(uv_text, [self.regex["idx_en_tddft"], self.regex["idx_imp_tddft"]])
            if uv_data:
                uv = np.array(uv_data, dtype=np.float64)

        if self.regex["s_ECD"] in self.fl:
            ecd_text = self.get_filtered_text(start=self.regex["s_ECD"], end="\n\n").splitlines()
            ecd_data = self.parse_table(ecd_text, [self.regex["idx_en_tddft"], self.regex["idx_imp_tddft"]])
            if ecd_data:
                ecd = np.array(ecd_data, dtype=np.float64)

        return uv, ecd

    def opt_done(self) -> bool:
        return len(re.findall(self.regex["opt_done"], self.fl)) >= 1

    def normal_termination(self) -> bool:
        return len(re.findall(self.regex["finish"], self.fl)) >= 1




if __name__ == '__main__': 
    from mock import MagicMock
    p = NWChemParser('water_nwchem.nwo', log=MagicMock())

    print(p.parse_geom())
    print(p.parse_B_m())
    print(p.opt_done())