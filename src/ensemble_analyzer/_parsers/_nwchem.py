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
        "s_freq": "NORMAL MODE EIGENVECTORS IN CARTESIAN COORDINATES",
        "s_IR": "Projected Infra Red Intensities",
        "geom_start": "Output coordinates in angstroms",
        "finish": "Total times",
        "opt_done": "Optimization converged",
        "ext": "log",
    }

    def __init__(self, output_name: str, log, conf=None) -> None:
        """Initialize NWChem parser.

        Args:
            output_name: Path to NWChem output file.
            log: Logger instance.
            conf: Optional Conformer object.
        """
        super().__init__(output_name, log, conf)
        self.regex = self.REGEX
        self.correct_exiting = self.normal_termination()
        if not self.correct_exiting:
            self.log.warning(self.skip_message)

    def parse_geom(self) -> np.ndarray:
        """Parse the final geometry from NWChem output.

        Returns:
            np.ndarray: Cartesian coordinates array of shape (N, 3).
        """
        fl = self.get_filtered_text(start=self.regex["geom_start"], end="Atomic Mass")

        pattern = r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
        coords = np.array(re.findall(pattern, fl, flags=re.MULTILINE), dtype=float)
        return coords

    def parse_energy(self) -> float:
        """Parse the final electronic energy from NWChem output.

        Returns:
            float: Energy in Hartree.
        """
        match = re.findall(self.regex["E"], self.fl)
        if not match:
            return 0.0
        return float(match[-1])

    def parse_B_m(self) -> Tuple[np.ndarray, np.ndarray]:
        """Parse rotational constants and dipole moment from NWChem output.

        Returns:
            Tuple[np.ndarray, np.ndarray]: (B vector, dipole moment vector).
        """
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
        """Parse vibrational frequencies and IR intensities from NWChem output.

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]:
                (frequencies, IR spectrum data, VCD spectrum data).
        """
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
        """Parse TD-DFT excited states (UV/ECD) from NWChem output.

        Returns:
            Tuple[np.ndarray, np.ndarray]: (UV data array, ECD data array).
        """
        uv = np.zeros(shape=(1, 2))
        ecd = np.zeros(shape=(1, 2))

        root_pat = r"(\d+)\s+singlet a\s+[-\d.]+\s+a\.u\.\s+([-\d.]+)\s+eV"
        os_pat = r"Dipole Oscillator Strength\s+([-\d.]+)"
        rs_pat = r"Rotatory Strength \(1E-40 esu\*\*2cm\*\*2\):\s+([-\d.]+)"

        uv_list = []
        ecd_list = []
        for chunk in self.fl.split("Root "):
            m = re.search(root_pat, chunk)
            if not m:
                continue
            ev = float(m.group(2))
            os_m = re.search(os_pat, chunk)
            rs_m = re.search(rs_pat, chunk)
            if os_m:
                uv_list.append([ev, float(os_m.group(1))])
            if rs_m:
                ecd_list.append([ev, float(rs_m.group(1))])

        if uv_list:
            uv = np.array(uv_list, dtype=np.float64)
        if ecd_list:
            ecd = np.array(ecd_list, dtype=np.float64)

        return uv, ecd

    def opt_done(self) -> bool:
        """Check if geometry optimisation converged successfully.

        Returns:
            bool: True if converged, False otherwise.
        """
        return len(re.findall(self.regex["opt_done"], self.fl)) >= 1

    def normal_termination(self) -> bool:
        """Check if the NWChem calculation terminated normally.

        Returns:
            bool: True if normal termination is detected.
        """
        return len(re.findall(self.regex["finish"], self.fl)) >= 1




if __name__ == '__main__': 
    from mock import MagicMock
    p = NWChemParser('water_nwchem.nwo', log=MagicMock())

    print(p.parse_geom())
    print(p.parse_B_m())
    print(p.opt_done())