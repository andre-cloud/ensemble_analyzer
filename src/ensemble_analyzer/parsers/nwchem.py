from ensemble_analyzer.parsers.base import BaseParser, register_parser

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

    def parse_geom(self) -> np.ndarray:
        """Parse the final geometry from NWChem output.

        Returns:
            np.ndarray: Cartesian coordinates array of shape (N, 3).
        """
        fl = self.get_filtered_text(start=self.regex["geom_start"], end="Atomic Mass")

        pattern = r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
        coords = np.array(re.findall(pattern, fl, flags=re.MULTILINE), dtype=float)
        return coords



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

    def parse_normal_modes(self, n_atoms: int) -> np.ndarray:
        if self.regex["s_freq"] not in self.fl:
            return np.empty((0, n_atoms, 3))
        try:
            section = self.fl.split(self.regex["s_freq"])[-1]
            section = section.split("\n\n\n")[0]
            lines = section.strip().split("\n")
            n_total = 3 * n_atoms
            max_col = -1
            for line in lines:
                stripped = line.strip()
                if not stripped or stripped.startswith("---"):
                    continue
                header_match = re.match(r'^\s*(\d+(?:\s+\d+)*)\s*$', stripped)
                if header_match:
                    cols = [int(x) for x in header_match.group(1).split()]
                    if cols and cols[-1] > max_col:
                        max_col = cols[-1]
            if max_col < 0:
                return np.empty((0, n_atoms, 3))
            matrix = np.zeros((n_total, max_col))
            current_columns = []
            for line in lines:
                stripped = line.strip()
                if not stripped or stripped.startswith("---"):
                    continue
                header_match = re.match(r'^\s*(\d+(?:\s+\d+)*)\s*$', stripped)
                if header_match:
                    current_columns = [int(x) for x in header_match.group(1).split()]
                    continue
                if stripped.startswith("P.Frequency"):
                    continue
                data_match = re.match(
                    r'^\s*(\d+)\s+((?:[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?\s*)+)', stripped
                )
                if data_match and current_columns:
                    row_idx = int(data_match.group(1)) - 1
                    coeff_str = data_match.group(2)
                    coeffs = [float(c) for c in re.findall(
                        r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', coeff_str
                    )]
                    for i, c in enumerate(coeffs):
                        if i < len(current_columns) and row_idx < n_total:
                            col = current_columns[i] - 1
                            if col < max_col:
                                matrix[row_idx, col] = c
            n_modes = matrix.shape[1]
            all_modes = np.zeros((n_modes, n_atoms, 3))
            for mode_idx in range(n_modes):
                vec = matrix[:n_total, mode_idx]
                all_modes[mode_idx] = vec.reshape(n_atoms, 3)
            return all_modes
        except Exception:
            return np.empty((0, n_atoms, 3))

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






if __name__ == '__main__': 
    from mock import MagicMock
    p = NWChemParser('water_nwchem.nwo', log=MagicMock())

    print(p.parse_geom())
    print(p.parse_B_m())
    print(p.opt_done())