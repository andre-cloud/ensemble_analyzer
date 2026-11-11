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
            self._freq_mask = mask.copy()
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
            if hasattr(self, '_freq_mask') and self._freq_mask is not None and len(self._freq_mask) == len(all_modes):
                all_modes = all_modes[self._freq_mask]
            return all_modes
        except Exception:
            return np.empty((0, n_atoms, 3))

    # --- ΔSCF support --------------------------------------------------------

    def parse_s2(self) -> float | None:
        """Last ground-state <S2> of an unrestricted DFT run, None if not printed."""
        values = re.findall(r"^\s+<S2> =\s+(-?\d+\.\d+)\s+\(Exact", self.fl, flags=re.MULTILINE)
        return float(values[-1]) if values else None

    def parse_n_electrons(self) -> tuple[int, int] | None:
        """(n_alpha, n_beta) of the last DFT run (valence electrons with ECPs)."""
        alpha = re.findall(r"Alpha electrons :\s+(\d+)", self.fl)
        beta = re.findall(r"Beta electrons :\s+(\d+)", self.fl)
        if not alpha or not beta:
            return None
        return int(alpha[-1]), int(beta[-1])

    def parse_excited_states(self, mult: int = 1) -> list:
        """TD-DFT roots with their orbital composition (1-based indices).
        For a closed-shell reference only the singlet roots are returned."""
        from ensemble_analyzer.protocol.deltascf import ExcitedState

        header_re = re.compile(
            r"^\s*Root\s+(\d+)\s+(singlet|triplet)?\s*\S*\s+(-?\d+\.\d+)\s+a\.u\.\s+(-?\d+\.\d+)\s+eV"
        )
        trans_re = re.compile(
            r"^\s*Occ\.\s+(\d+)\s+(?:(alpha|beta)\s+)?\S+\s+---\s+Virt\.\s+(\d+)\s+"
            r"(?:(?:alpha|beta)\s+)?\S+\s+(-?\d+\.\d+)(?:\s+([XY]))?"
        )
        states: dict = {}
        current = None
        for line in self.fl.splitlines():
            header = header_re.match(line)
            if header:
                kind = header.group(2)
                root = int(header.group(1))
                if (mult == 1 and kind == "triplet") or root in states:
                    current = None
                    continue
                current = ExcitedState(root=root, energy_ev=float(header.group(4)), transitions=[])
                states[root] = current
                continue
            if current is None:
                continue
            trans = trans_re.match(line)
            if trans and trans.group(5) != "Y":
                spin = trans.group(2) or "alpha"
                coef = float(trans.group(4))
                current.transitions.append((int(trans.group(1)), int(trans.group(3)), spin, coef * coef))
        return [states[k] for k in sorted(states)]

    def deltascf_status(self, opt: bool = False) -> str | None:
        """None if the run converged, otherwise the reason it failed."""
        if "Calculation failed to converge" in self.fl:
            return "SCF not converged"
        if "dft energy failed" in self.fl or "dft optimize failed" in self.fl:
            return "NWChem error"
        if not self.correct_exiting:
            return "NWChem did not terminate normally"
        if opt and not self.opt_done():
            return "geometry optimisation not converged"
        return None

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