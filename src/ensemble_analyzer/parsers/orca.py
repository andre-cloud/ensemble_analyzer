
from ensemble_analyzer.parsers.base import BaseParser, register_parser
from ensemble_analyzer.constants import FACTOR_EV_CM_1
from ensemble_analyzer.parsers.vibrations import orca_get_freq, orca_get_normal_modes, to_displacements

import re
import numpy as np
from typing import Tuple


@register_parser('orca')
class OrcaParser(BaseParser):
    """
    Parser implementation for ORCA output files.
    Supports ORCA 5 and ORCA 6 syntax variations (older versions use the ORCA 5 layout).
    """

    REGEX = {
        '6':{
        "B": r"Rotational constants in cm-1:\s*(-?\d+.\d*)\s*(-?\d+.\d*)\s*(-?\d+.\d*)",
        'units_B': 'cm-1',
        "m": r"Total Dipole Moment\s*:\s*([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)",
        "E": r"FINAL SINGLE POINT ENERGY.*?(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
        "s_IR": """Mode   freq       eps      Int      T**2         TX        TY        TZ
       cm**-1   L/(mol*cm) km/mol    a.u.
----------------------------------------------------------------------------
""",
        "s_VCD": """Mode   Freq    VCD-Intensity    
       (1/cm) (1E-44*esu^2*cm^2) 
---------------------------------""",
        "break": "\n\n",
        "idx_en_ir": 1,  # index for energy in the IR table in cm**-1
        "idx_imp_ir": 3,  # index for oscillator strength in the IR table
        "idx_en_vcd": 1,  # index for energy in the VCD table in cm**-1
        "idx_imp_vcd": 2,  # index for oscillator strength in the VCD table
        "s_freq": "VIBRATIONAL FREQUENCIES",
        "e_freq": "------------",
        "idx_freq": 1,  # index for frequency in frequency table
        # frequency / normal-mode blocks (same keys as frogchem regex_ORCA)
        "frequency": r'VIBRATIONAL FREQUENCIES([\S\s]*?)\n\n\n',
        "normal_modes": r'NORMAL MODES([\S\s]*?)\n\n\n',
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
        "s_IR": """Mode   freq       eps      Int      T**2         TX        TY        TZ
       cm**-1   L/(mol*cm) km/mol    a.u.
----------------------------------------------------------------------------
""",
        "s_VCD": None,
        "break": "\n\n",
        "idx_en_ir": 1,  # index for energy in the IR table in cm**-1
        "idx_imp_ir": 3,  # index for oscillator strength in the IR table
        "idx_imp_vcd": 2,  # index for oscillator strength in the VCD table
        "s_freq": "VIBRATIONAL FREQUENCIES",
        "e_freq": "------------",
        "idx_freq": 1,  # index for frequency in frequency table
        # frequency / normal-mode blocks (same keys as frogchem regex_ORCA)
        "frequency": r'VIBRATIONAL FREQUENCIES([\S\s]*?)\n\n\n',
        "normal_modes": r'NORMAL MODES([\S\s]*?)\n\n\n',
        "opt_done": "THE OPTIMIZATION HAS CONVERGED",
        "geom_start": """CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------""",
        "finish": "ORCA TERMINATED NORMALLY",
        "ext": "out",}
    }

    # --- TD-DFT spectra -------------------------------------------------
    # Titles of the tables (matched as a whole line, so "SOC CORRECTED ..." and
    # "... VELOCITY DIPOLE MOMENTS" blocks are never confused with these ones).
    TITLE_UV = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
    TITLES_ECD = (
        "CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS",  # ORCA >= 6
        "CD SPECTRUM",                                        # ORCA <= 5
    )
    SOC_PREFIX = "SOC CORRECTED "
    # ORCA >= 6 rows:  "  0-1A  ->  1-1A    2.778325   22408.7   446.3   0.079978963 ..."
    _ROW_EV = re.compile(
        r"^\s*\S+\s+->\s+\S+\s+(-?\d+\.\d+)\s+\S+\s+\S+\s+(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)"
    )
    # ORCA <= 5 rows:  "   1   43426.1    230.3   0.005999216 ..."  (energy in cm-1)
    _ROW_CM = re.compile(
        r"^\s*\d+\s+(-?\d+\.\d+)\s+\S+\s+(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)"
    )

    # --- TD-DFT energy summary -------------------------------------------
    # "    DE(CIS) =      0.065909138 Eh (Root  1)"
    # "    E(SOC CIS)    =      0.092111444 Eh (SOC Root  1)"
    _RE_EXCITATION = re.compile(
        r"(?:DE\(CIS\)|E\(SOC CIS\))\s*=\s*(-?\d+\.\d+)\s*Eh\s*\((?:SOC\s+)?Root\s+(\d+)\)"
    )
    _RE_FINAL = re.compile(r"FINAL SINGLE POINT ENERGY")

    _soc_notice_done = False

    def _init_regex(self) -> None:
        self.version = self.get_version()
        major = int(self.version)
        if major == 0:
            self.log.warning("ORCA version not detected in the output: assuming the ORCA 6 layout")
        self.regex = self.REGEX['5' if 0 < major < 6 else '6']

    def get_version(self) -> str:
        """
        Detect ORCA major version from output header.

        Returns:
            str: Major version (e.g. '5' or '6'), '0' when not found.
        """
        find = re.findall(r'Program Version (\d+)', self.fl)
        return find[0] if find else "0"

    def parse_energy(self, excited_state: bool = False) -> float:
        """
        Parse the final single point energy.

        In a TD-DFT/CIS run ORCA writes as ``FINAL SINGLE POINT ENERGY`` the energy
        of the followed root (E(SCF) + DE(CIS) [+ dispersion]). Unless
        ``excited_state`` is True, that excitation is removed so the ground-state
        energy (dispersion included) is returned, consistently with the other QM
        parsers. Outputs computed with ``iroot 0`` already report the ground state.

        Args:
            excited_state (bool): keep the energy of the followed excited root.

        Returns:
            float: Energy in Eh (0.0 if not found).
        """
        energy = super().parse_energy()
        if excited_state or not energy:
            return energy

        finals = [m.start() for m in self._RE_FINAL.finditer(self.fl)]
        if not finals:
            return energy
        lower = finals[-2] if len(finals) > 1 else -1
        upper = finals[-1]
        excitations = [
            m for m in self._RE_EXCITATION.finditer(self.fl)
            if lower < m.start() < upper
        ]
        if not excitations:
            return energy

        delta, root = excitations[-1].groups()
        if int(root) == 0:
            return energy
        self.log.debug(
            f"TD-DFT output: removing the excitation energy of root {root} "
            f"({float(delta):.6f} Eh) from the final energy"
        )
        return energy - float(delta)

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

        # freq (FROGChem parser): translations/rotations are printed as zeros and dropped here
        freq = self._all_freq()
        freq = freq[freq != 0]

        # IR
        if self.regex['s_IR'] not in self.fl:
            ir = np.zeros(shape=(1,2))
        else:
            ir_text = self.get_filtered_text(start=self.regex['s_IR'], end='\n\n').splitlines()
            ir = np.array(self.parse_table(ir_text, [self.regex['idx_en_ir'], self.regex['idx_imp_ir']]), dtype=np.float64)

        # VCD
        if self.regex['s_VCD'] is None or self.regex['s_VCD'] not in self.fl:
            vcd = np.zeros(shape=(1,2))
        else:
            vcd_text = self.get_filtered_text(start=self.regex['s_VCD'], end='\n\n').splitlines()
            vcd = np.array(self.parse_table(vcd_text, [self.regex['idx_en_vcd'], self.regex['idx_imp_vcd']]), dtype=np.float64)

        return freq, ir, vcd

    def _all_freq(self) -> np.ndarray:
        """All ORCA frequencies (cm-1), translations/rotations included as zeros."""
        if not re.search(self.regex['frequency'], self.fl):
            return np.array([])
        return orca_get_freq(self.fl, self.regex)[1]

    def parse_normal_modes(self, n_atoms: int) -> np.ndarray:
        """Cartesian displacement vectors of the vibrational modes (FROGChem parser).

        Returns:
            np.ndarray: (n_modes, n_atoms, 3), same order as the frequencies of parse_freq
            (translation/rotation modes, i.e. zero frequencies, are dropped). Empty
            (0, n_atoms, 3) array when the NORMAL MODES block is missing.
        """
        empty = np.empty((0, n_atoms, 3))
        try:
            matrix = orca_get_normal_modes(self.fl, self.regex, N=3 * n_atoms)
        except ValueError:
            return empty
        except Exception as err:
            self.log.warning(f"Normal modes could not be parsed: {err}")
            return empty
        modes = to_displacements(matrix, n_atoms)
        freq = self._all_freq()
        if len(freq) == len(modes):
            modes = modes[freq != 0]
        return modes

    def _parse_spectrum_block(self, title: str) -> np.ndarray:
        """
        Parse the last TD-DFT table whose title line is exactly ``title``.

        The layout is inferred from the table header: ORCA >= 6 tables have an
        energy column in eV, older ones only in cm-1 (converted to eV here).

        Returns:
            np.ndarray: (n_transitions, 2) array [energy (eV), intensity];
            empty (0, 2) when the table is missing.
        """
        pattern = re.compile(rf"^[ \t]*{re.escape(title)}[ \t]*$", re.MULTILINE)
        matches = list(pattern.finditer(self.fl))
        if not matches:
            return np.empty((0, 2))

        lines = self.fl[matches[-1].end():].splitlines()
        in_ev = "(eV)" in "\n".join(lines[:8])
        row_re = self._ROW_EV if in_ev else self._ROW_CM

        data = []
        for i, line in enumerate(lines):
            match = row_re.match(line)
            if match:
                data.append((float(match.group(1)), float(match.group(2))))
            elif data or i > 10:
                break

        spectrum = np.array(data, dtype=np.float64).reshape(-1, 2)
        if not in_ev and spectrum.size:
            spectrum[:, 0] *= FACTOR_EV_CM_1
        return spectrum

    def _parse_spectrum(self, title: str) -> np.ndarray:
        """SOC-corrected table when present (``dosoc true``), plain one otherwise."""
        soc = self._parse_spectrum_block(self.SOC_PREFIX + title)
        if soc.size:
            if not OrcaParser._soc_notice_done:
                self.log.info(
                    "ORCA SOC-corrected spectra found: using them. They are valid within the "
                    "computed roots only (increase nroots if needed); each triplet contributes "
                    "three sublevels."
                )
                OrcaParser._soc_notice_done = True
            return soc
        return self._parse_spectrum_block(title)

    def parse_tddft(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Parse TD-DFT excited states for UV and ECD spectra.

        Returns:
            Tuple: (UV data array, ECD data array), each [energy (eV), intensity].
            A (1, 2) zero array is returned for a missing spectrum.
        """
        uv = self._parse_spectrum(self.TITLE_UV)

        ecd = np.empty((0, 2))
        for title in self.TITLES_ECD:
            ecd = self._parse_spectrum(title)
            if ecd.size:
                break

        if uv.size == 0:
            uv = np.zeros(shape=(1, 2))
        if ecd.size == 0:
            ecd = np.zeros(shape=(1, 2))
        return uv, ecd
