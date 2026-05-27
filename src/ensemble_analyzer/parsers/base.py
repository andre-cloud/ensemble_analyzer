from abc import ABC, abstractmethod
from typing import Callable, Dict, List, Tuple
from ase import Atoms
from ensemble_analyzer.conformer.conformer import Conformer
from ensemble_analyzer.constants import ROT_CONST_FACTOR, CONVERT_B

import re
import numpy as np

def register_parser(name: str) -> Callable:
    """Decorator to register each parser in the global registry.

    Args:
        name: Parser identifier (e.g. 'gaussian', 'orca').

    Returns:
        Callable: Decorator that registers the class in PARSER_REGISTRY.
    """

    def decorator(cls: type) -> type:
        """Register a parser class.

        Args:
            cls: Parser class to register.

        Returns:
            type: The input class unchanged.
        """
        PARSER_REGISTRY[name.lower()] = cls
        return cls

    return decorator

class BaseParser(ABC):
    REGEX: Dict = {}

    def __init__(self, output_name: str, log, conf: Conformer = None) -> None:
        with open(output_name) as f:
            self.fl = f.read()
        self.log = log
        self.conf = conf
        self.skip_message = "ATTENTION: Calculation CRASHED, impossible parsing. Conformer will be deactivated and no longer considered"
        self._init_regex()
        self.correct_exiting = self.normal_termination()
        if not self.correct_exiting:
            self.log.warning(self.skip_message)

    def _init_regex(self) -> None:
        self.regex = self.REGEX

    @abstractmethod
    def parse_geom(self) -> np.ndarray:
        pass

    def parse_B_m(self) -> Tuple[np.ndarray, np.ndarray]:
        match_B = re.findall(self.regex['B'], self.fl)
        if match_B:
            B = np.array(match_B[-1], dtype=float)
            if self.regex['units_B'] != 'cm-1':
                B /= CONVERT_B[self.regex['units_B']]
        else:
            self.log.warning(f"\t{self.log.WARNING} B not found, calculating with ASE")
            B = self.calculate_B()

        dipole_text = self._get_dipole_text()
        match_M = re.findall(self.regex['m'], dipole_text)
        if match_M:
            M = np.array(match_M[-1], dtype=float)
        else:
            self.log.warning(f"\t{self.log.WARNING} M not found, storing a versor")
            M = np.array([1, 0, 0])
        return B, M

    def _get_dipole_text(self) -> str:
        return self.fl

    def parse_energy(self) -> float:
        match = re.findall(self.regex['E'], self.fl)
        if not match:
            return 0.0
        return float(match[-1])

    @abstractmethod
    def parse_freq(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        pass

    def parse_normal_modes(self, n_atoms: int) -> np.ndarray:
        return np.empty((0, n_atoms, 3))

    @abstractmethod
    def parse_tddft(self) -> Tuple[np.ndarray, np.ndarray]:
        pass

    def opt_done(self) -> bool:
        return len(re.findall(self.regex['opt_done'], self.fl)) >= 1

    def normal_termination(self) -> bool:
        return len(re.findall(self.regex['finish'], self.fl)) >= 1

    def get_filtered_text(self, start:str, end:str) -> str:
        """
        Extract a section of text between two delimiters.

        Args:
            start (str): Start delimiter.
            end (str): End delimiter.

        Returns:
            str: The extracted text block.
        """
        return self.fl.split(start)[-1].split(end)[0]
    

    def parse_table(self, table:list, list_index:list) -> List[List[str]]:
        """
        Parse a fixed-width or space-separated table from text lines.

        Args:
            table (List[str]): List of strings representing the table rows.
            list_index (List[int]): Indices of columns to extract.

        Returns:
            List[List[str]]: Extracted data matrix.
        """
        data = []
        for line in table: 
            if not line: 
                continue
            if '---' in line: 
                continue
            line_splitted = line.split()
            data.append([line_splitted[i] for i in list_index])
        
        return data

    def calculate_B(self) -> np.ndarray:
        """Compute principal rotational constants from conformer geometry.

        Uses the conformer's atomic positions and masses to calculate
        the three principal rotational constants.

        Returns:
            np.ndarray: Array of shape (3,) containing B_a, B_b, B_c in cm⁻¹.
            Falls back to [1, 0, 0] if conformer data is unavailable.
        """

        if self.conf is None or self.conf.last_geometry is None:
            self.log.warning("No conformer data for B calculation, returning default")
            return np.array([1.0, 0.0, 0.0])

        self.log.debug(f'{self.conf.last_geometry = }')
        atoms = Atoms(
            symbols="".join(tuple(self.conf.atoms)),
            positions=self.conf.last_geometry,
        )
        moments = atoms.get_moments_of_inertia()  # [amu·Å²]
        B_vec = np.divide(
            ROT_CONST_FACTOR, moments,
            out=np.zeros_like(moments),
            where=moments > 1e-6,
        )  # [cm⁻¹]
        return B_vec


PARSER_REGISTRY : Dict[str, BaseParser]= {}