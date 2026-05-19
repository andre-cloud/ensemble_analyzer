from abc import ABC, abstractmethod
from typing import Callable, Dict, List, Tuple
from ase import Atoms
from ensemble_analyzer._conformer.conformer import Conformer
from ensemble_analyzer.constants import ROT_CONST_FACTOR

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
    """
    Abstract Base Class for output parsers.
    
    This class defines the interface that any new QM software parser must implement
    to be compatible with Ensemble Analyzer.
    """

    def __init__(self, output_name: str, log: 'Logger', conf: Conformer) -> None:
        """Initialize the parser.

        Args:
            output_name: Path to the output file to parse.
            log: Logger instance for warnings and debug info.
        """

        with open(output_name) as f:
            self.fl = f.read()

        self.log = log
        self.conf = conf
        self.skip_message = "ATTENTION: Calculation CRASHED, impossible parsing. Conformer will be deactivated and no longer considered"
    
    @abstractmethod
    def parse_geom(self) -> np.ndarray:
        """Extract the final geometry from the output.

        Returns:
            np.ndarray: Array of shape (N_atoms, 3) containing Cartesian coordinates.
        """
        pass

    @abstractmethod
    def parse_B_m(self) -> Tuple[np.ndarray, np.ndarray]:
        """Extract Rotational Constants and Dipole Moment.

        Returns:
            Tuple[np.ndarray, np.ndarray]:
                - Rotational constants vector (B_vec).
                - Dipole moment vector (M_vec).
        """
        pass

    @abstractmethod
    def parse_energy(self) -> float:
        """Extract the final electronic energy.

        Returns:
            float: Electronic energy in Hartree (Eh).
        """
        pass
    
    @abstractmethod
    def parse_freq(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Extract vibrational frequencies and spectral data (IR/VCD).

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]:
                - Array of frequencies [cm^-1].
                - IR spectrum data (X, Y).
                - VCD spectrum data (X, Y).
        """
        pass

    @abstractmethod
    def parse_tddft(self) -> Tuple[np.ndarray, np.ndarray]:
        """Extract TD-DFT excited states data (UV/ECD).

        Returns:
            Tuple[np.ndarray, np.ndarray]:
                - UV spectrum data (Energy, Intensity).
                - ECD spectrum data (Energy, Rotational Strength).
        """
        pass

    @abstractmethod
    def opt_done(self) -> bool:
        """
        Check if geometry optimization converged successfully.

        Returns:
            bool: True if converged, False otherwise.
        """
        pass

    @abstractmethod
    def normal_termination(self) -> bool:
        """
        Check if the calculation terminated normally.

        Returns:
            bool: True if normal termination is detected.
        """
        pass
    
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