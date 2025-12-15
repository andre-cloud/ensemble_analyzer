from abc import ABC, abstractmethod
from typing import Dict, List

def register_parser(name):
    """Decorator to register each parser in the global registry."""

    def decorator(cls):
        PARSER_REGISTRY[name.lower()] = cls
        return cls

    return decorator

class BaseParser(ABC):
    """
    Abstract Base Class for output parsers.
    
    This class defines the interface that any new QM software parser must implement
    to be compatible with Ensemble Analyzer.
    """

    def __init__(self, output_name, log): 
        """
        Initialize the parser.

        Args:
            output_name (str): Path to the output file to parse.
            log (Logger): Logger instance for warnings and debug info.
        """

        with open(output_name) as f:
            self.fl = f.read()

        self.log = log

        self.skip_message = "ATTENTION: Calculation CRASHED, impossible parsing. Conformer will be deactivated and no longer considered"
    
    @abstractmethod
    def parse_geom(self):
        """
        Extract the final geometry from the output.

        Returns:
            np.ndarray: Array of shape (N_atoms, 3) containing Cartesian coordinates.
        """
        pass

    @abstractmethod
    def parse_B_m(self):
        """
        Extract Rotational Constants and Dipole Moment.

        Returns:
            Tuple[np.ndarray, np.ndarray]:
                - Rotational constants vector (B_vec).
                - Dipole moment vector (M_vec).
        """
        pass

    @abstractmethod
    def parse_energy(self):
        """
        Extract the final electronic energy.

        Returns:
            float: Electronic energy in Hartree (Eh).
        """
        pass
    
    @abstractmethod
    def parse_freq(self):
        """
        Extract vibrational frequencies and spectral data (IR/VCD).

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]:
                - Array of frequencies [cm^-1].
                - IR spectrum data (X, Y).
                - VCD spectrum data (X, Y).
        """
        pass

    @abstractmethod
    def parse_tddft(self):
        """
        Extract TD-DFT excited states data (UV/ECD).

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



PARSER_REGISTRY : Dict[str, BaseParser]= {}