
from dataclasses import dataclass

@dataclass
class Solvent:
    """
    Data container for implicit solvation settings.
    """

    solvent : str
    smd     : bool  = False 

    def __repr__(self) -> str:
        """
        Generate a string representation of the solvent model.

        Returns:
            str: e.g. "SMD(water)" or "CPCM(chloroform)".
        """
        
        if self.smd:
            return f"SMD({self.solvent})"
        elif self.solvent:
            return f"CPCM({self.solvent})"
        else:
            return "CPCM"