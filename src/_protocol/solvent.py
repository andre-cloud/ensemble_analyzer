
from dataclasses import dataclass

@dataclass
class Solvent:
    """
    Solvent class
    """

    solvent : str
    smd     : bool  = False 

    def __repr__(self):  
        if self.smd:
            return f"SMD({self.solvent})"
        elif self.solvent:
            return f"CPCM({self.solvent})"
        else:
            return "CPCM"