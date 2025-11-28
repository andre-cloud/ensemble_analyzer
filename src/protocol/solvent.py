



class Solvent:
    """
    Solvent class
    """

    def __init__(self, solv: dict):
        self.solvent = solv["solvent"]
        self.smd = solv["smd"]

    def __str__(self):  
        if self.smd:
            return f"SMD({self.solvent})"
        elif self.solvent:
            return f"CPCM({self.solvent})"
        else:
            return "CPCM"

    def __repr__(self): 
        return self.__str__()
