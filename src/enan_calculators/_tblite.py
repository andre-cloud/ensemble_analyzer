try:
    from tblite.ase import TBLite as TB
except ImportError:
    TB = None


def create_tblite_calc(charge, mult, method, solvent=None):
    if TB is None:
        raise ImportError("tblite module missing. Install via: pip install tblite")

    solv = None
    if solvent:
        solv = ("alpb", solvent)

    return TB(
        method=method,
        charge=charge,
        multiplicity=mult,
        solvation=solv,
    )
