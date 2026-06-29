from enan_calculators._models import get_models_dir


def get_ase_calculator(engine, charge, mult, method, solvent=None, **kwargs):
    engine = engine.lower()
    if engine == "aimnet":
        from enan_calculators._aimnet import create_aimnet_calc
        return create_aimnet_calc(charge, mult, method, solvent=solvent)
    elif engine == "tblite":
        from enan_calculators._tblite import create_tblite_calc
        return create_tblite_calc(charge, mult, method, solvent=solvent)
    elif engine == "uma":
        from enan_calculators._uma import create_uma_calc
        return create_uma_calc(charge, mult, method, solvent=solvent)
    elif engine == "fairchem":
        from enan_calculators._fairchem import create_fairchem_calc
        return create_fairchem_calc(charge, mult, method, solvent=solvent)
    elif engine == "mace":
        from enan_calculators._mace import create_mace_calc
        return create_mace_calc(charge, mult, method, solvent=solvent)
    elif engine == "orca":
        from enan_calculators._orca import create_orca_calc
        return create_orca_calc(
            charge, mult, method, kwargs.get("basis"),
            solvent=solvent, cpu=kwargs.get("cpu", 1),
            add_input=kwargs.get("add_input", ""),
            directory=kwargs.get("directory"),
        )
    elif engine == "gaussian":
        from enan_calculators._gaussian import create_gaussian_calc
        return create_gaussian_calc(
            charge, mult, method, kwargs.get("basis"),
            solvent=solvent, cpu=kwargs.get("cpu", 1),
            add_input=kwargs.get("add_input", ""),
            label=kwargs.get("label"),
        )
    elif engine == "nwchem":
        from enan_calculators._nwchem import create_nwchem_calc
        return create_nwchem_calc(
            charge, mult, method, kwargs.get("basis"),
            solvent=solvent, cpu=kwargs.get("cpu", 1),
            add_input=kwargs.get("add_input", ""),
            label=kwargs.get("label"),
            command=kwargs.get("command"),
        )
    elif engine == "skala":
        from enan_calculators._skala import create_skala_calc
        return create_skala_calc(
            charge, mult, method, kwargs.get("basis"), solvent=solvent,
        )
    raise ValueError(f"Unknown calculator engine: {engine}")


def convert_to_single_point(atoms):
    atoms_copy = atoms.copy()
    atoms_copy.calc = None
    return atoms_copy


__all__ = [
    "get_ase_calculator",
    "convert_to_single_point",
    "get_models_dir",
]
