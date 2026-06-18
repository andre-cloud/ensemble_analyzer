__all__ = [
    "CALCULATOR_REGISTRY", "BaseCalc", "BaseMlCalc",
    "ML_CALCULATORS", "register_calculator",
    "get_ase_calculator", "convert_to_single_point", "get_models_dir",
]

def __getattr__(name):
    if name in ("CALCULATOR_REGISTRY", "BaseCalc", "ML_CALCULATORS", "register_calculator"):
        from .base import CALCULATOR_REGISTRY, BaseCalc, ML_CALCULATORS, register_calculator
        g = globals()
        g["CALCULATOR_REGISTRY"] = CALCULATOR_REGISTRY
        g["BaseCalc"] = BaseCalc
        g["ML_CALCULATORS"] = ML_CALCULATORS
        g["register_calculator"] = register_calculator
        return g[name]
    if name == "BaseMlCalc":
        from ._ml_base import BaseMlCalc
        globals()["BaseMlCalc"] = BaseMlCalc
        return BaseMlCalc
    if name in ("get_ase_calculator", "convert_to_single_point", "get_models_dir"):
        from enan_calculators import get_ase_calculator, convert_to_single_point, get_models_dir
        g = globals()
        g["get_ase_calculator"] = get_ase_calculator
        g["convert_to_single_point"] = convert_to_single_point
        g["get_models_dir"] = get_models_dir
        return g[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
