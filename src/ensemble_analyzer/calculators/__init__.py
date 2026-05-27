import importlib
import pkgutil

_loaded = False

def _ensure_loaded():
    global _loaded
    if _loaded:
        return
    for module_info in pkgutil.iter_modules(__path__):
        if module_info.name not in ("base", "_ml_base"):
            importlib.import_module(f"{__name__}.{module_info.name}")
    _loaded = True

def __getattr__(name):
    if name in ("CALCULATOR_REGISTRY", "BaseCalc", "ML_CALCULATORS", "register_calculator"):
        _ensure_loaded()
        from .base import CALCULATOR_REGISTRY, BaseCalc, ML_CALCULATORS, register_calculator
        g = globals()
        g["CALCULATOR_REGISTRY"] = CALCULATOR_REGISTRY
        g["BaseCalc"] = BaseCalc
        g["ML_CALCULATORS"] = ML_CALCULATORS
        g["register_calculator"] = register_calculator
        return g[name]
    if name == "BaseMlCalc":
        _ensure_loaded()
        from ._ml_base import BaseMlCalc
        globals()["BaseMlCalc"] = BaseMlCalc
        return BaseMlCalc
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = ["CALCULATOR_REGISTRY", "BaseCalc", "BaseMlCalc", "ML_CALCULATORS", "register_calculator"]
