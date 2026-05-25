import importlib
import pkgutil

# Automatically import all modules in the directory
for module_info in pkgutil.iter_modules(__path__):
    if module_info.name not in ("base", "_ml_base"):  # avoid reloading base modules
        importlib.import_module(f"{__name__}.{module_info.name}")

# Expose the global registry and base classes
from .base import CALCULATOR_REGISTRY, BaseCalc, ML_CALCULATORS, register_calculator
from ._ml_base import BaseMlCalc

__all__ = ["CALCULATOR_REGISTRY", "BaseCalc", "BaseMlCalc", "ML_CALCULATORS", "register_calculator"]
