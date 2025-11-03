import importlib
import pkgutil

# Importa automaticamente tutti i moduli nella cartella
for module_info in pkgutil.iter_modules(__path__):
    if module_info.name != "base":  # evita di ricaricare base.py
        importlib.import_module(f"{__name__}.{module_info.name}")

# Espone il registro globale
from .base import CALCULATOR_REGISTRY, BaseCalc, register_calculator

__all__ = ["CALCULATOR_REGISTRY", "BaseCalc", "register_calculator"]