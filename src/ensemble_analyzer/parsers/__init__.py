import importlib
import pkgutil

# Automatically import all modules in the directory
for module_info in pkgutil.iter_modules(__path__):
    if module_info.name != "base":  # avoid reloading base.py
        importlib.import_module(f"{__name__}.{module_info.name}")

# Expose the global registry
from .base import PARSER_REGISTRY, BaseParser, register_parser

__all__ = ["PARSER_REGISTRY", "BaseParser", "register_parser"]