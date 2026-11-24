from .base import BaseParser
from .gaussian import GaussianParser
from .nwchem import NWChemParser
from .orca import OrcaParser

PARSER_REGISTRY: dict[str, type[BaseParser]] = {
    "orca": OrcaParser,
    "gaussian": GaussianParser,
    "nwchem": NWChemParser,
}

__all__ = ["PARSER_REGISTRY", "BaseParser"]
