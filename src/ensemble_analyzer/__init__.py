from .api import run, EnsembleResult, load_ensemble, save_ensemble
from .launch import main
from .conformer.conformer import Conformer
from .protocol.protocol import Protocol, load_protocol
from ._managers.calculation_config import CalculationConfig
from .ensemble_io import read_ensemble, save_snapshot
import importlib.metadata

try:
    __version__ = importlib.metadata.version("ensemble-analyzer")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"

__all__ = [
    "run",
    "main",
    "EnsembleResult",
    "Conformer",
    "Protocol",
    "CalculationConfig",
    "load_protocol",
    "load_ensemble",
    "save_ensemble",
    "read_ensemble",
    "save_snapshot",
    "__version__",
]
