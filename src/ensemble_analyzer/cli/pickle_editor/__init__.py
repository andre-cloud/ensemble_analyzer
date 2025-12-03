"""
Matplotlib Pickle Editor Package

Module for interactive and batch editing of serialized matplotlib figures.
"""

from .core import MatplotlibPickleEditor, PickleSecurityError
from .tui import InteractiveTUI

__version__ = "1.0.0"
__all__ = ["MatplotlibPickleEditor", "PickleSecurityError", "InteractiveTUI"]