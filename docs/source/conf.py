# Configuration file for the Sphinx documentation builder.
import os
import sys
from unittest.mock import MagicMock

# -- Path setup --------------------------------------------------------------
# docs/source/conf.py -> ../../src
sys.path.insert(0, os.path.abspath('../../src'))


# -- Project information -----------------------------------------------------
project = 'Ensemble Analyzer'
copyright = '2025, Andrea Pellegrini, Paolo Righi'
author = 'Andrea Pellegrini, Paolo Righi'


# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'myst_parser',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']


# -- MyST Parser configuration -----------------------------------------------
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "colon_fence",
]


# -- ROBUST MOCKING CONFIGURATION --------------------------------------------

# 1. Pre-inject modules to fix "constants.py" math errors
# We assign real numbers to constants so (h * c) works during import
MOCK_CONSTANTS = MagicMock()
MOCK_CONSTANTS.h = 1.0
MOCK_CONSTANTS.c = 1.0
MOCK_CONSTANTS.electron_volt = 1.0
MOCK_CONSTANTS.Boltzmann = 1.0
MOCK_CONSTANTS.N_A = 1.0
MOCK_CONSTANTS.R = 1.0
# Mock dictionary for physical_constants lookups
MOCK_CONSTANTS.physical_constants = {
    "joule-hartree relationship": (1.0, '', 0.0),
    "atomic mass constant": (1.0, '', 0.0)
}
sys.modules['scipy.constants'] = MOCK_CONSTANTS

# 2. Pre-inject modules to fix "from submodule import Class" errors
sys.modules['matplotlib.figure'] = MagicMock()
sys.modules['ase.atoms'] = MagicMock()
sys.modules['ase.calculators.gaussian'] = MagicMock()
sys.modules['ase.calculators.orca'] = MagicMock()

# 3. Standard mocking for the rest
autodoc_mock_imports = [
    'numpy',
    'scipy',
    'matplotlib',
    'ase',
    'numba',
    'sklearn',
    'InquirerPy',
    'tabulate',
    'mock'
]