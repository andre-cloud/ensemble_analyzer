# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Ensemble Analyzer'
copyright = '2025, Andrea Pellegrini, Paolo Righi'
author = 'Andrea Pellegrini, Paolo Righi'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

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
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']


myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "colon_fence",
]

import mock
import os, sys
sys.path.insert(0, os.path.abspath('../src'))

class Mock(mock.MagicMock):
    @classmethod
    def __getattr__(cls,name): 
        return mock.MagicMock()
    
MOCK_MODULES = [
    'numpy', 
    'scipy', 
    'scipy.spatial',
    'scipy.spatial.distance',
    'scipy.optimize',
    'scipy.constants',
    'scipy.interpolate',
    'matplotlib', 
    'matplotlib.pyplot', 
    'ase', 
    'ase.io',
    'ase.units',
    'numba', 
    'sklearn',
    'sklearn.cluster'
]

sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)