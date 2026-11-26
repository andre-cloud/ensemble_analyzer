# Configuration file for the Sphinx documentation builder.
import os
import sys

# -- Path setup --------------------------------------------------------------
# docs/source/conf.py -> ../../src
sys.path.insert(0, os.path.abspath('../../src'))


# -- Project information -----------------------------------------------------
project = 'Ensemble Analyzer'
copyright = '2026, Andrea Pellegrini, Paolo Righi, Andrea Mazzanti, Michele Mancinelli'
author = r'Andrea Pellegrini,* Paolo Righi, Andrea Mazzanti, Michele Mancinelli'


# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'myst_parser',
    'sphinx_rtd_theme',
    'sphinx_design',
    'sphinx_copybutton',
]

templates_path = ['_templates']
exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

latex_engine = 'lualatex'

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']


# -- MyST Parser configuration -----------------------------------------------
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "colon_fence",
    "substitution"
]
