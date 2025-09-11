# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))

project = 'QuLTRA'
copyright = '2025, Simona Zaccaria, Antonio Gnudi'
author = 'Simona Zaccaria, Antonio Gnudi'
release = 'v0.1.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'nbsphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
]


templates_path = ['_templates']
exclude_patterns = []

# Aggiungi questa configurazione per autodoc
autodoc_default_options = {
    'members': True,
    'show-inheritance': True
}

napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_examples = True
napoleon_use_param = True
napoleon_use_rtype = True


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
