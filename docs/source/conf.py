# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information



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

# Abilita supporto JavaScript (importante per GitHub Pages)
html_js_files = [
    'https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.0/jquery.min.js'
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
