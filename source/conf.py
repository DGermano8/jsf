# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath('../jsf/'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Jump-Switch-Flow'
copyright = '2023, Domenic Germano'
author = 'Domenic Germano'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
# html_theme = 'classic'
html_theme = 'bizstyle'
# html_theme = 'sphinxdoc'
# html_theme = 'nature'
# html_theme = 'pyramid'
# html_theme = 'scrolls'
# html_theme = 'agogo'
# html_theme = 'traditional'
# html_theme = 'haiku'

html_static_path = ['_static']
html_css_files = ['custom.css',]
