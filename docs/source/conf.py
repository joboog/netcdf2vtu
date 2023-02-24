# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('../../src/netcdf2vtu'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'netcdf2vtu'
copyright = '2023, Helmholtz Centre for Environmental Research - UFZ (www.ufz.de)'
author = 'Johannes Boog'
release = '0.1'
version = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx_rtd_theme',
              'myst_parser'
             # 'sphinx.ext.autosummary',
            #'sphinx_autodoc_typehints'
]

# autosummaries from source-files
#autosummary_generate = True

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
