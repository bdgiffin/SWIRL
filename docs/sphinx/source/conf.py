# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SWIRL'
copyright = '2024, Brian Doran Giffin'
author = 'Brian Doran Giffin'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Further details on how to include a bibliography may be found here:
#   https://chiplicity.readthedocs.io/en/latest/Using_Sphinx/UsingBibTeXCitationsInSphinx.html
extensions = ['sphinxcontrib.bibtex']
bibtex_bibfiles = ['_static/references.bib']
bibtex_default_style = 'unsrt'

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
