
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import sphinx_rtd_theme
import sys
import os
sys.path.insert(0, os.path.abspath('../../translatome/src/'))

# pip install sphinx
# pip install myst_parser
# pip install sphinx_rtd_theme
# pip install sphinxcontrib.bibtex
#
# to generate doc for python files:
# cd docs/
# sphinx-apidoc -f -o source/generated/ ../translatome/src/

# -- Project information -----------------------------------------------------

project = 'TRAIN'
copyright = '2023, Julie Ripoll, Céline Mandier, Fati Chen, Eric Rivals'
author = 'Julie Ripoll'

# The full version, including alpha/beta/rc tags
release = '0.1.0'


# -- General configuration ---------------------------------------------------

# file extensions of source files. Sphinx considers the files with this suffix
# as sources.
# The value can be a dictionary mapping file extensions to file types.
source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'restructuredtext',
    '.md': 'markdown',
}

# Files with a suffix that is not in the dictionary will be parsed with the
# default reStructuredText parser.
#source_parsers = {'.md': 'recommonmark.parser.CommonMarkParser'}

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['myst_parser',
              'sphinx.ext.duration',
              'sphinxcontrib.bibtex',
              'sphinx.ext.doctest',
              'sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx.ext.napoleon']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for bibtex bibliography------------------------------------------

# define bib file
bibtex_bibfiles = ['../../report/src/pol_seq_references.bib']

# encoding, default: utf-8-sig, or latin
bibtex_encoding = 'utf-8'

# style supported: plain, unsrt and unsrtalpha
# import _templates.template_style as template
bibtex_default_style = 'unsrtalpha'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#


html_theme = 'sphinx_rtd_theme'  # sphinx_rtd_theme or bizstyle

# options for the html theme
html_theme_options = {
    'canonical_url': '',
    'analytics_id': '',
    'logo_only': False,
    'display_version': True,  # version number is shown at the top of the sidebar
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    # Toc options
    'includehidden': True,
    # With this enabled, navigation entries are not expandable
    'collapse_navigation': True,
    # Scroll the navigation with the main page content as you scroll the page
    'sticky_navigation': True,
    'navigation_depth': -1,  # Set this to -1 to allow unlimited depth
    'titles_only': False
}

html_sidebars = {
    '**': [
        'globaltoc.html',
    ]
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
