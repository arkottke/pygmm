# Configuration file for the Sphinx documentation builder.
import os
import sys
from datetime import date

# Add the src directory to the path
sys.path.insert(0, os.path.abspath("../src"))

import pygmm

project = "pyGMM"
copyright = f"2024-{date.today().year}, Albert Kottke"
author = "Albert Kottke"
release = pygmm.__version__

# -- General configuration ---------------------------------------------------

extensions = [
    "matplotlib.sphinxext.plot_directive",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinxcontrib.bibtex",
    "sphinxext.opengraph",
    "sphinx_copybutton",
    "sphinx_design",
    "myst_parser",
    "nbsphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Bibliography configuration
bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"

# Autosummary configuration
autosummary_generate = True

# Napoleon configuration for NumPy-style docstrings
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# Intersphinx configuration
intersphinx_mappings = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

# MyST configuration
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "deflist",
    "fieldlist",
    "html_admonition",
    "html_image",
    "colon_fence",
    "smartquotes",
    "replacements",
    "linkify",
    "strikethrough",
    "substitution",
    "tasklist",
]

# MyST parsing configuration
myst_heading_anchors = 3

# -- Options for HTML output -------------------------------------------------
pygments_style = "friendly"

html_theme = "furo"
html_static_path = ["_static"]
html_title = f"{project} v{release}"
html_css_files = [
    "custom.css",
]

# Furo theme options
html_theme_options = {
    "sidebar_hide_name": True,
    "navigation_with_keys": True,
    "top_of_page_button": "edit",
    "source_repository": "https://github.com/arkottke/pygmm/",
    "source_branch": "main",
    "source_directory": "docs/",
}

# Show typehints as content of the function or method
autodoc_typehints = "description"
autodoc_member_order = "bysource"
autodoc_default_options = {
    "members": True,
    "inherited-members": True,
    "private-members": False,
    "special-members": "__init__",
    "undoc-members": True,
    "exclude-members": "__weakref__",
    "no-index": False,
}

# Notebook execution configuration
nbsphinx_execute = "never"  # Don't execute notebooks during build
nbsphinx_allow_errors = True
