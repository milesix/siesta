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


# The master toctree document.
master_doc = 'index'

# -- Project information -----------------------------------------------------


# Add any paths that contain templates here, relative to this directory.
#templates_path = ['_templates']

project = 'Siesta Documentation'
full_title = project
copyright = '2023, The Siesta Group'
author = 'The Siesta Group'
version = "0.2"
release = version


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
]

#
# -- Options for autosectionlabel ---------------------------------------------
#

autosectionlabel_prefix_document = True
autosectionlabel_maxdepth = 2

#
# -- Options for intersphinx --------------------------------------------------
#

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

#
# -- Options for the theme ----------------------------------------------------
#

html_theme = 'sphinx_rtd_theme'

extlinks = {
    'doi': ('https://doi.org/%s', 'doi:'),
}

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', "env"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTML output -------------------------------------------------

html_title = full_title
htmlhelp_basename = "SiestaDocs"


# Enable labeling for figures
numfig = True

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']

#
# -- Options for latex output -------------------------------------------------
#

latex_documents = [
    (master_doc, htmlhelp_basename + ".tex", full_title, author, "manual"),
]

#
# -- Options for manual page output -------------------------------------------
#

man_pages = [
    (master_doc, htmlhelp_basename, full_title, [author], 1)
]

#
# -- Options for TODOs --------------------------------------------------------
#

todo_include_todos = True
todo_link_only = True
todo_emit_warnings = False
