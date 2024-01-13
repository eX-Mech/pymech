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
import os
import sys
from datetime import date
from importlib.metadata import metadata

sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------

project = "pymech"
_meta = metadata(project)
_today = date.today()

author = _meta["Author"]
copyright = f"2016-{_today.year}"  # ", {author}"
master_doc = "index"

# The full version, including alpha/beta/rc tags
version = _meta.get("Version")
release = ".".join(version.split(".")[:3])

_py_min_version = _meta.get("Requires-Python").split(">=")[-1]

rst_prolog = f"""
.. |author| replace:: {author}
.. |today| replace:: {_today}
.. |py_min_version| replace:: {_py_min_version}
"""

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be

# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    # for ipython3 Pygments lexer
    # https://nbsphinx.readthedocs.io/en/0.9.1/installation.html#Pygments-Lexer-for-Syntax-Highlighting # noqa
    "IPython.sphinxext.ipython_console_highlighting",
    "myst_nb",
    "sphinx.ext.autodoc",
    # 'sphinx.ext.autosummary',
    # "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    #  "sphinx.ext.coverage",
    #  "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",  # Numpy-style docstrings
    #  "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinx_inline_tabs",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
    ".myst.md": "myst-nb",
}
nb_execution_mode = "cache"
nb_execution_in_temp = True
nb_execution_raise_on_error = True
nb_execution_show_tb = True
nb_execution_timeout = 600
nb_merge_streams = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "asv_bench/.asv"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

# -- Other options -----------------------------------------------------------
#  autosummary_generate = True

autodoc_default_options = {
    "members": True,
}

autodoc_mock_imports = ["tvtk", "pymech._version"]

# -- Options for Intersphinx -------------------------------------------------

intersphinx_mapping = dict(
    python=("https://docs.python.org/3", None),
    nek=("https://nek5000.github.io/NekDoc", None),
    xr=("https://xarray.pydata.org/en/stable/", None),
)

# -- MyST options ------------------------------------------------------------

myst_heading_anchors = 2
myst_enable_extensions = ["amsmath", "dollarmath", "colon_fence"]


def asv_publish():
    import sys
    from unittest.mock import patch

    cur_dir = os.getcwd()
    old_argv = sys.argv

    try:
        from asv.main import main

        os.chdir("asv_bench")
        sys.argv = ("asv", "publish")

        with patch("sys.exit"):
            main()
    finally:
        sys.argv = old_argv
        os.chdir(cur_dir)


asv_publish()
