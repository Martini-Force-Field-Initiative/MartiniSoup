import os
import sys

# Ensure the package is importable when building locally without install
sys.path.insert(0, os.path.abspath("../src"))

# Mock heavy dependencies so autodoc can import modules without a full install
autodoc_mock_imports = [
    "tqdm",
    "scipy",
    "MDAnalysis",
    "freud",
    "pandas",
    "numpy",
    "lmfit",
    "requests",
    "tidynamics",
    "matplotlib",
    "uncertainties",
]

project = "MartiniSoup"
author = "Chris Brasnett"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
]

autodoc_member_order = "bysource"
napoleon_numpy_docstring = True
napoleon_google_docstring = False

html_theme = "furo"
html_title = "MartiniSoup"
