# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "twirl"
copyright = "2023, Lionel Garcia"
author = "Lionel Garcia"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
extensions = [
    "myst_nb",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
]

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]

source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}

root_doc = "index"

html_title = "twirl"
html_theme_options = {
    "repository_url": "https://github.com/lgrcia/twirl",
    "use_repository_button": True,
}

nb_render_image_options = {"align": "center"}

html_css_files = ["style.css"]
autodoc_typehints = "signature"
autoclass_content = "both"
