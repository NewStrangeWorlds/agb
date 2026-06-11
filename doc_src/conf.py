# Configuration file for the Sphinx documentation builder.
#
# Full documentation of the AGB stationary dust-driven wind model.
# Build with:  make html   (output goes to ../doc_built/html)

# -- Project information -----------------------------------------------------

project = "AGB Wind Model"
copyright = "Daniel Kitzmann, Joachim Stock"
author = "Daniel Kitzmann, Joachim Stock"
release = "dev"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.mathjax",      # render LaTeX maths in HTML
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosectionlabel",
]

# Allow :ref: to a section by its title, namespaced by document to avoid clashes.
autosectionlabel_prefix_document = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# A short numbering of the displayed equations.
math_number_all = False
math_eqref_format = "Eq. {number}"

todo_include_todos = True

# -- Options for HTML output -------------------------------------------------

html_theme = "furo"
html_theme_options = {
    "navigation_with_keys": True,
}
html_static_path = ["_static"]
html_title = "AGB Dust-Driven Wind Model"
html_logo = '_static/agb-icon-512.png'

# -- Options for LaTeX / PDF output ------------------------------------------

latex_elements = {
    "papersize": "a4paper",
    "pointsize": "11pt",
}
