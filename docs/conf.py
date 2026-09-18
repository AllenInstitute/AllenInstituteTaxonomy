"""Sphinx configuration for the Allen Institute Taxonomy documentation."""

project = "Allen Institute Taxonomy"
copyright = "2026, Allen Institute"
author = "Allen Institute"

extensions = ["myst_parser"]
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
suppress_warnings = ["myst.header"]

html_theme = "sphinx_rtd_theme"
html_static_path = []

myst_heading_anchors = 3
