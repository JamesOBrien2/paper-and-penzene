# Sphinx configuration for penzene.readthedocs.io. Pages are Markdown (MyST); the Furo theme
# carries the app's lab notebook colours.
project = "Penzene"
author = "James O'Brien"
copyright = "James O'Brien"

extensions = ["myst_parser", "sphinx_copybutton"]
source_suffix = {".md": "markdown"}
exclude_patterns = ["_build", "shortcuts.md"]
myst_heading_anchors = 3
myst_enable_extensions = ["colon_fence"]

html_theme = "furo"
html_title = "Penzene"
html_logo = "_static/logo.svg"
html_favicon = "_static/logo.svg"
html_static_path = ["_static"]
html_css_files = ["penzene.css"]
html_theme_options = {
    "source_repository": "https://github.com/JamesOBrien2/penzene/",
    "source_branch": "main",
    "source_directory": "docs/",
    "light_css_variables": {
        "color-brand-primary": "#0F6E56",
        "color-brand-content": "#0F6E56",
        "color-brand-visited": "#0F6E56",
        "color-background-primary": "#FBF8F1",
        "color-background-secondary": "#F3EFE4",
        "color-background-hover": "#E1F5EE",
        "color-background-border": "#E4E1D6",
        "color-foreground-primary": "#2C2C2A",
        "color-foreground-secondary": "#5F5E5A",
        "color-sidebar-background": "#F3EFE4",
        "color-sidebar-item-background--hover": "#E1F5EE",
        "color-code-background": "#FFFFFF",
        "color-admonition-background": "#FFFFFF",
        "color-admonition-title--note": "#0F6E56",
        "color-admonition-title-background--note": "#E1F5EE",
    },
    "dark_css_variables": {
        "color-brand-primary": "#5DCAA5",
        "color-brand-content": "#5DCAA5",
        "color-brand-visited": "#5DCAA5",
        "color-background-primary": "#22211F",
        "color-background-secondary": "#1B1A18",
        "color-background-hover": "#0B3B30",
        "color-background-border": "#444441",
        "color-foreground-primary": "#F1EFE8",
        "color-foreground-secondary": "#B4B2A9",
        "color-sidebar-background": "#1B1A18",
        "color-sidebar-item-background--hover": "#0B3B30",
        "color-code-background": "#2C2C2A",
        "color-admonition-background": "#2C2C2A",
        "color-admonition-title--note": "#5DCAA5",
        "color-admonition-title-background--note": "#0B3B30",
    },
}
copybutton_prompt_text = r"\$ |>>> "
copybutton_prompt_is_regexp = True
