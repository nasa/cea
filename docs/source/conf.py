# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

from docutils import nodes
from sphinx.directives import ObjectDescription

# Avoid initializing the native library during docs builds.
os.environ.setdefault("CEA_SKIP_INIT", "1")

sys.path.insert(0, os.path.abspath("../../source/bind/python/"))
import cea

project = 'CEA'
copyright = ''
author = 'Mark Leader'
version = '3.0'
release = '3.0.4'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "breathe",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.githubpages",
    "sphinx.ext.mathjax",
    # "sphinx.ext.autosectionlabel",
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

breathe_projects = {
    'cea': '../doxygen/xml',
}
breathe_default_project = 'cea'
breathe_domain_by_extension = {
    'py': 'py',
    'cs': 'cs',
    'f': 'fortran',
    'F': 'fortran',
    'for': 'fortran',
    'F90': 'fortran',
    'f90': 'fortran',
}

# Breathe doesn't register Fortran ``type`` compounds by default; map them to
# the standard C++ class directive so derived types render without errors.
try:
    from breathe.renderer.sphinxrenderer import (
        BaseObject,
        DomainDirectiveFactory,
        SphinxRenderer,
    )

    DomainDirectiveFactory.cpp_classes.setdefault(
        "type", DomainDirectiveFactory.cpp_classes["class"]
    )

    class FortranObject(BaseObject, ObjectDescription):
        """Minimal directive that renders Fortran signatures verbatim."""

        has_content = True

        def handle_signature(self, sig, signode):
            text = " ".join(sig.split())
            signode += nodes.Text(text)
            return text

        def needs_arglist(self):
            return False

    _FORTRAN_FILE_SUFFIXES = (".f", ".for", ".f90", ".f95", ".f03", ".f08")

    DomainDirectiveFactory.fortran_classes = {
        "function": (FortranObject, "function"),
        "variable": (FortranObject, "variable"),
        "type": (FortranObject, "type"),
        "interface": (FortranObject, "type"),
        "class": (FortranObject, "type"),
        "struct": (FortranObject, "type"),
        "namespace": (FortranObject, "module"),
    }

    _original_domain_directive_create = DomainDirectiveFactory.create

    def _fortran_domain_directive(domain, args, _orig=_original_domain_directive_create):
        if domain == "fortran":
            kind = args[0]
            mapping = getattr(DomainDirectiveFactory, "fortran_classes", {})
            if kind in mapping:
                cls, name = mapping[kind]
                directive_args = [f"{domain}:{name}"] + args[1:]
                return cls(*directive_args)
        return _orig(domain, args)

    DomainDirectiveFactory.create = staticmethod(_fortran_domain_directive)

    _original_get_domain = SphinxRenderer.get_domain

    def _is_fortran_context(context):
        node_stack = getattr(context, "node_stack", None)
        if not node_stack:
            return False
        node = node_stack[0]
        language = getattr(node, "language", None)
        if isinstance(language, str) and language.lower().startswith("fortran"):
            return True
        location = getattr(node, "location", None)
        if location:
            source = getattr(location, "file", "") or ""
            _, ext = os.path.splitext(source)
            if ext.lower() in _FORTRAN_FILE_SUFFIXES:
                return True
        return False

    def _fortran_get_domain(self):
        context = getattr(self, "context", None)
        if context and _is_fortran_context(context):
            return "fortran"
        domain = _original_get_domain(self)
        return domain

    SphinxRenderer.get_domain = _fortran_get_domain

    _original_handle_declaration = SphinxRenderer.handle_declaration

    def _fortran_safe_handle_declaration(self, *args, **kwargs):
        if self.get_domain() == "fortran" and kwargs.get("display_obj_type"):
            kwargs = dict(kwargs)
            kwargs["display_obj_type"] = None
        return _original_handle_declaration(self, *args, **kwargs)

    SphinxRenderer.handle_declaration = _fortran_safe_handle_declaration
except Exception:
    pass

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = "python"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

# -----------------------------------------------------------------------------
# Document both class docstring and init docstring
autoclass_content = "both"

# Order members by source and not alphabetical order
autodoc_member_order = "bysource"

# html_theme_options = {
#     "prev_next_buttons_location": "bottom",  # or "top", "both", "sidebar", "None"
# }
