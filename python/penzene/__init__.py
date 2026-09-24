"""Penzene: ChemDraw-style 2D chemical structures, drawn by the same engine as the app.

    import penzene as pz
    doc = pz.from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")   # caffeine
    doc.export("caffeine.svg")
"""
import os as _os

# A wheel bundles Qt's plugins next to this file; point Qt there before it loads.
_plugins = _os.path.join(_os.path.dirname(__file__), "plugins")
if _os.path.isdir(_plugins):
    _os.environ.setdefault("QT_PLUGIN_PATH", _plugins)

from ._penzene import Atom, Bond, Document, __version__, drawing_styles, from_json, from_smiles, read  # noqa: E402

__all__ = ["Atom", "Bond", "Document", "__version__", "drawing_styles", "from_json", "from_smiles", "read"]
