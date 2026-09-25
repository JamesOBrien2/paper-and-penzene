# Python

```sh
pip install penzene
```

Publication-quality 2D chemical structures from Python, drawn by the same engine as the [Penzene](https://github.com/JamesOBrien2/penzene) desktop app.

```python
import penzene as pz

doc = pz.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
doc.clean()
doc.export("aspirin.svg")        # identical to the app's export (also .png, .pdf)
doc.save("aspirin.penz")         # open and keep editing in the app
doc.formula, doc.mw               # ('C9H8O4', 180.16…)
doc                               # renders inline in Jupyter

# Build with hotkeys: from H2N-CH3, "42n152o" makes Ala-Ala
m = pz.Document()
n = m.add_atom("N"); c = m.add_atom("C", x=14.4); m.add_bond(n, c)
m.hotkeys(c, "42n152o")
```

Also: `pz.read(path)` for `.penz`, MOL, ChemDraw `.cdxml`/`.cdx` and more; `to_smiles()`, `to_molblock()`, `to_inchi()`, `to_inchikey()`, `to_svg()`, `to_png(dpi)`; `doc.style = "RSC"` (or `"ACS 1996"`, `"JDP"`); abbreviations such as `add_atom("Boc")`.

Penzene brings its own copy of Qt. Loading it into a process that already has PySide6 or PyQt (for example Jupyter with a Qt event loop) may conflict; the plain Jupyter kernel is fine.

Every function and property: [Python API reference](python-api.md).

## Stability

From 1.0 the module follows [semantic versioning](https://semver.org):

- **The public API is what `penzene.__all__` lists**: `Document`, `Atom`, `Bond`, `read`,
  `from_smiles`, `from_json`, `drawing_styles` and `__version__`, with their methods and properties.
  Names that start with an underscore, such as `penzene._penzene`, can change in any release.
- **Minor releases (1.1, 1.2…) only add.** New functions, methods and optional arguments; existing
  calls keep working and return the same kinds of values.
- **Removing or changing something takes a deprecation first.** It keeps working for at least one
  minor release while raising a `DeprecationWarning` that names the replacement, and goes in the
  next major release. The [changelog](https://github.com/JamesOBrien2/penzene/blob/main/CHANGELOG.md)
  lists each deprecation.
- **Drawings may look slightly different between releases** as the layout and rendering improve.
  The chemistry (SMILES, InChI, formula) of a document doesn't change.

The package ships type stubs (`py.typed`), so editors and type checkers see every signature and
docstring.
