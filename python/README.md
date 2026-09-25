# penzene

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

Also: `pz.read(path)` for `.penz`, `.mol` and ChemDraw `.cdxml`; `to_smiles()`, `to_molblock()`, `to_inchi()`, `to_inchikey()`, `to_svg()`, `to_png(dpi)`; `doc.style = "RSC"` (or `"ACS 1996"`, `"JDP"`); abbreviations such as `add_atom("Boc")`.

Penzene brings its own copy of Qt. Loading it into a process that already has PySide6 or PyQt (for example Jupyter with a Qt event loop) may conflict; the plain Jupyter kernel is fine.

GPL-3.0.
