# Python API reference

Generated from the `penzene` module's docstrings. A guide with examples: [Python](python.md).

## Functions

### `from_smiles(smiles: str) -> Document`

A structure from SMILES, laid out in 2D.

### `read(path: str) -> Document`

Open a file: .penz, MOL, ChemDraw .cdxml/.cdx, .rxn, a Penzene SVG/PNG, or every record of an SDF, .smi or .inchi file laid out as a grid.

### `from_json(json: str) -> Document`

A Document from .penz JSON (as to_json writes).

### `drawing_styles() -> list[str]`

The drawing style names Document.style takes.


## Document

A drawing: molecules, arrows and text, as the app holds it.

### `__init__(self) -> None`

An empty drawing.

### `add_atom(self, label: str = 'C', x: float = 0.0, y: float = 0.0) -> int`

Add an atom: an element, an abbreviation (OMe, Boc…) or a SMILES fragment. Returns its index.

### `add_bond(self, a: int, b: int, order: int = 1) -> None`

Bond atoms a and b (order 1, 2 or 3).

### `Document.atoms`

The atoms, in order (their index is what add_bond and hotkeys take)

### `Document.bonds`

The bonds

### `clean(self, atoms: collections.abc.Sequence[int] = []) -> None`

Lay out afresh with RDKit; with atoms, only the molecules containing them.

### `Document.exact_mass`

Monoisotopic mass

### `export(self, path: str) -> None`

Write .svg, .png or .pdf, exactly as the app exports.

### `Document.formula`

Molecular formula, Hill order, all molecules together

### `hotkeys(self, atom: int, keys: str, bond: int | None = None) -> tuple[int, int]`

Type hotkeys with this atom as the hotspot (or atom=-1, bond=i for a bond). Returns the final (atom, bond) hotspot.

### `Document.mw`

Molecular weight

### `save(self, path: str) -> None`

Write .penz (full fidelity), ChemDraw .cdxml or .cdx, or MOL for any other extension.

### `Document.style`

Drawing style: 'ACS 1996' (default), 'JDP' or 'RSC'

### `to_inchi(self) -> str`

Standard InChI.

### `to_inchikey(self) -> str`

Standard InChIKey.

### `to_json(self) -> str`

The .penz document (JSON).

### `to_molblock(self) -> str`

MDL MOL (V2000), with the drawing's wedges.

### `to_png(self, dpi: float = 300) -> bytes`

PNG bytes at this resolution.

### `to_smiles(self) -> str`

Canonical SMILES ("" if the drawing isn't valid chemistry).

### `to_svg(self) -> str`

SVG, as the app exports it (the drawing embedded, so it reopens editable).


## Atom

An atom of a Document (read-only; edit through Document).

### `Atom.__init__`



### `Atom.charge`

Formal charge

### `Atom.symbol`

Element symbol, or the abbreviation (Boc, OMe…) if it has one

### `Atom.x`

Points, x right

### `Atom.y`

Points, y down


## Bond

A bond of a Document (read-only).

### `Bond.__init__`



### `Bond.a`

Index of the first atom (a wedge starts here)

### `Bond.b`

Index of the second atom

### `Bond.order`

1, 2 or 3

