"""Penzene: 2D chemical structures, drawn by the same engine as the app."""

from collections.abc import Sequence


class Atom:
    """An atom of a Document (read-only; edit through Document)."""

    @property
    def symbol(self) -> str:
        """Element symbol, or the abbreviation (Boc, OMe…) if it has one"""

    @property
    def x(self) -> float:
        """Points, x right"""

    @property
    def y(self) -> float:
        """Points, y down"""

    @property
    def charge(self) -> int:
        """Formal charge"""

    def __repr__(self) -> str: ...

class Bond:
    """A bond of a Document (read-only)."""

    @property
    def a(self) -> int:
        """Index of the first atom (a wedge starts here)"""

    @property
    def b(self) -> int:
        """Index of the second atom"""

    @property
    def order(self) -> int:
        """1, 2 or 3"""

    def __repr__(self) -> str: ...

class Document:
    """A drawing: molecules, arrows and text, as the app holds it."""

    def __init__(self) -> None:
        """An empty drawing."""

    @property
    def atoms(self) -> list[Atom]:
        """The atoms, in order (their index is what add_bond and hotkeys take)"""

    @property
    def bonds(self) -> list[Bond]:
        """The bonds"""

    @property
    def style(self) -> str:
        """Drawing style: 'ACS 1996' (default), 'JDP' or 'RSC'"""

    @style.setter
    def style(self, arg: str, /) -> None: ...

    def add_atom(self, label: str = 'C', x: float = 0.0, y: float = 0.0) -> int:
        """
        Add an atom: an element, an abbreviation (OMe, Boc…) or a SMILES fragment. Returns its index.
        """

    def add_bond(self, a: int, b: int, order: int = 1) -> None:
        """Bond atoms a and b (order 1, 2 or 3)."""

    def hotkeys(self, atom: int, keys: str, bond: int | None = None) -> tuple[int, int]:
        """
        Type hotkeys with this atom as the hotspot (or atom=-1, bond=i for a bond). Returns the final (atom, bond) hotspot.
        """

    def clean(self, atoms: Sequence[int] = []) -> None:
        """
        Lay out afresh with RDKit; with atoms, only the molecules containing them.
        """

    def to_smiles(self) -> str:
        """Canonical SMILES ("" if the drawing isn't valid chemistry)."""

    def to_molblock(self) -> str:
        """MDL MOL (V2000), with the drawing's wedges."""

    def to_inchi(self) -> str:
        """Standard InChI."""

    def to_inchikey(self) -> str:
        """Standard InChIKey."""

    def to_json(self) -> str:
        """The .penz document (JSON)."""

    def to_svg(self) -> str:
        """
        SVG, as the app exports it (the drawing embedded, so it reopens editable).
        """

    def to_png(self, dpi: float = 300) -> bytes:
        """PNG bytes at this resolution."""

    @property
    def formula(self) -> str:
        """Molecular formula, Hill order, all molecules together"""

    @property
    def mw(self) -> float:
        """Molecular weight"""

    @property
    def exact_mass(self) -> float:
        """Monoisotopic mass"""

    def save(self, path: str) -> None:
        """
        Write .penz (full fidelity), ChemDraw .cdxml or .cdx, or MOL for any other extension.
        """

    def export(self, path: str) -> None:
        """Write .svg, .png or .pdf, exactly as the app exports."""

    def __repr__(self) -> str: ...

def from_smiles(smiles: str) -> Document:
    """A structure from SMILES, laid out in 2D."""

def from_json(json: str) -> Document:
    """A Document from .penz JSON (as to_json writes)."""

def read(path: str) -> Document:
    """
    Open a file: .penz, MOL, ChemDraw .cdxml/.cdx, .rxn, a Penzene SVG/PNG, or every record of an SDF, .smi or .inchi file laid out as a grid.
    """

def drawing_styles() -> list[str]:
    """The drawing style names Document.style takes."""
