"""The Python module: the #47 snippet as written, plus the rest of the API."""
import os
import tempfile

import penzene as pz

CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
out = tempfile.mkdtemp()

doc = pz.from_smiles(CAFFEINE)
doc.clean()
doc.save(os.path.join(out, "caffeine.penz"))
doc.export(os.path.join(out, "caffeine.svg"))
assert "<svg" in doc._repr_svg_()
assert "<path" in doc.to_svg(), "labels must render as glyph outlines"
assert doc.to_png(dpi=72)[:8] == b"\x89PNG\r\n\x1a\n"
assert doc.formula == "C8H10N4O2" and abs(doc.mw - 194.19) < 0.01
assert doc.to_inchikey() == "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
assert len(doc.atoms) == 14 and doc.atoms[0].symbol == "C"

again = pz.read(os.path.join(out, "caffeine.penz"))
assert again.to_smiles() == doc.to_smiles()
assert pz.from_json(doc.to_json()).to_smiles() == doc.to_smiles()
assert pz.read(os.path.join(os.environ["PENZENE_TEST_DATA"], "caffeine.mol")).formula == "C8H10N4O2"

# The hotkey builder: ChemDraw's cheat-sheet dipeptide, 42n152o from H2N-CH3.
m = pz.Document()
a = m.add_atom("N")
b = m.add_atom("C", x=14.4)
m.add_bond(a, b)
m.hotkeys(b, "42n152o")
assert m.formula == "C6H12N2O3", m.formula  # Ala-Ala

m2 = pz.Document()
a = m2.add_atom("N")
m2.hotkeys(a, "42n152o")  # the issue's snippet, from a lone N: must just work
assert m2.to_smiles()

boc = pz.Document()
boc.add_atom("Boc")
assert boc.atoms[0].symbol == "Boc" and boc.formula == "C5H10O2"

doc.style = "RSC"
assert doc.style == "RSC" and "RSC" in pz.drawing_styles()
for bad in (lambda: pz.from_smiles("C1CC"), lambda: m.add_atom("notachem!!"), lambda: setattr(doc, "style", "Comic")):
    try:
        bad()
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")
print("python ok", pz.__version__)
