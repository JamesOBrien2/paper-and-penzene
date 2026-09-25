# File formats

| Format | Open | Save | Notes |
|---|---|---|---|
| Penzene `.penz` | ✓ | ✓ | everything: molecules, arrows, text, shapes, styles |
| MDL MOL / SDF | ✓ | ✓ | V2000 and V3000; a multi-record SDF opens as a grid |
| SMILES `.smi`, InChI `.inchi` | ✓ | | one per line, as a grid; paste either as text too |
| MDL Rxnfile `.rxn` | ✓ | ✓ | and reaction SMILES by copy and paste |
| ChemDraw XML `.cdxml` | ✓ | ✓ | molecules, arrows (curved too), text, lines, boxes, ovals, lone pairs |
| ChemDraw `.cdx` | ✓ | ✓ | molecules only, and not in the Windows build yet |
| SVG, PNG, PDF | Penzene's own SVG and PNG | export | |

## The .penz format

A `.penz` file is JSON: `{"format": "penzene", "version": 1, "atoms": [...], "bonds": [...], ...}`,
with arrows, texts, fills, brackets and drawing settings alongside. It's written to be read back
exactly. Its full, versioned description comes with Penzene 1.0.

## ChemDraw

CDXML import keeps what the file draws: molecules with their labels and abbreviations, arrows, text,
lines, boxes, ovals and lone-pair symbols. Pasting from ChemDraw works on macOS (and on Windows
where ChemDraw offers CDXML). Penzene doesn't import ChemDraw's orbitals or TLC plates.
