# File formats

| Format | Open | Save | Notes |
|---|---|---|---|
| Penzene `.penz` | ✓ | ✓ | everything: molecules, arrows, text, shapes, styles |
| MDL MOL / SDF | ✓ | ✓ | V2000 and V3000; a multi-record SDF opens as a grid |
| SMILES `.smi`, InChI `.inchi` | ✓ | | one per line, as a grid; paste either as text too |
| MDL Rxnfile `.rxn` | ✓ | ✓ | and reaction SMILES by copy and paste |
| ChemDraw XML `.cdxml` | ✓ | ✓ | molecules, arrows (curved too), text, lines, boxes, ovals, lone pairs |
| ChemDraw `.cdx` | ✓ | ✓ | as CDXML, on every platform |
| SVG, PNG, PDF | Penzene's own | export | the drawing rides along, so they reopen editable |

## The .penz format

A `.penz` file is JSON: `{"format": "penzene", "version": 1, "atoms": [...], "bonds": [...], ...}`,
with arrows, texts, fills, brackets and drawing settings alongside. It's written to be read back
exactly. The [JSON Schema](_static/penz.schema.json) describes every field. Here is the smallest
drawing, a single bond:

```json
{
    "format": "penzene",
    "version": 1,
    "atoms": [{"x": 0, "y": 0, "z": 6}, {"x": 14.4, "y": 0, "z": 8}],
    "bonds": [{"a": 0, "b": 1, "order": 1}]
}
```

Coordinates are points, x to the right and y down, with bonds 14.4 points long. Atoms are
referred to by their index in `atoms`. A field left out takes its default: carbon, no charge,
a single bond and so on.

### Versions

Every Penzene release reads the files of every earlier one. Examples from each release live in
`tests/data/penz/`, and the tests open all of them.

- **Readers ignore fields they don't know.** A newer Penzene can add a field (say, a new arrow
  setting), and an older one still opens the file and draws what it understands.
- **New fields are optional**, and leaving one out means what files meant before it existed.
- **`version` changes only for a breaking change**: a field that changes meaning, or one that can't
  be left out. A reader refuses a version it doesn't know rather than guess.

## ChemDraw

CDXML import keeps what the file draws: molecules with their labels and abbreviations, arrows, text,
lines, boxes, ovals and lone-pair symbols. Pasting from ChemDraw works on macOS (and on Windows
where ChemDraw offers CDXML). The file's label size relative to its bond length comes with it, so
crowded drawings with small labels look as drawn; choosing a drawing style goes back to that style's
own. Penzene doesn't import ChemDraw's orbitals or TLC plates.
