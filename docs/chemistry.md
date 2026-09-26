# Chemistry

Penzene keeps the drawing and hands its chemistry to RDKit. Everything below updates as you draw.

```{image} _static/properties-light.png
:alt: The Properties panel for aspirin: formula, masses, cLogP, TPSA and the Lipinski and Veber checks
:class: shot only-light
```

```{image} _static/properties-dark.png
:alt: The Properties panel for aspirin: formula, masses, cLogP, TPSA and the Lipinski and Veber checks, dark theme
:class: shot only-dark
```

::::{container} features-grid

```{feature} scale
:title: Formula and mass
The status bar shows the formula, molecular weight and exact mass of the selection, or of the whole
drawing. **View → Properties Panel** (Ctrl+I) adds elemental analysis, cLogP, TPSA, H-bond donors
and acceptors, rotatable bonds, heavy atoms, and Lipinski and Veber checks, with **Copy as Text**
for a supporting-information table. **File → Export Descriptors…** writes them as a CSV, one row
per molecule ([descriptor tables](cli.md#descriptor-tables)).
```

```{feature} atom
:title: Implicit hydrogens and valence
Labels get their hydrogens (OH, NH₂). An atom with too many bonds is drawn in red.
```

```{feature} rotate-3d
:title: Stereochemistry
Wedges and hashes set stereocentres; E/Z comes from the drawing. **View → Show Stereo Labels**
shows CIP (R)/(S) and (E)/(Z).
```

```{feature} checklist
:title: Check Structure
The Structure menu lists valence errors, unknown labels, overlapping atoms, stereocentres without a
wedge, and wedges on atoms that aren't stereocentres. Click one to select it.
```

```{feature} wand
:title: Clean
Ctrl+Shift+K lays a structure out afresh, keeping its stereo and bond styles.
```

```{feature} letter-h
:title: Hydrogens
**Structure → Add or Remove Explicit Hydrogens**. **View → Carbon Labels** and **Show Implicit
Hydrogens** change how they're shown.
```

```{feature} hexagon
:title: Aromatic circles
The View menu draws benzene-like rings with a circle, for the whole drawing or for selected rings.
```

```{feature} hash
:title: Numbers
**View → Atom Numbers** shows atom indices; `'` on an atom sets its reaction atom-map number, which
goes into SMILES and MOL.
```

```{feature} world-search
:title: Names (online)
**File → Import Name from PubChem** turns a name into a structure, and **Structure → Name from
PubChem** copies a structure's IUPAC name. Both look the compound up on PubChem, so they need an
internet connection and only work for compounds PubChem knows.
```

```{feature} copy
:title: Copy as
SMILES, InChI, InChIKey and reaction SMILES (Edit menu).
```

::::

:::{note}
**What isn't chemistry.** Interaction and partial bonds, lone pairs, δ, brackets, arrows, shapes and
text are drawn only. Free-text labels are generic atoms (`*` in SMILES).
:::
