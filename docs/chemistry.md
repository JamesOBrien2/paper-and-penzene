# Chemistry

Penzene keeps the drawing and hands its chemistry to RDKit. Everything below updates as you draw.

- **Formula and mass.** The status bar shows the formula, molecular weight and exact mass of the
  selection, or of the whole drawing. **View → Properties Panel** (Ctrl+I) adds elemental analysis,
  cLogP, TPSA, H-bond donors and acceptors, rotatable bonds, heavy atoms, and Lipinski and Veber
  checks, with **Copy as Text** for a supporting-information table.
- **Implicit hydrogens and valence.** Labels get their hydrogens (OH, NH₂). An atom with too many
  bonds is drawn in red.
- **Stereochemistry.** Wedges and hashes set stereocentres; E/Z comes from the drawing.
  **View → Show Stereo Labels** shows CIP (R)/(S) and (E)/(Z).
- **Check Structure** (Structure menu) lists valence errors, unknown labels, overlapping atoms,
  stereocentres without a wedge, and wedges on atoms that aren't stereocentres. Click one to select it.
- **Clean** (Ctrl+Shift+K) lays a structure out afresh, keeping its stereo and bond styles.
- **Hydrogens.** Structure → Add or Remove Explicit Hydrogens. View → Carbon Labels and Show Implicit
  Hydrogens change how they're shown.
- **Aromatic circles** (View menu) draw benzene-like rings with a circle, for the whole drawing or
  for selected rings.
- **Numbers.** View → Atom Numbers shows atom indices; `'` on an atom sets its reaction atom-map
  number, which goes into SMILES and MOL.
- **Names (online).** File → Import Name from PubChem turns a name into a structure, and
  Structure → Name from PubChem copies a structure's IUPAC name. Both look the compound up on
  PubChem, so they need an internet connection and only work for compounds PubChem knows.
- **Copy as** SMILES, InChI, InChIKey and reaction SMILES (Edit menu).
- **What isn't chemistry.** Interaction and partial bonds, lone pairs, δ, brackets, arrows, shapes
  and text are drawn only. Free-text labels are generic atoms (`*` in SMILES).
