# Drawing

## Tools and the hotspot

The palette on the left holds the tools: select, eraser, bonds (single, double, triple, wedge,
hash, interaction, partial), chain, rings, atoms, charges, colour, arrows, lines and shapes, and
text. Hover over a tool to see its key.

Point at an atom or bond and it becomes the **hotspot**, marked in green. It stays put when the mouse
moves off, so you can keep typing. Keys act on the hotspot: `1` sprouts a bond, `2`
a carbonyl, `a` a phenyl, `O` an OMe, and `+`/`−` change the charge. The arrow keys walk the
molecule atom → bond → atom. Every key is listed in [Keyboard shortcuts](keys.md) and under
**Help → Keyboard Shortcuts** (F1).

With no hotspot (press Esc), keys pick tools: `x` bond, `X` chain, `j` benzene, `e` arrow, `t` text,
Space select.

## Atoms and labels

- Press **Enter**, `=` or `t` on an atom (or click it with the Text tool) to type a label: an element
  (`Br`, `NH2`), an abbreviation (`OMe`, `Boc`, `TBS`, `CO2Me`, …), a SMILES fragment (drawn out), or
  any text (`R1`, `X`, `MgEt`). Text that isn't chemistry is drawn as written and counts as a generic atom.
- `x` and `r` on an atom label it X and R. `R1`, `R2`… export as MDL R-groups.
- **Structure → Expand Abbreviations** draws them out in full.
- Atom Properties (`/` on an atom) sets the charge, map number, lone pairs, radical electrons and δ±.
- `:` cycles lone pairs, `*` toggles a radical dot, and the atom's context menu adds δ+ or δ−.
  Radicals are chemistry (one H fewer each); lone pairs and δ are only drawn.
- `.` adds an attachment point. `j` and `J` bond an η⁵-cyclopentadienyl or η⁶-benzene through the
  ring's centre.

## Bonds

- `1 2 3` set the order; `2` on a double bond moves its second line to the other side, and `l c r`
  place it left, centred or right.
- `w` / `h` wedge and hash (again to flip); `b` bold, `d` dashed, `y` wavy.
- `i` makes a bond an **interaction** (dotted: hydrogen bonds, contacts, coordination), and `p` / `P`
  a **partial** bond (dashed single or solid + dashed: bonds forming or breaking in a transition
  state). Neither counts as a covalent bond, so a drawn transition state keeps its reactants'
  formula.
- `f` brings a bond to the front: bonds it crosses are drawn with a gap.

## Rings and templates

Ring tools draw 3- to 8-membered rings and benzene: click empty space, an atom (spiro or attached),
or a bond (fused). On a bond, `a z v 4–8 9 0` fuse rings, `9` and `0` chair cyclohexanes.

**View → Templates** opens the template library: amino acids, sugars, nucleobases, common scaffolds,
and Haworth, Fischer and Newman projections. Haworth and Fischer drawings give the right stereo
when converted to SMILES. **Structure → Save Selection as Template** adds your own.

## Selecting and arranging

- Drag to select; double-click selects a whole molecule. Drag a selection to move it, and drop an
  atom on another to merge them.
- The arrow keys nudge a selection (Shift for 10 points). Ctrl+arrow duplicates it across the next
  reaction arrow.
- **Rotate:** drag the round handle above a selection (Shift for 15° steps, Ctrl for 45°), Alt+drag,
  or Alt+←/→ in 15° steps. **Out of the page:** Shift+Alt+drag or
  Shift+Alt+arrows turn it in 3D and keep its stereo; **Structure → Turn Over** flips it 180°.
- **Stretch and squash:** drag the handles around a selection (corners scale, edges stretch), or use
  **Structure → Transform** for exact values.
- **Structure → Align and Distribute**, **Flip**, and **Arrange Scheme** (lines up a reaction scheme,
  with agents centred over their arrows).
- **Structure → Brackets** puts square or round brackets, with a subscript such as *n*, around the
  selected atoms.

## Arrows, text, shapes and colour

- Arrows: reaction, equilibrium, resonance, retrosynthesis, and curved (electron pair) and fishhook
  (single electron) arrows. Click a curved arrow again to flip its curve.
- Text: formulas get subscripts automatically (H2O → H₂O).
- Lines (solid or dashed), boxes, rounded boxes and ellipses for grouping; Shift draws a square or
  circle.
- The colour tool paints atoms, bonds, arrows and text; the fill tool shades a ring.
