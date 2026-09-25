# Keyboard shortcuts compared with ChemDraw

Penzene follows ChemDraw's hotkeys. This table checks every entry in ChemDraw 22's
`hotkeys.xml` (the file ChemDraw itself reads), plus the keys on its *Hotkeys
Cheat Sheet*, against Penzene. The full list of Penzene's keys is under
Help → Keyboard Shortcuts (F1).

**Status:** ✓ same · ≈ same idea, different detail · ✗ not in Penzene

## Atom hotspot

| Key | ChemDraw | Penzene | |
|---|---|---|---|
| c n/w o/q s p f l/C b i h | C N O S P F Cl Br I H | same | ✓ |
| B S L | B, Si, Li | same | ✓ |
| m e P A | Me, Et, Ph, Ac | same | ✓ |
| O N F E Z | OMe, NO₂, CF₃, CO₂Me, N₃ | same | ✓ |
| y, H, Q, M | Boc, Cbz, Fmoc, MgBr | same (Y also gives Boc) | ✓ |
| x r | X, R | same (generic atoms, #156) | ✓ |
| d | D (deuterium) | H: Penzene has no isotopes yet | ≈ |
| ! | query atom `?` | none (no query atoms) | ✗ |
| 0 1 | bond, cyclic / linear mode | same | ✓ |
| 2 | carbonyl | same | ✓ |
| 3 a | benzene | same | ✓ |
| 4 5 | wedged / hashed bond | same | ✓ |
| 6 7 u v | cyclohexane, cyclopentane, cyclobutane, cyclopropane | same | ✓ |
| 8 9 z | methylidene, dimethyl, alkyne | same | ✓ |
| k K | sulfonyl, t-Bu | same | ✓ |
| . | attachment point | same (wavy bond to a bare point) | ✓ |
| j J | η⁵-cyclopentadienyl, η⁶-benzene | same, flat ring bonded through its centre (ChemDraw draws it in perspective) | ≈ |
| ' | atom number | atom-map number (next free); View → Atom Numbers shows indices | ≈ |
| + − | charge | same | ✓ |
| g | select the atom | same | ✓ |
| = | nickname dialog | label editor (abbreviations, SMILES, any text) | ≈ |
| / ? | properties dialog | same | ✓ |
| > | sprout with labels | none | ✗ |

## Bond hotspot

| Key | ChemDraw | Penzene | |
|---|---|---|---|
| 1 2 3 | single, double, triple | same (2 on a double swaps its side) | ✓ |
| w / h W | wedged / hashed wedge | same | ✓ |
| H | hashed (not tapered) | hashed wedge | ≈ |
| b B | bold, double bold | same | ✓ |
| d D | dashed, dashed double | same | ✓ |
| y | wavy | same | ✓ |
| l c r | double bond left / centre / right | same | ✓ |
| f | bring to front | same (crossed bonds get a gap) | ✓ |
| a z | fuse benzene / cyclopentadiene | same | ✓ |
| v 4–8 | fuse 3–8 ring | same | ✓ |
| 9 0 | fuse chair cyclohexane | same | ✓ |
| g | select the bond | same (its two atoms) | ✓ |
| — | — | i, p, P: interaction and partial (TS) bonds | Penzene only |

## No hotspot (tool keys)

| Key | ChemDraw | Penzene | |
|---|---|---|---|
| Space | marquee | select tool | ✓ |
| x X | bond, chain | same | ✓ |
| j J | benzene, cyclopentadiene tools | j benzene; J not a tool | ≈ |
| e | arrows | same | ✓ |
| t | text | same | ✓ |
| T | brackets | Structure → Brackets (no tool) | ≈ |
| E | chemical symbols | Electrons and δ on the atom menu, `:` `*` keys | ≈ |
| g G | TLC plate, orbitals | none (outside Penzene's scope) | ✗ |

## Selection and moves (cheat sheet)

| Keys | ChemDraw | Penzene | |
|---|---|---|---|
| Enter | selected molecule → hotspot | same | ✓ |
| Space / Tab | hotspot → its molecule | Space, same (Tab moves focus) | ≈ |
| ←↑→↓, Shift | nudge 1 / 10 | same, with a selection; otherwise they walk the hotspot | ✓ |
| Ctrl+arrows | copy across a reaction arrow | same | ✓ |
| Alt+arrows | rotate | same (15°) | ✓ |
| Shift+Alt+arrows | 3D rotation | same (15°, keeps stereo) | ✓ |

## Menu shortcuts that differ

Penzene uses the platform's standard keys (Ctrl on Windows and Linux, ⌘ on macOS)
for File and Edit commands, as ChemDraw does. ChemDraw's own extras (⌘D copy as
CDXML, ⌥⌘C copy as SMILES, ⌥⌘P paste as SMILES) aren't bound: Copy already
carries SMILES, MOL and (where RDKit can write it, not yet on Windows, #185) CDX together,
and Paste reads any of them.
