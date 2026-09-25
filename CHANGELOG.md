# Changelog

What's new in each release. The app shows the newest section once after an update.

## Unreleased

- Documentation site: [penzene.readthedocs.io](https://penzene.readthedocs.io).
- Help → Check for Updates; an optional weekly check (off unless turned on in Preferences).
- A "What's new" window after an update.
- ChemDraw files: plain lines are no longer imported twice; lone pairs come through.

## 0.8.0 (2026-09-25)

- Template library (View → Templates): amino acids, sugars, nucleobases, scaffolds, and your own.
- Haworth, Fischer and Newman projections; Haworth and Fischer give the right stereo.
- Interaction bonds (hydrogen bonds, contacts) and partial bonds for transition states.
- Lone pairs, radicals, δ+ and δ−, and brackets with a subscript.
- Lines, boxes, rounded boxes and ellipses.
- Rotate in 3D (Shift+Alt+drag), keeping stereo; stretch and squash with handles; Structure → Transform.
- Structure → Arrange Scheme tidies a reaction scheme.
- R-groups (R1, R2…) and generic atoms (X, Ar); attachment points and η-bonded rings.
- Atom Properties; bring a bond to the front so crossings show a gap.
- Arrow keys nudge a selection; Space and Enter move between the hotspot and its molecule.

## 0.7.0 (2026-09-24)

- Exported SVG and PNG reopen as editable drawings.
- Reactions: copy as reaction SMILES, open and save .rxn.
- Multi-record SDF, .smi and .inchi files open as a grid; MOL V3000; paste InChI.
- Save as ChemDraw CDXML (and CDX); paste from ChemDraw.
- Export scale and margin; page mode (A4, Letter, journal columns); printing.

## 0.6.0 (2026-09-24)

- Stereo labels: CIP (R)/(S) and (E)/(Z).
- Check Structure explains valence, stereo and label problems.
- Properties panel: cLogP, TPSA, H-bond donors and acceptors, rotatable bonds, Lipinski, Veber.
- Carbon and hydrogen display options; add or remove explicit hydrogens.
- Aromatic circles; atom numbers and reaction atom-map numbers.
- X and R on an atom, and free-text atom labels.
- Name to structure and structure to name via PubChem (online, on request).
- The version in the window title; `pip install penzene` from PyPI.

## 0.5.0 (2026-09-24)

- Preferences; recent files, autosave and crash recovery.
- Right-click menus; flip, align and distribute; colour atoms, bonds and text.
- Drag an atom onto another to merge; drag to size rings.
- Per-style bond lengths; fixes to Clean and to CDXML import.

## 0.4.0 (2026-09-24)

- The Python package: `import penzene`, with wheels for macOS, Linux and Windows.
- `penzene --render` for batch figures from the command line.

## 0.3.0 (2026-09-24)

- ACS, RSC and JDP drawing styles; ring fill; dark mode and Catppuccin themes.
- A new tool palette with a periodic-table picker.

## 0.2.0 (2026-09-24)

- Reaction arrows and text; abbreviations (Me, Ph, Boc…).
- Formula and molecular weight; InChI and InChIKey; CDXML import.
- Windows installer and an Intel Mac build.

## 0.1.0 (2026-09-23)

- The first release: draw with bonds, chains, rings and keyboard shortcuts; undo and redo.
- MOL, SDF, SMILES and .penz; Clean; implicit hydrogens.
- SVG, PNG and PDF export; copy as image, MOL and SMILES.
