<p align="center"><img src="resources/logo.svg" width="160" alt="Penzene logo: a benzene ring drawn in blue ink"></p>

<h1 align="center">Penzene</h1>

An open-source, native desktop chemical structure editor — a free alternative to ChemDraw.
Runs on macOS, Linux and Windows. Built with C++20, Qt 6 and [RDKit](https://www.rdkit.org).

> Status: beta. **Documentation: [penzene.readthedocs.io](https://penzene.readthedocs.io)**.
> See the [roadmap](https://github.com/users/JamesOBrien2/projects/1) and [issues](https://github.com/JamesOBrien2/penzene/issues).

## Install

Download the installer for your system from the [latest release](https://github.com/JamesOBrien2/penzene/releases/latest)
(macOS `.dmg`, Windows `setup.exe` or `.zip`, Linux `.AppImage`), or `pip install penzene` for the Python package.

The installers aren't signed with a paid certificate, so the first time you open Penzene:
- **macOS:** if it "can't be opened", go to System Settings → Privacy & Security and click **Open Anyway**.
- **Windows:** if SmartScreen appears, click **More info → Run anyway**.
- **Linux:** `chmod +x penzene-linux-x86_64.AppImage`, then run it.

<p align="center"><img src="docs/screenshot.png" width="720" alt="Penzene main window showing aspirin"></p>

## Features

- Draw atoms, bonds (single/double/triple, wedge/hash), chains and rings with the mouse
- **ChemDraw-style hotkeys**: point at an atom or bond and type. `1111` draws a chain,
  `2` sprouts a carbonyl, `a` a phenyl, `O` an OMe. The arrow keys walk the molecule
  (Help → Keyboard Shortcuts)
- Open/save `.penz`, MOL and SDF; open ChemDraw `.cdxml` (molecules, arrows, text) and `.cdx` (molecules); paste or import SMILES; Clean structure (RDKit)
- Reaction, equilibrium, resonance, retrosynthesis, curved and fishhook arrows; text with automatic formula subscripts (H2O → H₂O)
- Live formula, MW and exact mass for the selection; copy as InChI / InChIKey
- Export SVG, PNG and PDF; copy as image + MOL + SMILES
- Abbreviations (Me, OMe, CO2Me, Boc, TBS, Ts, Bpin…) drawn as labels and expanded for chemistry; Structure → Expand draws them out
- Ring fill: click inside a ring with the fill tool to shade it
- Themes: follow the OS, Light, Dark, or Catppuccin Latte / Frappé / Macchiato / Mocha (View → Theme); exports always stay black on clear
- ACS 1996 drawing style, implicit hydrogens and valence warnings

## Command line

Render without opening a window, for scripts, notebooks and batch figures:

```sh
penzene --render aspirin.mol aspirin.svg                    # one file
penzene --render "CC(=O)Oc1ccccc1C(=O)O" aspirin.pdf        # a SMILES string
penzene --render library.sdf hits.smi --out figs --format png --drawing-style RSC --clean
```

Inputs can be SMILES strings or `.smi`, `.sdf`, `.mol`, `.penz` and `.cdxml` files. In a `.smi` or `.sdf`, every record becomes its own file, named after the record's name. The paths written are printed one per line. On macOS the binary is `Penzene.app/Contents/MacOS/penzene`.

## Python

The same engine as the app, from Python 3.12+ on macOS, Linux (glibc 2.17+) and Windows:

```sh
pip install penzene
```

```python
import penzene as pz
doc = pz.from_smiles("CC(=O)Oc1ccccc1C(=O)O")   # aspirin
doc.clean(); doc.export("aspirin.svg")          # identical to the app's export
doc.formula, doc.mw                             # ('C9H8O4', 180.16)
doc                                             # renders inline in Jupyter
```

More in [python/README.md](python/README.md).

## Build

Requires [pixi](https://pixi.sh). It fetches Qt, RDKit and the toolchain from conda-forge.

```sh
pixi run build   # build/bin/penzene
pixi run test
pixi run run     # launch the app
pixi run install-app   # macOS: self-contained app in ~/Applications
```

## Roadmap

| Milestone | Highlights |
|---|---|
| v0.1–v0.3 | Drawing and ChemDraw hotkeys, MOL/SDF/`.penz`, SMILES, Clean, SVG/PNG/PDF; schemes (arrows, text), abbreviations, formula/MW/InChI, CDXML import; themes, ring fill, ACS/JDP/RSC styles |
| v0.4 | Python package (`import penzene`), batch command line |
| v0.5 | Polish and correctness: flip/align, context menus, colouring, preferences, autosave |
| v0.6 | Chemistry: stereo labels, structure checks, properties panel, name ↔ structure |
| v0.7 | Interop: editable exports, reactions, more formats, paste from ChemDraw |
| v0.8 | Templates, brackets and electron dots, shapes, projections |
| v0.9 | Beta: documentation site, signed installers, PyPI, accessibility |
| v1.0 | Stable `.penz` format and Python API |

Details and progress on the [project board](https://github.com/users/JamesOBrien2/projects/1).

## Acknowledgements

UX and tool set inspired by [Ketcher](https://github.com/epam/ketcher) (Apache-2.0).
Chemistry by RDKit (BSD-3). GUI by Qt (LGPL-3.0).

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
