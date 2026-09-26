<p align="center"><img src="resources/logo.svg" width="160" alt="Penzene logo: a benzene ring drawn in teal ink"></p>

<h1 align="center">Penzene</h1>

An open-source, native desktop chemical structure editor.
Runs on macOS, Linux and Windows. Built with C++20, Qt 6 and [RDKit](https://www.rdkit.org).

> **Documentation: [penzene.readthedocs.io](https://penzene.readthedocs.io)**.
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

- **Keyboard drawing**: point at an atom or bond and type. `1111` draws a chain, `2` a carbonyl,
  `a` a phenyl, `O` an OMe, and the arrow keys walk the molecule.
- **Works with ChemDraw**: opens `.cdxml` and `.cdx` drawings with their arrows, text and shapes,
  and pastes from ChemDraw.
- **Reaction schemes and mechanisms**: reaction, equilibrium and retrosynthesis arrows, curved and
  fishhook arrows, lone pairs, radicals and partial charges.
- **Chemistry built in**: formula, mass, cLogP, TPSA and drug-likeness checks as you draw, stereo
  labels, name ↔ structure, and descriptor tables for whole datasets.
- **Figures that stay editable**: SVG, PNG and PDF exports carry the drawing, so they open or paste
  back into Penzene. Copy pastes sharp vector images into Word, PowerPoint and Keynote.
- **Figures at final size**: lay a figure out on a page or a journal's column width.
- **Scriptable**: the same engine from the command line and from Python.

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

Plans and progress are on the [project board](https://github.com/users/JamesOBrien2/projects/1).

## Acknowledgements

UX and tool set inspired by [Ketcher](https://github.com/epam/ketcher) (Apache-2.0).
Chemistry by RDKit (BSD-3). GUI by Qt (LGPL-3.0).

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
