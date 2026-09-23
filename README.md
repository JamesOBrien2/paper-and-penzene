# Paper & Penzene

An open-source, native desktop chemical structure editor — a free alternative to ChemDraw.
Runs on macOS, Linux and Windows. Built with C++20, Qt 6 and [RDKit](https://www.rdkit.org).

> Status: early development. See the [roadmap](https://github.com/users/JamesOBrien2/projects) and
> [issues](https://github.com/JamesOBrien2/paper-and-penzene/issues).

## Build

Requires [pixi](https://pixi.sh). It fetches Qt, RDKit and the toolchain from conda-forge.

```sh
pixi run build   # build/bin/penzene
pixi run test
pixi run run     # launch the app
```

## Roadmap

| Milestone | Highlights |
|---|---|
| v0.1 | Draw atoms/bonds/rings, undo, MOL/SDF + `.penz`, SMILES paste, Clean, SVG/PNG/PDF export, clipboard |
| v0.2 | Dark/light mode, Catppuccin themes, ring fill, text, arrows, formula/MW, InChI, CDXML import |
| Later | Structure → name, name → structure, style presets, templates, printing |

## Acknowledgements

UX and tool set inspired by [Ketcher](https://github.com/epam/ketcher) (Apache-2.0).
Chemistry by RDKit (BSD-3). GUI by Qt (LGPL-3.0).

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
