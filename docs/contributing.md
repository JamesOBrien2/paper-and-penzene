# Contributing

## Building

Requires [pixi](https://pixi.sh), which fetches Qt, RDKit and the toolchain from conda-forge.

```sh
pixi run build         # build/bin/penzene
pixi run test          # C++ tests, the CLI and the Python package
pixi run run           # launch the app
pixi run install-app   # macOS: a self-contained app in ~/Applications
```

`python docs/generate.py` rebuilds the generated reference pages (keys, Python API); a test fails
if they're out of date.

- Pick an issue from the project board; comment to claim it.
- Branch `feat/<issue>-<slug>` or `fix/<issue>-<slug>`; open a PR with `Closes #<issue>`.
- `pixi run test` must pass; CI builds macOS, Linux and Windows.
- Keep it small: one issue per PR, no speculative abstractions.

## Translating

Every interface string goes through `tr()`. A translation is one Qt Linguist file,
`translations/penzene_<lang>.ts` (for example `penzene_de.ts`):

```sh
pixi run lupdate -ts translations/penzene_de.ts
```

collects the strings (and keeps what's already translated). Translate them in Qt Linguist or any
`.ts` editor, rebuild, and pick the language under **Preferences → Language**. The build embeds
every `.ts` file it finds, so a pull request with just that file is all it takes. See
[translations/README.md](https://github.com/JamesOBrien2/penzene/blob/main/translations/README.md).

## Code map

| File | Role |
|---|---|
| `src/Document.*` | Plain atom/bond model + `.penz` JSON |
| `src/Chem.*` | The only RDKit bridge (SMILES, MOL, coordinates) |
| `src/Canvas.*` | Drawing surface and tools |
| `src/MainWindow.*` | Menus, files, clipboard, export |

## Project board views

Status columns are Todo / In Progress / Done. Views (set up once, by hand):
**Features** (`label:feature`), **Bugs** (`label:bug`), **Roadmap** (table, group by Milestone).
