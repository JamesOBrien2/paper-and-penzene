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
`pixi run screenshots` retakes the documentation's screenshots, light and dark, into `docs/_static`;
run it when the interface changes.

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

## Changelog

`CHANGELOG.md` has one section per release, newest first, under an `Unreleased` section that
collects changes as they merge. Write each entry for users, not developers:

- A main feature is a highlight: `- **Title**: what it does. <!-- icon: name -->`. The icon is a
  file in `resources/whatsnew/` (Tabler icons); the app's What's New window shows highlights as
  cards.
- Anything smaller is a plain bullet.

The release workflows take each release's notes from its section, so the GitHub release says the
same as the app.

## Releasing

1. Everything in the milestone is merged, or moved to a later one.
2. The checks pass on `main`: the nightly fuzz and AddressSanitizer runs, and a code review of
   the changes since the last release.
3. Bump `project(penzene VERSION X.Y.Z)` in `CMakeLists.txt`.
4. In `CHANGELOG.md`, rename `Unreleased` to `X.Y.Z (YYYY-MM-DD)`.
5. Add `<release version="X.Y.Z" date="YYYY-MM-DD"/>` at the top of the releases in
   `packaging/linux/io.github.jamesobrien2.penzene.metainfo.xml`.
6. Save a drawing with the new version as `tests/data/penz/vX.Y.Z.penz`, so later releases are
   tested against it (see [the .penz format](files.md)).
7. Merge that, then tag it: `git tag vX.Y.Z && git push origin vX.Y.Z`.
8. Check the release once the workflows finish:
   - It has 9 assets: 2 macOS `.dmg`, the Windows `.zip` and installer, the Linux AppImage, and
     4 wheels.
   - The new version is on [PyPI](https://pypi.org/project/penzene/).
   - The documentation built on Read the Docs.
9. Close the milestone and update the project board.
