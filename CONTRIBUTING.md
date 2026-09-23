# Contributing

- Pick an issue from the project board; comment to claim it.
- Branch `feat/<issue>-<slug>` or `fix/<issue>-<slug>`; open a PR with `Closes #<issue>`.
- `pixi run test` must pass; CI builds macOS, Linux and Windows.
- Keep it small: one issue per PR, no speculative abstractions.

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
