# Command line

`penzene --render` draws structures without opening a window, for scripts, notebooks and batch
figures:

```sh
penzene --render aspirin.mol aspirin.svg                    # one file
penzene --render "CC(=O)Oc1ccccc1C(=O)O" aspirin.pdf        # a SMILES string
penzene --render library.sdf hits.smi --out figs --format png --drawing-style RSC --clean
```

Inputs can be SMILES strings or `.smi`, `.sdf`, `.inchi`, `.mol`, `.penz`, `.rxn` and `.cdxml` files.
Each record of a `.smi`, `.sdf` or `.inchi` becomes its own file, named after the record. The
paths written are printed one per line.

| Option | |
|---|---|
| `--out DIR` | write into this folder (otherwise give one input and one output path) |
| `--format svg\|png\|pdf` | with `--out` (default svg) |
| `--drawing-style NAME` | `ACS 1996`, `JDP` or `RSC` |
| `--clean` | lay each structure out afresh |
| `--version`, `--help` | |

## Descriptor tables

`penzene --descriptors` writes a CSV with one row per structure, for building datasets:

```sh
penzene --descriptors library.sdf --out library.csv
penzene --descriptors "CC(=O)Oc1ccccc1C(=O)O" --columns name,smiles,clogp,tpsa
```

It takes the same inputs as `--render`, and writes to standard output unless given `--out FILE.csv`.
`--columns` picks columns, in the order given. In the app, **File → Export Descriptors…** writes the
same table for each molecule in the selection, or on the page.

A structure that can't be read or isn't valid chemistry keeps its row, with the reason under `error`.
Values come from RDKit, as in the Properties panel:

| Column | |
|---|---|
| `id` | row number, from 1 |
| `name` | the record's name: an SDF title, the second column of a `.smi` line, or the file name |
| `smiles` | canonical SMILES |
| `inchikey` | standard InChIKey |
| `formula` | molecular formula, Hill order |
| `mw`, `exact_mass` | average molecular weight (2 decimals) and monoisotopic mass (4 decimals) |
| `clogp` | Wildman–Crippen logP |
| `tpsa` | topological polar surface area, Å² (Ertl; N and O only) |
| `hbd`, `hba` | Lipinski hydrogen-bond donors and acceptors |
| `rotatable_bonds` | RDKit's default rotatable-bond count |
| `heavy_atoms` | non-hydrogen atoms |
| `aromatic_rings` | aromatic rings in the smallest set of smallest rings |
| `lipinski_violations` | how many of MW > 500, logP > 5, HBD > 5, HBA > 10 |
| `veber` | `true` if rotatable bonds ≤ 10 and TPSA ≤ 140 |
| `error` | empty, `unreadable` or `not valid chemistry` |

On macOS the program is `Penzene.app/Contents/MacOS/penzene`.
