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

On macOS the program is `Penzene.app/Contents/MacOS/penzene`.
