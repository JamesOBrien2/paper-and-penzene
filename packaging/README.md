# Packaging

Files prepared for distribution channels outside this repository. Each is submitted to its
channel by hand, with the owner's approval.

| Channel | File | Submitted? |
|---|---|---|
| conda-forge | `conda-forge/recipe.yaml`: the Python module, built against conda-forge's RDKit and Qt | not yet |

Test a recipe locally with rattler-build (a variants file gives what conda-forge's pinning
normally supplies, e.g. `c_stdlib` on macOS):

```sh
pixi exec rattler-build build -r packaging/conda-forge/recipe.yaml -c conda-forge -m variants.yaml
```
