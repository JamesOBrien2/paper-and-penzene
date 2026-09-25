# Penzene

::::{container} hero

:::{container} tagline
An open-source, native desktop chemical structure editor, and a free alternative to ChemDraw.
It runs on macOS, Linux and Windows, and does its chemistry with [RDKit](https://www.rdkit.org).
:::

```{image} _static/screenshot-light.png
:alt: Penzene with aspirin
:class: shot only-light
```

```{image} _static/screenshot-dark.png
:alt: Penzene with aspirin, dark theme
:class: shot only-dark
```

::::

The same drawing from a script:

```sh
penzene --render "CC(=O)Oc1ccccc1C(=O)O" aspirin.svg
```

::::{container} features

:::{container} card
**Draw fast**

With the mouse or [the keyboard](keys.md): point at an atom and type `2` for a carbonyl or `a` for a phenyl.
:::

:::{container} card
**Chemistry built in**

Live formula and mass, stereo labels, a properties panel, structure checks, and names from PubChem.
:::

:::{container} card
**Figures for papers**

ACS, RSC and JDP styles, reaction schemes, electron-pushing arrows, templates, and exports that reopen as editable drawings.
:::

:::{container} card
**Works with ChemDraw**

Open and save CDXML and CDX, and paste from ChemDraw.
:::

:::{container} card
**Scriptable**

The same engine from [Python](python.md) and the [command line](cli.md).
:::

:::{container} card
**Free and open**

GPL-3.0, with [the source on GitHub](https://github.com/JamesOBrien2/penzene).
:::

::::

Start with [Installing](install.md), then [Drawing](drawing.md).

```{toctree}
:caption: Getting started
:hidden:

install
drawing
keys
```

```{toctree}
:caption: User guide
:hidden:

chemistry
figures
files
```

```{toctree}
:caption: Scripting
:hidden:

cli
python
python-api
```

```{toctree}
:caption: Project
:hidden:

contributing
Changelog <https://github.com/JamesOBrien2/penzene/blob/main/CHANGELOG.md>
```
