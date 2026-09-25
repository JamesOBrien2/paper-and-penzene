# Translations

Penzene's interface strings all go through `tr()`. Each `penzene_<lang>.ts` here (Qt Linguist
format, e.g. `penzene_de.ts`, `penzene_pt_BR.ts`) is compiled into the app at build time, and the
language then appears under **Edit → Preferences → Language**.

To start or update a translation, from the repository root:

```sh
pixi run lupdate -ts translations/penzene_de.ts
```

This collects every string from `src/` into the file (new strings are added, existing
translations kept). Translate them in Qt Linguist or any `.ts` editor, rebuild with `pixi run build`, and choose the language in Preferences.
Open a pull request with the `.ts` file; nothing else needs to change.
