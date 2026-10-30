# exploreDE: instructions for an AI assistant

The instructions for helping a user build an exploreDE `.h5ad` from their own
data live in [`inst/skills/exploreDE-byod/SKILL.md`](inst/skills/exploreDE-byod/SKILL.md).
From an installed package:

```r
system.file("skills/exploreDE-byod/SKILL.md", package = "exploreDE")
```

The human-oriented guide is `vignette("bring-your-own-data", package = "exploreDE")`.

The file formats are specified in `inst/spec/`:

- [`exploreDE-h5ad-1.0.md`](inst/spec/exploreDE-h5ad-1.0.md): The `.h5ad` layout, for writing files outside R.
- [`exploreDE-session-2.md`](inst/spec/exploreDE-session-2.md): The saved-session `.json` file.
