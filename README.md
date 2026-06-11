# Atmospheric Model for Asymptotic Giant Branch Stars
#### Authors: Daniel Kitzmann, Joachim Stock ####

This is the documentation branch for the AGB stationary dust-driven wind model.

The documentation is written in [Sphinx](https://www.sphinx-doc.org).

## Layout

- `doc_src/`   — reStructuredText sources (version controlled)
- `doc_built/` — generated HTML/PDF output (git-ignored)

## Building

Requires Sphinx and the Furo theme:

```
pip install sphinx furo
```

Then build the HTML:

```
cd doc_src
make html
```

The rendered site is written to `../doc_built/html`; open
`doc_built/html/index.html` in a browser.

## Contents

The documentation covers both the **configuration / usage** (installation,
running, the full `config.toml` reference, file formats) and the **theory and
modelling** (radiative transfer, radiative-equilibrium temperature correction,
gas chemistry, Gail & Sedlmayr dust formation, opacities, wind hydrodynamics,
and the numerical robustness methods), plus an implementation guide to the
source tree.
