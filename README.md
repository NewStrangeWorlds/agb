# Atmospheric Model for Asymptotic Giant Branch Stars
#### Authors: Daniel Kitzmann, Joachim Stock ####


# Overview #

This is a self-consistent model for the **stationary, spherically symmetric,
dust-driven wind of a carbon-rich asymptotic giant branch (AGB) star**. Carbon
that is not locked into CO condenses into amorphous-carbon grains in the cool,
extended atmosphere; the radiative acceleration on these grains, once it exceeds
gravity, drives a slow, massive outflow. The model solves the coupled
radiation–chemistry–dust–hydrodynamics problem for the steady-state structure
and the mass-loss rate.

The code couples:

* a **spherical, frequency-dependent radiative transfer** solver (Rybicki–Hummer
  variable-Eddington-factor moment method);
* a **radiative-equilibrium temperature** determination (Unsöld–Lucy correction
  and a full-linearisation Newton corrector);
* **equilibrium gas-phase chemistry** (FastChem);
* **carbon-dust nucleation and growth** following Gail & Sedlmayr (moment method
  with self-consistent, differential carbon depletion);
* **stationary wind hydrodynamics** (Melia `Φ` transform with a mass-loss
  eigenvalue, plus an optional Henyey global-Newton structure solver),

iterated to a self-consistent stationary outflow.


# Requirements #

* A C++17 compiler (GCC or Clang) with **OpenMP**.
* **CMake** ≥ 3.10.
* Network access at configure time — the build automatically fetches its
  dependencies.

The following are downloaded and built by CMake (`FetchContent` /
`ExternalProject`), so no manual installation is needed:

| Dependency | Version | Purpose |
|------------|---------|---------|
| [FastChem](https://github.com/exoclime/fastchem) | 4.0.2 | equilibrium gas-phase chemistry |
| [LX-MIE](https://github.com/daniel-kitzmann/LX-MIE) | pinned | Mie theory for the dust opacities |
| [Eigen](https://gitlab.com/libeigen/eigen) | 3.4.0 | linear algebra (Henyey solver) |
| [toml++](https://github.com/marzer/tomlplusplus) | 3.4.0 | parsing the `config.toml` |
| [CppAD](https://github.com/coin-or/CppAD) | 20250000.3 | automatic differentiation (Henyey Jacobian) |

In addition the model needs, at run time:

* an **opacity database** of cross sections (HELIOS-K format), and
* a **dust refractive-index** file (e.g. amorphous carbon).

> **Note:** Do **not** add `-march=native` to the build flags — FastChem's
> bundled Eigen breaks under it. The Release flags are deliberately limited to
> `-O3 -DNDEBUG -funroll-loops` (and `-ffast-math` is omitted) to keep the
> radiative-transfer results bit-reproducible.


# Building #

Standard out-of-source CMake build:

```bash
mkdir build
cd build
cmake ..
make -j
```

The first `cmake ..` clones and builds the dependencies; subsequent builds reuse
them. The executable is written into the source directory as `agb_wind`.

A built-in self-test validates the Henyey structure solver (automatic-
differentiation Jacobian, Newton/eigenvalue solve, grey bootstrap) without
running a full model:

```bash
./agb_wind --selftest
```


# Running #

Pass the path to a **model folder** containing a `config.toml`:

```bash
./agb_wind IRC10216/
```

All input and output paths in the configuration are interpreted relative to the
model folder. Parallelism is controlled through `OMP_NUM_THREADS`; results are
bit-reproducible run-to-run for a fixed thread count.

A typical model folder:

```
IRC10216/
├── config.toml                 # model configuration (required)
├── fastchem_parameters.dat     # FastChem element-abundance file
├── atmosphere_converged.dat    # starting structure (or use "grey")
├── atmosphere.dat              # output: converged structure
├── spectrum.dat                # output: emergent spectrum
├── dust.dat                    # output: dust distribution
└── hydro.dat                   # output: wind structure
```


# Configuration #

The model is configured by a single [TOML](https://toml.io) file,
`config.toml`, in the model folder. It is organised into blocks — `[star]`,
`[model]`, `[spectral_grid]`, `[radiative_transfer]`, `[temperature]`,
`[hydrodynamics]`, `[grid]`, `[dust]`, `[movable_grid]`, `[opacity]`,
`[output]` — each with sensible built-in defaults. A complete, annotated example
ships in `IRC10216/config.toml`, and every key is documented in the
[full documentation](#documentation).

The starting model is either a path to an existing structure file or the literal
string `"grey"`, which builds a hydrostatic Lucy grey structure on a logarithmic
radial grid with no input file.


# Output #

Depending on the `[output]` block, the model writes:

* the converged **atmospheric structure** (radius, density, pressure, gas/dust
  temperature, velocity);
* the emergent **spectrum** (wavelength vs. flux);
* the **dust distribution** (degree of condensation, grain number density and
  radius, nucleation rate, growth timescale, moments K0…K3);
* the **wind structure** (velocity, sound speed, radiative acceleration α, the
  Melia `Φ` variable).

File formats are documented in the [documentation](#documentation).


# Documentation #

Extensive documentation — covering both the configuration/usage and the
theory/modelling of every module — is maintained as a Sphinx project on the
`docs` branch of the companion `agb_docs` repository. Build it with:

```bash
cd doc_src
make html        # output in ../doc_built/html
```

(requires `pip install sphinx furo`).


# References #

* **Kitzmann, D.**, diploma thesis — the spherical variable-Eddington-factor
  radiative-transfer scheme and the full-linearisation temperature correction.
* **Gail, H.-P. & Sedlmayr, E.** — dust nucleation and the moment method.
* **Winters, J. M.**, PhD thesis — the coupled stationary-wind scheme and the
  differential carbon depletion.
* **Melia, F. (1988)** — the `Φ` transform removing the sonic-point singularity.
* **Dominik, C. (1990)** — the Henyey-type relaxation of the wind structure.


# Status #

Research software under active development. The default configuration —
Melia `Φ` shooting solver with an eigenvalue mass-loss rate, Taylor
radiative-transfer discretisation, and the two-phase Unsöld–Lucy →
linearisation temperature corrector — converges the IRC+10216 model and
reproduces its observed mass-loss rate.
