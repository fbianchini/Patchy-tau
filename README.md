# PatchyTau

Theory power spectra for **patchy reionization** and their cross-correlation
with **CMB lensing**, implemented in Julia.

> These spectra were used in
> [Bianchini & Millea (2022), arXiv:2210.10893](https://arxiv.org/abs/2210.10893).

During the epoch of reionization the free-electron field is spatially
inhomogeneous ("patchy"), which imprints fluctuations in the Thomson optical
depth along each line of sight. `PatchyTau` models this with an ionized-bubble
halo model and computes:

- `Cℓττ` — the optical-depth angular auto-spectrum,
- `Cℓτϕ` — the optical-depth × CMB-lensing-potential cross-spectrum,
- `Cℓκκ` — the CMB lensing convergence spectrum (used as a cross-check),

together with the underlying 3-D free-electron power spectra `Pee`, `Pδe`.

The physics follows the bubble model of
[Wang & Hu (2006)](https://arxiv.org/abs/astro-ph/0607652); the τ–ϕ
cross-spectrum is of the kind discussed in
[Feng & Holder (2018)](https://arxiv.org/abs/1808.01592).

This package is a refactor of the original
`cmb_tau_lensing_theory_spectra.ipynb` notebook into reusable, documented
source files.

## Installation

The package is not registered; install it locally. From the repository root:

```julia
using Pkg
Pkg.activate(".")
Pkg.add(["QuadGK", "Trapz", "Dierckx", "MatterPower"])  # core dependencies
```

or, equivalently, run the helper script which also fills in the correct
dependency UUIDs in `Project.toml`:

```bash
julia scripts/setup.jl              # core dependencies only
julia scripts/setup.jl --with-plots # also installs PyPlot, JLD2, CMBLensing
```

The `--with-plots` extras (`PyPlot`, `JLD2`, `CMBLensing`) are needed only by
the reproduction scripts in `scripts/`, not by the core package.

## Quick start

```julia
using PatchyTau

c = Cosmology()                         # fiducial parameters + precomputed tables

# Reionization history
ionized_fraction(8; zre=10, Δz=4)       # mean ionized fraction x̄ₑ(z)
optical_depth(c, 30; zre=10, Δz=4)      # homogeneous optical depth τ̄(z)

# Angular spectra (R̄ = median bubble radius [Mpc], σlnR = log-width)
clττ = Cℓττ(c; R̄=5, σlnR=log(2), zre=10, Δz=4)   # returns a Spline1D in ℓ
clττ(1000)                                          # evaluate at ℓ = 1000

ℓ, clτϕ = Cℓτϕ(c; R̄=5, σlnR=log(2), zre=10, Δz=4) # (multipoles, Cℓ)
```

Run the bundled tour with:

```bash
julia --project=. examples/quickstart.jl
```

## Package layout

| File | Contents |
|------|----------|
| `src/constants.jl`        | physical constants and unit conversions |
| `src/cosmology.jl`        | `Cosmology` type, `hubble`, `comoving_distance`, `growth_factor`, `matter_power` |
| `src/reionization.jl`     | `ionized_fraction`, `optical_depth` |
| `src/bubbles.jl`          | bubble-size distribution, real/Fourier filters, `sigma2`, `mass_threshold` |
| `src/electron_spectra.jl` | `Pee`, `Pee1h`, `Pee2h`, `Pδe` |
| `src/angular_spectra.jl`  | `Wκ`, `Cℓκκ`, `Cℓττ`, `Cℓτϕ` |
| `scripts/`                | dependency setup + notebook reproduction (spectra files, paper figure) |
| `examples/`               | minimal usage example |

## Model parameters

The main physical knobs, passed as keyword arguments to the spectrum routines:

| Symbol | Keyword | Meaning |
|--------|---------|---------|
| R̄      | `R̄`     | median ionized-bubble radius [Mpc] |
| σ_lnR  | `σlnR`  | log-width of the bubble-size distribution |
| z_re   | `zre`   | midpoint redshift of reionization |
| Δz     | `Δz`    | width of the reionization transition |
| b      | `b`     | linear bias of the ionized regions |

Cosmological parameters are set when constructing `Cosmology(; As, ns, Ωm, Ωb,
h₀, …)`; the defaults reproduce the fiducial model of the original notebook.

## Reproducing the notebook outputs

With the plotting extras installed:

```bash
julia --project=. scripts/generate_spectra.jl   # writes data/theory_spectra_{optimistic,pessimistic}.jld2
julia --project=. scripts/make_paper_plot.jl     # writes theory_spectra_paper.pdf
```

## Performance

- A `Cosmology` precomputes the linear growth solution, a χ↔z lookup spline
  and a 2-D `P(k, z)` spline, so construction is the main up-front cost; reuse
  the same object across calls.
- `Cℓττ` caches the bubble-averaged variance `⟨σ²_R⟩(z)` once per redshift (it
  is `k`-independent) rather than recomputing it for every multipole — the
  dominant cost. The `σ²b` keyword on `Pee`/`Pee1h`/`P̃δδ` exposes the same
  trick for custom loops.
- Further opt-in speedup: since `P(k, z) = D(z)² f(k)`, every variance integral
  scales as `D(z)²`, so `⟨σ²_R⟩(z) = (D(z)/D(z₀))² ⟨σ²_R⟩(z₀)`. You can compute
  the integral once and rescale across redshift via `growth_factor`.
- The `j`-loops in `Cℓττ`/`Cℓτϕ` are independent across redshift columns and
  parallelise cleanly with `Threads.@threads` (start Julia with `-t auto`).

## Notes

- The original notebook (`cmb_tau_lensing_theory_spectra.ipynb`) is kept in the
  repository as a reference; the theory line for `Cℓτϕ` was derived by hand and
  may contain errors (see the notebook's own caveat).
