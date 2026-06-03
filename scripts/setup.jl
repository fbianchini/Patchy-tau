#!/usr/bin/env julia
#
# One-off dependency installer.
#
# Run from the repository root:
#
#     julia scripts/setup.jl
#
# This activates the package environment and `Pkg.add`s every dependency,
# which writes the correct UUIDs into Project.toml/Manifest.toml. The
# `--with-plots` flag additionally installs the (heavier) packages needed by
# the reproduction scripts under `scripts/`.

using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

# Core runtime dependencies of the PatchyTau package.
core = ["QuadGK", "Trapz", "Dierckx", "MatterPower"]

# Extra packages used only by the example / figure scripts.
plotting = ["PyPlot", "JLD2", "CMBLensing"]

deps = core
if "--with-plots" in ARGS
    append!(deps, plotting)
end

@info "Installing dependencies" deps
Pkg.add(deps)
Pkg.instantiate()
@info "Done. Try:  julia --project=. examples/quickstart.jl"
