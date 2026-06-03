#!/usr/bin/env julia
#
# Reproduce the "optimistic" and "pessimistic" theory-spectra files from the
# original notebook and save them as JLD2.
#
# Requires the plotting/extras packages (see `scripts/setup.jl --with-plots`):
# this script uses CMBLensing for the ℓ-extrapolation helper and JLD2 to save.
#
#     julia --project=. scripts/generate_spectra.jl

using PatchyTau
using CMBLensing: extrapolate_Cℓs, ℓ⁴
using JLD2

const OUTDIR = joinpath(@__DIR__, "..", "data")
mkpath(OUTDIR)

c = Cosmology()

"""
    bandpowers(c; R̄, σlnR, zre, Δz)

Build the ϕϕ / τϕ / ττ spectra on ℓ = 2:8000 for a given bubble + reionization
model, matching the construction in the notebook.
"""
function bandpowers(c; R̄, σlnR, zre, Δz)
    ℓ, clτϕ = Cℓτϕ(c; R̄=R̄, σlnR=σlnR, zre=zre, Δz=Δz)
    Cℓτϕ_ = extrapolate_Cℓs(2:8000, ℓ, clτϕ)
    Cℓττ_ = extrapolate_Cℓs(2:8000, 1:6000, Cℓττ(c; R̄=R̄, σlnR=σlnR, zre=zre, Δz=Δz)(1:6000))
    # GetCℓκκ returns the convergence spectrum on ℓ = 2:2000; convert κκ -> ϕϕ
    # via Cℓϕϕ = 4 Cℓκκ / ℓ⁴.
    Cℓϕϕ_ = extrapolate_Cℓs(2:8000, 1:1999, Cℓκκ(c, lensing_kernel(c))) / ℓ⁴ * 4
    return Cℓϕϕ_, Cℓτϕ_, Cℓττ_
end

# NB: use lowercase locals so we don't shadow the `Cℓ*` package functions that
# `bandpowers` calls; save them back under the canonical `Cℓ*` keys.
@info "Optimistic model (R̄=5, σlnR=log2, zre=10, Δz=4)"
clϕϕ, clτϕ, clττ = bandpowers(c; R̄=5, σlnR=log(2), zre=10, Δz=4)
jldsave(joinpath(OUTDIR, "theory_spectra_optimistic.jld2"); Cℓϕϕ=clϕϕ, Cℓτϕ=clτϕ, Cℓττ=clττ)

@info "Pessimistic model (R̄=1, σlnR=0.1, zre=7, Δz=1)"
clϕϕ, clτϕ, clττ = bandpowers(c; R̄=1, σlnR=0.1, zre=7, Δz=1)
jldsave(joinpath(OUTDIR, "theory_spectra_pessimistic.jld2"); Cℓϕϕ=clϕϕ, Cℓτϕ=clτϕ, Cℓττ=clττ)

@info "Saved spectra to $OUTDIR"
