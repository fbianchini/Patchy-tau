#!/usr/bin/env julia
#
# Minimal, plotting-free tour of the package. Run with:
#
#     julia --project=. examples/quickstart.jl

using PatchyTau

# Fiducial background cosmology (precomputes growth + matter power spline).
c = Cosmology()

println("Derived quantities")
println("  χ(z=1089) = ", round(c.χlss, digits=1), " Mpc")
println("  σ8-ish: √σ²(R=8 Mpc, z=0) = ", round(√(sigma2(c, 0.0, 8.0)), digits=4))

# Reionization history.
println("\nReionization (zre=10, Δz=4)")
for z in (6, 8, 10, 12)
    xe = ionized_fraction(z; zre=10, Δz=4)
    println("  x̄ₑ(z=$z) = ", round(xe, digits=3))
end
println("  τ̄(z=30)  = ", round(optical_depth(c, 30; zre=10, Δz=4), digits=4))

# Angular spectra for the "optimistic" bubble model.
println("\nAngular spectra (R̄=5, σlnR=log2, zre=10, Δz=4)")
clττ = Cℓττ(c; R̄=5, σlnR=log(2), zre=10, Δz=4)
ℓ, clτϕ = Cℓτϕ(c; R̄=5, σlnR=log(2), zre=10, Δz=4)
for L in (100, 500, 1000, 2000)
    println("  ℓ=$L  Cℓττ = ", clττ(L))
end
println("  Cℓτϕ sampled at ℓ = ", Int.(ℓ[1:5]), " ...")
