# Projected (angular) power spectra in the Limber approximation:
#   Cℓκκ : CMB lensing convergence auto-spectrum (sanity check)
#   Cℓττ : patchy optical-depth auto-spectrum
#   Cℓτϕ : optical-depth x CMB-lensing-potential cross-spectrum

"""
    DEFAULT_ℓS

Multipoles at which `Cℓττ` / `Cℓτϕ` are evaluated by default before splining /
returning.
"""
const DEFAULT_ℓS = [10, 50, 100, 150, 200, 250, 300, 350, 400, 500, 700, 1000,
                    1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000]

"""
    Wκ(c, z)

CMB lensing convergence kernel W^κ(z) (dimensionless).
"""
Wκ(c::Cosmology, z::Real) =
    1.5 * c.Ωm / c_kms * (100c.h₀)^2 / hubble(c, z) * (1 + z) *
    comoving_distance(c, z) * (1 - comoving_distance(c, z) / c.χlss)

"""
    lensing_kernel(c)

Return the CMB lensing convergence kernel as a single-argument closure `z -> W^κ(z)`,
suitable for passing to [`Cℓκκ`](@ref).
"""
lensing_kernel(c::Cosmology) = z -> Wκ(c, z)

"""
    Cℓκκ(c, W₁, W₂=nothing; ℓmax=2000, zmin=1e-4, zmax=1089,
         kmin=1e-4, kmax=1e2, npts=100, limb_fact=0.5)

Limber-approximation angular auto/cross-spectrum of two convergence-like
kernels `W₁(z)`, `W₂(z)` (functions of redshift). With a single kernel the
auto-spectrum is returned. Result is a vector over `ℓ = 2:ℓmax`.
"""
function Cℓκκ(c::Cosmology, W₁, W₂=nothing; ℓmax=2000, zmin=1e-4, zmax=1089,
              kmin=1e-4, kmax=1e2, npts=100, limb_fact=0.5)
    W₂ === nothing && (W₂ = W₁)

    ℓrange = 2:ℓmax
    χmin = comoving_distance(c, zmin)
    χmax = comoving_distance(c, zmax)
    χs = collect(range(χmin, χmax, length=npts))
    zs = c.χtoz(χs)
    dz = diff(zs); insert!(dz, 1, 0)

    fact = [hubble(c, zs[j]) / χs[j]^2 / c_kms for j in eachindex(zs)]
    kern = [W₁(z) * W₂(z) for z in zs]

    # K[i,j] = (ℓᵢ + limb_fact) / χⱼ   ->   wavenumber for each (ℓ, z) pair
    K = [(ℓ + limb_fact) / χ for ℓ in ℓrange, χ in χs]

    Mᵢⱼ = zeros(length(ℓrange), length(zs))   # P_δδ(ℓᵢ/χⱼ; zⱼ)
    for j in eachindex(zs)
        Mᵢⱼ[:, j] = c.pkz_spl(K[:, j], fill(zs[j], length(ℓrange)))
    end

    Mᵢⱼ[K .< kmin] .= 0
    Mᵢⱼ[K .≥ kmax] .= 0

    vⱼ = kern .* fact .* dz
    return Mᵢⱼ * vⱼ
end

"""
    Cℓττ(c; ℓs=DEFAULT_ℓS, zmin=1e-4, zmax=20, kmin=1e-4, kmax=1e1,
         npts=25, limb_fact=0.5, R̄=5, σlnR=log(2), zre=10.0, Δz=2.0, b=6)

Patchy optical-depth angular auto-spectrum,

    Cℓττ = ∫ dχ (σT² n_p0² / a⁴ χ²) Pee(k = ℓ/χ, χ),

returned as a `Dierckx.Spline1D` in ℓ so it can be evaluated at arbitrary
multipoles.
"""
function Cℓττ(c::Cosmology; ℓs=nothing, zmin=1e-4, zmax=20, kmin=1e-4, kmax=1e1,
              npts=25, limb_fact=0.5, R̄=5, σlnR=log(2), zre=10.0, Δz=2.0, b=6)
    np₀ = (1 - 0.75c.Yp) * c.ρc₀ * c.Ωb / mp   # m⁻³

    χmin = comoving_distance(c, zmin)
    χmax = comoving_distance(c, zmax)
    χs = collect(range(χmin, χmax, length=npts))
    zs = c.χtoz(χs)
    dχ = diff(χs); insert!(dχ, 1, 0)

    ℓrange = ℓs === nothing ? DEFAULT_ℓS : ℓs

    kern = [(1 + zs[j])^4 / χs[j]^2 for j in eachindex(zs)]
    K = [(ℓ + limb_fact) / χ for ℓ in ℓrange, χ in χs]

    Mᵢⱼ = zeros(length(ℓrange), length(χs))   # Pee(ℓᵢ/χⱼ; χⱼ)
    for j in eachindex(zs)
        Mᵢⱼ[:, j] = [Pee(c, K[i, j], zs[j]; R̄=R̄, σlnR=σlnR, zre=zre, Δz=Δz, b=b)
                     for i in eachindex(ℓrange)]
    end

    Mᵢⱼ[K .< kmin] .= 0
    Mᵢⱼ[K .≥ kmax] .= 0

    vⱼ = kern .* dχ
    cl = (Mᵢⱼ * vⱼ) .* (np₀^2 * σT^2 * Mpc2m^2)

    return Spline1D(collect(Float64.(ℓrange)), cl)
end

"""
    Cℓτϕ(c; ℓs=DEFAULT_ℓS, zmin=1e-4, zmax=20, kmin=1e-4, kmax=1e1,
         npts=25, limb_fact=0.5, R̄=5, σlnR=log(2), zre=10.0, Δz=2.0, b=6)

Cross-spectrum of the patchy optical depth with the CMB lensing potential,

    Cℓτϕ = (3 H₀² Ωm σT n_p0 / c ℓ²) ∫ (dχ/a³) (χ* - χ)/(χ* χ) Pδe(k = ℓ/χ, χ).

Returns the tuple `(ℓs, Cℓ)` of evaluation multipoles and the cross-spectrum
(already divided by ℓ²).
"""
function Cℓτϕ(c::Cosmology; ℓs=nothing, zmin=1e-4, zmax=20, kmin=1e-4, kmax=1e1,
              npts=25, limb_fact=0.5, R̄=5, σlnR=log(2), zre=10.0, Δz=2.0, b=6)
    np₀ = (1 - 0.75c.Yp) * c.ρc₀ * c.Ωb / mp   # m⁻³
    fact = 3 * (100c.h₀)^2 * c.Ωm * σT * np₀ / c_kms^2

    χmin = comoving_distance(c, zmin)
    χmax = comoving_distance(c, zmax)
    χs = collect(range(χmin, χmax, length=npts))
    zs = c.χtoz(χs)
    dχ = diff(χs); insert!(dχ, 1, 0)

    ℓrange = ℓs === nothing ? DEFAULT_ℓS : ℓs

    kern = [(1 + zs[j])^3 * (c.χlss - χs[j]) / (c.χlss * χs[j]) for j in eachindex(zs)]
    K = [(ℓ + limb_fact) / χ for ℓ in ℓrange, χ in χs]

    Mᵢⱼ = zeros(length(ℓrange), length(χs))   # Pδe(ℓᵢ/χⱼ; χⱼ)
    for j in eachindex(zs)
        Mᵢⱼ[:, j] = [Pδe(c, K[i, j], zs[j]; R̄=R̄, σlnR=σlnR, zre=zre, Δz=Δz, b=b)
                     for i in eachindex(ℓrange)]
    end

    Mᵢⱼ[K .< kmin] .= 0
    Mᵢⱼ[K .≥ kmax] .= 0

    vⱼ = kern .* dχ
    cl = (Mᵢⱼ * vⱼ) .* (fact * Mpc2m)

    ℓout = collect(Float64.(ℓrange))
    return ℓout, cl ./ ℓout .^ 2
end
