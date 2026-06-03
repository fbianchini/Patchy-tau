# 3-D power spectra of the free-electron (ionization) field and its
# cross-spectrum with the matter density.
#
# Notation follows Wang & Hu (2006), astro-ph/0607652:
#   Pee  = electron-electron power spectrum  (1-halo + 2-halo)
#   Pδe  = density-electron cross spectrum
# All spectra are in Mpc³ and take `k` in Mpc⁻¹.

"""
    P̃δδ(c, k, z; R̄=5, σlnR=log(2), σ²b=nothing)

Interpolating term between the small- and large-scale limits of the 1-halo
electron spectrum (Wang & Hu 2006, Eq. 17).

`σ²b` is the bubble-averaged density variance ⟨σ²_R⟩(z); it is `k`-independent,
so when looping over many `k` at fixed `z` it can be precomputed once with
[`sigma2_bubble`](@ref) and passed in to avoid recomputation. If `nothing`
(default) it is computed internally.
"""
function P̃δδ(c::Cosmology, k::Real, z::Real; R̄=5, σlnR=log(2), σ²b=nothing)
    σ²b === nothing && (σ²b = sigma2_bubble(c, z; R̄=R̄, σlnR=σlnR))
    pkz = matter_power(c, k, z)
    vb = mean_bubble_volume(R̄, σlnR)
    num = pkz * σ²b * vb
    den = √(pkz^2 + (σ²b * vb)^2)
    num / den
end

"""
    Pee1h(c, k, z; R̄=5, σlnR=log(2), zre=7.0, Δz=1.0, σ²b=nothing)

One-halo (shot-noise + Poisson) contribution to the electron power spectrum
(Wang & Hu 2006, Eq. 17). See [`P̃δδ`](@ref) for the precomputed `σ²b` argument.
"""
function Pee1h(c::Cosmology, k::Real, z::Real; R̄=5, σlnR=log(2), zre=7.0, Δz=1.0, σ²b=nothing)
    xe = ionized_fraction(z; zre=zre, Δz=Δz)
    xe * (1 - xe) * (filter_W2_quad(k; R̄=R̄, σlnR=σlnR) * mean_bubble_volume(R̄, σlnR) +
                     P̃δδ(c, k, z; R̄=R̄, σlnR=σlnR, σ²b=σ²b))
end

"""
    Pee2h(c, k, z; R̄=5, σlnR=log(2), zre=7.0, Δz=1.0, b=6.0)

Two-halo (clustering) contribution to the electron power spectrum
(Wang & Hu 2006, Eq. 31), with linear bias `b`.
"""
function Pee2h(c::Cosmology, k::Real, z::Real; R̄=5, σlnR=log(2), zre=7.0, Δz=1.0, b=6.0)
    xe = ionized_fraction(z; zre=zre, Δz=Δz)
    xH = 1 - xe
    W = filter_W_quad(k; R̄=R̄, σlnR=σlnR)
    (xH * log(xH) * b * W - xe)^2 * matter_power(c, k, z)
end

"""
    Pee(c, k, z; R̄=5, σlnR=log(2), zre=7.0, Δz=1.0, b=6.0, σ²b=nothing)

Total electron-electron power spectrum, `Pee1h + Pee2h`. See [`P̃δδ`](@ref) for
the precomputed `σ²b` argument.
"""
function Pee(c::Cosmology, k::Real, z::Real; R̄=5, σlnR=log(2), zre=7.0, Δz=1.0, b=6.0, σ²b=nothing)
    Pee1h(c, k, z; R̄=R̄, σlnR=σlnR, zre=zre, Δz=Δz, σ²b=σ²b) +
    Pee2h(c, k, z; R̄=R̄, σlnR=σlnR, zre=zre, Δz=Δz, b=b)
end

"""
    Pδe(c, k, z; R̄=5, σlnR=log(2), zre=7.0, Δz=1.0, b=6.0)

Density-electron cross power spectrum.
"""
function Pδe(c::Cosmology, k::Real, z::Real; R̄=5, σlnR=log(2), zre=7.0, Δz=1.0, b=6.0)
    xe = ionized_fraction(z; zre=zre, Δz=Δz)
    xH = 1 - xe
    W = filter_W_quad(k; R̄=R̄, σlnR=σlnR)
    (xe - xH * log(xH) * b * W) * matter_power(c, k, z)
end
