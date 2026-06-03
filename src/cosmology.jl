# Background cosmology: parameters, expansion history, linear growth and the
# (no-wiggle) linear matter power spectrum.
#
# A `Cosmology` bundles the input parameters together with everything that is
# expensive to compute once and reused many times: the growth-factor ODE
# solution, the χ <-> z lookup spline and a 2-D spline of the matter power
# spectrum P(k, z). Construct one with the keyword constructor and pass it as
# the first argument to the physics routines.

"""
    Cosmology(; kwargs...)

Container for the background cosmology plus precomputed lookup tables.

Keyword arguments (defaults reproduce the fiducial model of the original
notebook):

- `As`, `ns`, `kpivot` : primordial amplitude, tilt and pivot scale [Mpc⁻¹]
- `Ωm`, `Ωb`           : matter and baryon density parameters
- `h₀`                 : dimensionless Hubble parameter (H₀ = 100 h₀ km/s/Mpc)
- `Yp`                 : primordial helium mass fraction
- `ΩΛ`                 : dark-energy density (defaults to a flat universe, 1 - Ωm)

The grids used to build the matter-power spline can be tuned via `kmin`, `kmax`,
`dlogk` and the number/edges of the redshift grid (`zgrid`).
"""
struct Cosmology{G,S1,S2}
    # --- input parameters ---
    As::Float64
    ns::Float64
    kpivot::Float64
    Ωm::Float64
    ΩΛ::Float64
    Ωb::Float64
    Ωk::Float64
    h₀::Float64
    Yp::Float64
    # --- derived quantities ---
    ωm::Float64          # physical matter density Ωm h²
    fb::Float64          # baryon fraction Ωb / Ωm
    ρc₀::Float64         # critical density today [kg/m³]
    χlss::Float64        # comoving distance to last-scattering [Mpc]
    # --- precomputed lookups ---
    growth::G            # linear growth ODE solution, growth(a)[1] = D(a)
    χtoz::S1             # spline mapping comoving distance [Mpc] -> redshift
    ks::Vector{Float64}  # wavenumber grid [Mpc⁻¹] backing pkz_spl
    zs::Vector{Float64}  # redshift grid backing pkz_spl
    pkz_spl::S2          # 2-D spline P(k, z) [Mpc³]
end

function Cosmology(; As=2.46e-9, ns=0.96, kpivot=0.002,
                     Ωm=0.23 + 0.044, Ωb=0.044, h₀=0.704, Yp=0.24,
                     ΩΛ=nothing,
                     zgrid=vcat(range(0, 10, length=100), range(10.1, 1100, length=100)),
                     kmin=1e-4, kmax=100.0, dlogk=0.01)

    ΩΛ === nothing && (ΩΛ = 1 - Ωm)
    Ωk = 1 - Ωm - ΩΛ
    ωm = Ωm * h₀^2
    fb = Ωb / Ωm
    ρc₀ = 3 * (100h₀ * km2Mpc)^2 / (8π * G)   # kg/m³

    # Local closures over the parameters so we can build the lookups before the
    # struct exists.
    Hloc(z) = 100h₀ * √(Ωm * (1 + z)^3 + ΩΛ)                                  # km/s/Mpc
    χloc(z) = QuadGK.quadgk(zed -> 1 / Hloc(zed), 0, z, rtol=1e-3)[1] * c_kms  # Mpc
    χlss = χloc(1089)

    zg = collect(zgrid)
    χtoz = Spline1D(χloc.(zg), zg, k=1)

    growth = MatterPower.setup_growth(Ωm, ΩΛ)

    pkloc(k, red) = begin
        kovh = k / h₀
        growth(1 / (1 + red))[1]^2 *
        As * (k / kpivot)^(ns - 1) *
        (2 * kovh^2 * 2998^2 / 5 / Ωm)^2 *
        MatterPower.t_nowiggle(k, ωm, fb)^2 *
        2π^2 / k^3
    end

    ks = exp.(log(kmin):dlogk:log(kmax))
    zs = zg
    pkz = zeros(length(ks), length(zs))
    for i in eachindex(zs)
        pkz[:, i] = pkloc.(ks, zs[i])
    end
    pkz_spl = Spline2D(ks, zs, pkz)

    return Cosmology(As, ns, kpivot, Ωm, ΩΛ, Ωb, Ωk, h₀, Yp,
                     ωm, fb, ρc₀, χlss, growth, χtoz, collect(ks), zs, pkz_spl)
end

"""
    hubble(c, z)

Hubble parameter H(z) in km/s/Mpc.
"""
hubble(c::Cosmology, z::Real) = 100c.h₀ * √(c.Ωm * (1 + z)^3 + c.ΩΛ)

"""
    comoving_distance(c, z)

Line-of-sight comoving distance χ(z) in Mpc.
"""
comoving_distance(c::Cosmology, z::Real) =
    QuadGK.quadgk(zed -> 1 / hubble(c, zed), 0, z, rtol=1e-3)[1] * c_kms

"""
    growth_factor(c, z)

Linear growth factor D(z) (un-normalised, as returned by `MatterPower`).
"""
growth_factor(c::Cosmology, z::Real) = c.growth(1 / (1 + z))[1]

"""
    matter_power(c, k, z)

Linear (no-wiggle) matter power spectrum P(k, z) in Mpc³ evaluated directly
(not from the spline). `k` is in Mpc⁻¹.
"""
function matter_power(c::Cosmology, k::Real, z::Real)
    kovh = k / c.h₀
    growth_factor(c, z)^2 *
    c.As * (k / c.kpivot)^(c.ns - 1) *
    (2 * kovh^2 * 2998^2 / 5 / c.Ωm)^2 *
    MatterPower.t_nowiggle(k, c.ωm, c.fb)^2 *
    2π^2 / k^3
end
