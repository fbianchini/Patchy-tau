# Reionization history: the mean ionized fraction and the homogeneous optical
# depth that follows from it.

"""
    ionized_fraction(z; zre=7.0, Δz=1.0, f=1.0)

Mean ionized fraction x̄ₑ(z) for a `tanh` reionization history with midpoint
redshift `zre`, width `Δz` and asymptotic ionized fraction `f`.
"""
function ionized_fraction(z::Real; zre=7.0, Δz=1.0, f=1.0)
    y(x) = (1 + x)^1.5
    yre = y(zre)
    Δy = 1.5 * sqrt(1 + zre) * Δz
    f / 2 * (1 + tanh((yre - y(z)) / Δy))
end

"""
    optical_depth(c, z; zre=7.0, Δz=1.0, f=1.0, npts=100)

Homogeneous Thomson optical depth τ̄(z) integrated from 0 to `z`, for the
reionization history set by `zre`, `Δz`, `f`. Uses trapezoidal integration on
`npts` points.
"""
function optical_depth(c::Cosmology, z::Real; zre=7.0, Δz=1.0, f=1.0, npts=100)
    ρb₀ = c.ρc₀ * c.Ωb                       # kg/m³
    fact = σT * (1 - 0.75c.Yp) * ρb₀ / mp    # m⁻¹
    zarr = range(0, z, length=npts)
    dτdz = [ionized_fraction(zz; zre=zre, Δz=Δz, f=f) * (1 + zz)^2 / hubble(c, zz)
            for zz in zarr]
    trapz(zarr, dτdz) * fact * c_kms * Mpc2m
end
