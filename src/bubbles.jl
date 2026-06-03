# Ionized-bubble model: a log-normal distribution of bubble radii, the
# associated real- and Fourier-space filters, and the variance of the density
# field smoothed on bubble scales.
#
# Reference: Wang & Hu (2006), astro-ph/0607652.

"""
    bubble_pdf(R; R̄=5.0, σlnR=log(2))

Log-normal probability distribution P(R) of ionized-bubble radii (`R` in Mpc),
with median radius `R̄` and log-width `σlnR`.
"""
bubble_pdf(R::Real; R̄=5.0, σlnR=log(2)) =
    exp(-((log(R / R̄)^2) / (2σlnR^2))) / sqrt(2π * σlnR^2) / R

"""
    bubble_volume(R)

Volume of a spherical bubble of radius `R`.
"""
bubble_volume(R::Real) = 4π * R^3 / 3

"""
    mean_bubble_volume(R̄=5.0, σlnR=log(2))

First moment ⟨V_b⟩ of the bubble volume under the log-normal distribution.
"""
mean_bubble_volume(R̄=5.0, σlnR=log(2)) = 4π * R̄^3 * exp(4.5σlnR^2) / 3

"""
    mean_bubble_volume2(R̄=5.0, σlnR=log(2))

Second moment ⟨V_b²⟩ of the bubble volume under the log-normal distribution.
"""
mean_bubble_volume2(R̄=5.0, σlnR=log(2)) = (4π)^2 * R̄^6 * exp(18σlnR^2) / 9

"""
    window(k, R)

Fourier transform of a real-space spherical top-hat of radius `R` at
wavenumber `k`.
"""
window(k::Real, R::Real) = 3 * (k * R)^(-3) * (sin(k * R) - (k * R) * cos(k * R))

"""
    filter_W(k; R̄=5.0, σlnR=log(2), Rmax=50, npts=300)

Volume-weighted mean of the top-hat window over the bubble distribution,
⟨W_R⟩(k), via trapezoidal integration in R.
"""
function filter_W(k::Real; R̄=5.0, σlnR=log(2), Rmax=50, npts=300)
    R = range(1e-3, Rmax, length=npts)
    dWᵣdR = [bubble_pdf(r; R̄=R̄, σlnR=σlnR) * bubble_volume(r) * window(k, r) for r in R]
    trapz(R, dWᵣdR) / mean_bubble_volume(R̄, σlnR)
end

"""
    filter_W_quad(k; R̄=5.0, σlnR=log(2), χmax=Inf)

Same as [`filter_W`](@ref) but evaluated with adaptive quadrature in the
dimensionless variable χ = R / R̄.
"""
function filter_W_quad(k::Real; R̄=5.0, σlnR=log(2), χmax=Inf)
    κ = k * R̄
    dWᵣdχ(χ) = χ^2 * exp(-(log(χ)^2) / (2σlnR^2)) * (sin(κ * χ) - κ * χ * cos(κ * χ)) / (κ * χ)^3
    QuadGK.quadgk(dWᵣdχ, 0, χmax, rtol=1e-2)[1] * (3exp(-9σlnR^2 / 2) / (√(2π) * σlnR))
end

"""
    filter_W2(k; R̄=5.0, σlnR=log(2), Rmax=300, npts=100)

Volume-squared-weighted mean of the squared window over the bubble
distribution, ⟨W²_R⟩(k), via trapezoidal integration in R.
"""
function filter_W2(k::Real; R̄=5.0, σlnR=log(2), Rmax=300, npts=100)
    R = range(1e-3, Rmax, length=npts)
    dW²ᵣdR = [bubble_pdf(r; R̄=R̄, σlnR=σlnR) * bubble_volume(r)^2 * window(k, r)^2 for r in R]
    trapz(R, dW²ᵣdR) / mean_bubble_volume(R̄, σlnR)^2
end

"""
    filter_W2_quad(k; R̄=5.0, σlnR=log(2), χmax=Inf)

Same as [`filter_W2`](@ref) but evaluated with adaptive quadrature in the
dimensionless variable χ = R / R̄.
"""
function filter_W2_quad(k::Real; R̄=5.0, σlnR=log(2), χmax=Inf)
    κ = k * R̄
    dW²ᵣdχ(χ) = χ^5 * exp(-(log(χ)^2) / (2σlnR^2)) * (sin(κ * χ) - κ * χ * cos(κ * χ))^2 / (κ * χ)^6
    QuadGK.quadgk(dW²ᵣdχ, 0, χmax, rtol=1e-2)[1] * (9exp(-9σlnR^2) / (√(2π) * σlnR))
end

"""
    mass_threshold(z; ωh²=0.143, Tvir=5e4)

Minimum halo mass (in M⊙) able to collapse and host ionizing sources at
redshift `z`, for a virial temperature `Tvir` [K].
"""
mass_threshold(z::Real; ωh²=0.143, Tvir=5e4) =
    (10 / (1 + z))^1.5 * (0.15 / ωh²)^2 * (Tvir / 1e4)^1.5 * (1 / 1.1)^1.5 * 1e8

"""
    M2R(M, ρ̄₀)

Lagrangian radius enclosing mass `M` at mean density `ρ̄₀`.
"""
M2R(M::Real, ρ̄₀::Real) = ∛((3M) / (4π * ρ̄₀))

"""
    sigma2(c, z, R; npts=50)

Variance σ²(R, z) of the linear density field smoothed with a top-hat of
radius `R` [Mpc], by trapezoidal integration of the matter power spectrum.
"""
function sigma2(c::Cosmology, z::Real, R::Real; npts=50)
    ks = range(1e-4, 20 / R, length=npts)
    integ = map(ks) do k
        x = k * R
        W = (3 / x) * (sin(x) / x^2 - cos(x) / x)
        W^2 * matter_power(c, k, z) * k^2 / 2 / π^2
    end
    trapz(ks, integ)
end

"""
    sigma2_bubble(c, z; R̄=5.0, σlnR=log(2), Rmax=300, npts=50)

Bubble-distribution-averaged density variance ⟨σ²_R⟩(z), weighted by V_b².
"""
function sigma2_bubble(c::Cosmology, z::Real; R̄=5.0, σlnR=log(2), Rmax=300, npts=50)
    R = range(1e-3, Rmax, length=npts)
    dsigma2 = [bubble_pdf(r; R̄=R̄, σlnR=σlnR) * bubble_volume(r)^2 * sigma2(c, z, r) for r in R]
    trapz(R, dsigma2) / mean_bubble_volume(R̄, σlnR)^2
end

"""
    sigma2_bubble_meerburg(c, z; R̄=5.0, σlnR=log(2), npts=50)

Alternative estimate of ⟨σ²_R⟩(z) obtained by integrating P(k) against the
Fourier-space bubble filter ⟨W²_R⟩(k).
"""
function sigma2_bubble_meerburg(c::Cosmology, z::Real; R̄=5.0, σlnR=log(2), npts=50)
    ks = range(1e-4, 20, length=npts)
    integ = [filter_W2_quad(k; R̄=R̄, σlnR=σlnR) * matter_power(c, k, z) * k^2 / 2 / π^2 for k in ks]
    trapz(ks, integ)
end
