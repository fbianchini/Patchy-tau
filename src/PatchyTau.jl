"""
    PatchyTau

Theory power spectra for patchy reionization: the optical-depth auto-spectrum
`Cℓττ`, its cross-spectrum with CMB lensing `Cℓτϕ`, and the underlying 3-D
free-electron power spectra, built on an ionized-bubble halo model.

Quick start:

```julia
using PatchyTau
c = Cosmology()                       # fiducial parameters
optical_depth(c, 30; zre=10, Δz=4)    # τ̄ out to z = 30
clττ = Cℓττ(c; R̄=5, σlnR=log(2), zre=10, Δz=4)
clττ(1000)                            # evaluate at ℓ = 1000
ℓ, clτϕ = Cℓτϕ(c; R̄=5, σlnR=log(2), zre=10, Δz=4)
```
"""
module PatchyTau

using QuadGK
using Trapz
using Dierckx
import MatterPower

include("constants.jl")
include("cosmology.jl")
include("reionization.jl")
include("bubbles.jl")
include("electron_spectra.jl")
include("angular_spectra.jl")

# Cosmology / background
export Cosmology, hubble, comoving_distance, growth_factor, matter_power

# Reionization
export ionized_fraction, optical_depth

# Bubble model
export bubble_pdf, bubble_volume, mean_bubble_volume, mean_bubble_volume2,
       window, filter_W, filter_W_quad, filter_W2, filter_W2_quad,
       mass_threshold, M2R, sigma2, sigma2_bubble, sigma2_bubble_meerburg

# 3-D electron spectra
export Pee, Pee1h, Pee2h, Pδe, P̃δδ

# Angular spectra
export Wκ, lensing_kernel, Cℓκκ, Cℓττ, Cℓτϕ, DEFAULT_ℓS

end # module
