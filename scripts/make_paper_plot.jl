#!/usr/bin/env julia
#
# Reproduce the three-panel paper figure from the notebook:
#   (1) optical-depth history τ̄(z)
#   (2) L² Cℓττ for two bubble models
#   (3) L³ Cℓτϕ for two bubble models
#
# Requires PyPlot + CMBLensing (see `scripts/setup.jl --with-plots`):
#
#     julia --project=. scripts/make_paper_plot.jl

using PatchyTau
using PyPlot
using CMBLensing: extrapolate_Cℓs

c = Cosmology()

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 15
rcParams["mathtext.fontset"] = "stix"
rcParams["font.family"] = "STIXGeneral"
rcParams["axes.linewidth"] = 1.2
rcParams["axes.labelsize"] = 15
rcParams["axes.titlesize"] = 18
rcParams["xtick.labelsize"] = 10
rcParams["ytick.labelsize"] = 10
rcParams["xtick.major.size"] = 4
rcParams["ytick.major.size"] = 4
rcParams["xtick.minor.size"] = 3
rcParams["ytick.minor.size"] = 3
rcParams["legend.fontsize"] = 14
rcParams["legend.frameon"] = "False"
rcParams["xtick.major.width"] = 1
rcParams["ytick.major.width"] = 1
rcParams["xtick.minor.width"] = 1
rcParams["ytick.minor.width"] = 1

figure(figsize=(12, 4))

# --- Panel 1: optical-depth history ---
subplot(131)
zz = 0:0.01:20
plot(zz, [optical_depth(c, z; zre=10, Δz=4) for z in zz], label=L"(z_{\rm re}, \Delta z)=(10,4)")
plot(zz, [optical_depth(c, z; zre=7, Δz=1) for z in zz], label=L"(z_{\rm re}, \Delta z)=(7,1)")
axhspan(0.058 - 0.012, 0.058 + 0.012, alpha=0.2, label="Planck18")
axhspan(0.058 - 2 * 0.012, 0.058 + 2 * 0.012, alpha=0.1)
xlabel(L"z", size=20)
ylabel(L"\bar{\tau}", size=20)
legend()

# --- Panel 2: optical-depth auto-spectrum ---
subplot(132)
loglog(extrapolate_Cℓs(1:3000, 1:3000, collect(1:3000) .^ 2 .* Cℓττ(c; R̄=5, σlnR=log(2), zre=10, Δz=4)(1:3000)),
       label=L"(\bar{R},\sigma_{\ln R})=(5,\ln(2))")
loglog(extrapolate_Cℓs(1:3000, 1:3000, collect(1:3000) .^ 2 .* Cℓττ(c; R̄=1, σlnR=0.1, zre=7, Δz=1)(1:3000)),
       label=L"(\bar{R},\sigma_{\ln R})=(1,0.1)")
xlabel("Multipole L", size=20)
ylabel(L"L^2C_L^{ττ}", size=20)
xlim(10, 3000)
ylim(1e-9, 3e-4)
grid()
legend()

# --- Panel 3: optical-depth x lensing cross-spectrum ---
subplot(133)
ℓ, clτϕ_opt = Cℓτϕ(c; R̄=5, σlnR=log(2), zre=10, Δz=4)
_, clτϕ_pes = Cℓτϕ(c; R̄=1, σlnR=0.1, zre=7, Δz=1)
plot(extrapolate_Cℓs(2:3000, ℓ, ℓ .^ 3 .* clτϕ_opt), label=L"(\bar{R},\sigma_{\ln R})=(5,\ln(2))")
plot(extrapolate_Cℓs(2:3000, ℓ, ℓ .^ 3 .* clτϕ_pes), label=L"(\bar{R},\sigma_{\ln R})=(1,0.1)")
xlabel("Multipole L", size=20)
ylabel(L"L^3C_L^{τϕ}", size=20)
xlim(10, 3000)
ylim(0, 1.75e-6)
grid()
legend()

tight_layout()
savefig(joinpath(@__DIR__, "..", "theory_spectra_paper.pdf"))
@info "Wrote theory_spectra_paper.pdf"
