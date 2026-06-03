using PatchyTau
using Test

@testset "PatchyTau" begin
    c = Cosmology()

    @testset "background" begin
        @test hubble(c, 0) ≈ 100c.h₀
        @test comoving_distance(c, 0) ≈ 0 atol = 1e-6
        @test c.χlss > 13000           # ~14 Gpc to last scattering
        @test matter_power(c, 0.1, 0) > 0
    end

    @testset "reionization" begin
        # Monotonically decreasing with redshift, bounded in [0, 1].
        @test ionized_fraction(0; zre=10, Δz=4) ≈ 1 atol = 1e-3
        @test 0 ≤ ionized_fraction(10; zre=10, Δz=4) ≤ 1
        @test ionized_fraction(20; zre=10, Δz=4) < ionized_fraction(10; zre=10, Δz=4)
        # Optical depth grows with z and is of order the Planck value by z ~ 30.
        @test optical_depth(c, 5; zre=10, Δz=4) < optical_depth(c, 30; zre=10, Δz=4)
        @test 0.02 < optical_depth(c, 30; zre=10, Δz=4) < 0.2
    end

    @testset "bubble model" begin
        @test bubble_volume(2) ≈ 4π * 8 / 3
        @test mean_bubble_volume(5, log(2)) > 0
        @test sigma2(c, 0, 8) > 0
        @test filter_W(1.0; R̄=5) isa Real
        @test filter_W2_quad(1.0; R̄=5) ≥ 0
    end

    @testset "electron spectra" begin
        @test Pee(c, 0.1, 9; R̄=5, σlnR=log(2), zre=10, Δz=4) > 0
        @test Pee1h(c, 0.1, 9; R̄=5) ≥ 0
        @test Pδe(c, 0.1, 9; R̄=5, zre=10, Δz=4) isa Real
    end

    @testset "angular spectra" begin
        clκκ = Cℓκκ(c, lensing_kernel(c); ℓmax=200, npts=20)
        @test length(clκκ) == 199
        @test all(clκκ .≥ 0)

        clττ = Cℓττ(c; R̄=5, σlnR=log(2), zre=10, Δz=4, npts=15)
        @test clττ(1000) > 0

        ℓ, clτϕ = Cℓτϕ(c; R̄=5, σlnR=log(2), zre=10, Δz=4, npts=15)
        @test length(ℓ) == length(clτϕ)
    end
end
