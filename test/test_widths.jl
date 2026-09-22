using Test

@testset "Gaussian_peaks_width" begin
    θ = collect(LinRange(deg2rad(5.0), deg2rad(60.0), 100))
    U, V, W = 0.0001, -0.00005, 0.00001

    result = Gaussian_peaks_width.(θ, U, V, W)
    @test length(result) == length(θ)
    @test all(result .> 0)
    @test result[end] > result[1]
end

@testset "Lorentzian_peaks_width" begin
    θ = collect(LinRange(deg2rad(5.0), deg2rad(60.0), 100))
    K, ϵ, λ, D = 0.9, 0.001, 1.5418, 500.0

    result = Lorentzian_peaks_width.(θ, K, ϵ, λ, D)
    @test length(result) == length(θ)
    @test all(result .> 0)
    @test result[end] > result[1]
end

@testset "Scherrer size broadening" begin
    # Strain off: FWHM (radians of 2θ) = Kλ/(D cos θ), with D converted nm → Å
    K, λ, D_nm = 0.9, 1.5418, 50.0
    θ = deg2rad(20.0)
    expected = 0.9 * 1.5418 / (500.0 * cos(deg2rad(20.0)))   # ≈ 2.953e-3 rad
    @test Lorentzian_peaks_width(θ, K, 0.0, λ, D_nm) ≈ expected rtol=1e-12
    @test Lorentzian_peaks_width(θ, K, 0.0, λ, D_nm) ≈ 2.953e-3 rtol=1e-3
end

@testset "Stokes-Wilson strain broadening" begin
    # Size term negligible (huge D): FWHM (radians of 2θ) = 4ε tan θ
    θ = deg2rad(30.0)
    @test Lorentzian_peaks_width(θ, 0.9, 0.002, 1.5418, 1e12) ≈ 4 * 0.002 * tan(θ) rtol=1e-6
end

@testset "Caglioti gives FWHM" begin
    θ = deg2rad(25.0)
    U, V, W = 1e-4, -2e-5, 1e-5
    @test Gaussian_peaks_width(θ, U, V, W) ≈ sqrt(U * tan(θ)^2 + V * tan(θ) + W)
end
