using Test

@testset "pseudo_Voigt_peak" begin
    θ = collect(LinRange(0.0, 2.0, 1000))
    θ₀ = 1.0
    A = 1.0
    w_L = 0.01
    w_G = 0.005

    result = pseudo_Voigt_peak(θ, θ₀, A, w_L, w_G)
    @test all(result .>= 0)
    peak_idx = argmax(result)
    @test abs(θ[peak_idx] - θ₀) < 0.01

    result_norm = pseudo_Voigt_peak(θ, θ₀, A, w_L, w_G; normalize=true)
    @test maximum(result_norm) ≈ 1.0 atol=1e-6

    @test result[1] ≈ 0.0 atol=1e-10
    @test result[end] ≈ 0.0 atol=1e-10

    @test_throws ArgumentError pseudo_Voigt_peak(θ, θ₀, -1.0, w_L, w_G)
    @test_throws ArgumentError pseudo_Voigt_peak(θ, θ₀, A, -0.1, w_G)
    @test_throws ArgumentError pseudo_Voigt_peak(θ, θ₀, A, w_L, -0.1)
end

@testset "Voigt_peak" begin
    θ = collect(LinRange(0.0, 2.0, 1000))
    θ₀ = 1.0
    A = 1.0
    w_L = 0.01
    w_G = 0.005

    result = Voigt_peak(θ, θ₀, A, w_L, w_G)
    @test all(result .>= 0)
    peak_idx = argmax(result)
    @test abs(θ[peak_idx] - θ₀) < 0.01

    result_norm = Voigt_peak(θ, θ₀, A, w_L, w_G; normalize=true)
    @test maximum(result_norm) ≈ 1.0 atol=1e-6

    @test result[1] ≈ 0.0 atol=1e-10
    @test result[end] ≈ 0.0 atol=1e-10
end

@testset "peak_fwhm" begin
    @test peak_fwhm(0.01, 0.005) ≈ 0.5346*0.01 + sqrt(0.2166*0.01^2 + 0.005^2) atol=1e-10
    @test peak_fwhm(0.0, 0.01) ≈ 0.01 atol=1e-10
    @test peak_fwhm(0.01, 0.0) ≈ 0.5346*0.01 + sqrt(0.2166)*0.01 atol=1e-10

    w_L_vec = [0.01, 0.02]
    w_G_vec = [0.005, 0.01]
    result = peak_fwhm(w_L_vec, w_G_vec)
    @test length(result) == 2
    @test result[1] ≈ peak_fwhm(0.01, 0.005) atol=1e-10
    @test result[2] ≈ peak_fwhm(0.02, 0.01) atol=1e-10

    @test_throws DimensionMismatch peak_fwhm([0.01], [0.005, 0.01])
end
