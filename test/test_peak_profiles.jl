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
end

# Numerical FWHM of a sampled single peak, by linear interpolation of the
# half-maximum crossings on each side of the maximum
function measured_fwhm(x, y)
    i₀ = argmax(y)
    half = y[i₀] / 2
    i = i₀
    while y[i] > half; i -= 1; end
    x_left = x[i] + (half - y[i]) / (y[i+1] - y[i]) * (x[i+1] - x[i])
    j = i₀
    while y[j] > half; j += 1; end
    x_right = x[j-1] + (half - y[j-1]) / (y[j] - y[j-1]) * (x[j] - x[j-1])
    return x_right - x_left
end

@testset "profile FWHM matches peak_fwhm" begin
    x = collect(LinRange(0.0, 2.0, 200_001))
    x₀ = 1.0
    for (w_L, w_G) in [(0.01, 0.005), (0.005, 0.01), (0.01, 0.01), (0.002, 0.02)]
        f = peak_fwhm(w_L, w_G)
        @test measured_fwhm(x, Voigt_peak(x, x₀, 1.0, w_L, w_G)) ≈ f rtol=0.01
        @test measured_fwhm(x, pseudo_Voigt_peak(x, x₀, 1.0, w_L, w_G)) ≈ f rtol=0.01
    end
end

@testset "profile symmetry" begin
    # Odd number of points, centre on the middle point
    x = collect(LinRange(0.0, 2.0, 2001))
    x₀ = 1.0
    for peak in (Voigt_peak, pseudo_Voigt_peak)
        y = peak(x, x₀, 1.0, 0.01, 0.005)
        @test y ≈ reverse(y) rtol=1e-10
    end
end

@testset "profile area" begin
    # A large cutoff keeps the Lorentzian tails, which the default 5·FWHM cuts
    x = collect(LinRange(0.0, 2.0, 200_001))
    dx = x[2] - x[1]
    x₀, A = 1.0, 3.0
    for peak in (Voigt_peak, pseudo_Voigt_peak)
        y = peak(x, x₀, A, 0.005, 0.005; cutoff_sigma=100.0)
        @test sum(y) * dx ≈ A rtol=0.01
    end
end
