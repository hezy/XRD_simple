using Test
using Random

const DATA_TOML = joinpath(@__DIR__, "..", "data.toml")

@testset "sum_peaks" begin
    x = collect(LinRange(deg2rad(10.0), deg2rad(120.0), 1000))
    x_list = [0.6, 1.0, 1.6]
    mult = [1, 1, 1]
    w_L = fill(0.01, 3)
    w_G = fill(0.005, 3)

    result = sum_peaks(x, x_list, mult, w_L, w_G)
    @test length(result) == length(x)
    @test all(result .>= 0)
    @test sum(result) > 0

    single = sum_peaks(x, [x_list[1]], [1], [0.01], [0.005])
    @test sum(result) > sum(single)

    # Doubling the multiplicity doubles the contribution of that peak
    double = sum_peaks(x, [x_list[1]], [2], [0.01], [0.005])
    @test sum(double) ≈ 2 * sum(single)

    # Each peak uses its own widths, evaluated at its centre
    @test sum_peaks(x, [x_list[1]], [1], [0.01], [0.005]) ≈
          pseudo_Voigt_peak(x, x_list[1], 1.0, 0.01, 0.005)

    @test_throws DimensionMismatch sum_peaks(x, x_list, [1, 1], w_L, w_G)
    @test_throws DimensionMismatch sum_peaks(x, x_list, mult, [0.01], w_G)
end

@testset "compute_peak_widths" begin
    _, peak_width, _ = read_xrd_config(DATA_TOML)
    θ_B = deg2rad.([10.0, 30.0, 50.0])

    w_L, w_G = compute_peak_widths(θ_B, peak_width, 1.5418)
    @test length(w_L) == length(θ_B)
    @test length(w_G) == length(θ_B)
    @test all(w_L .> 0)
    @test all(w_G .> 0)
end

@testset "compute_xrd_pattern" begin
    instrument, peak_width, _ = read_xrd_config(DATA_TOML)
    two_θ = collect(LinRange(instrument["two_theta_min"], instrument["two_theta_max"], 2000))
    λ = 1.5418
    a = 3.352
    max_hkl_sq = bragg_max_hkl_sq(a, λ)
    indices, multiplicities = Miller_indices("SC", max_hkl_sq)

    y = compute_xrd_pattern(two_θ, indices, multiplicities, λ, a, peak_width)
    @test length(y) == length(two_θ)
    @test all(y .>= 0)
    @test sum(y) > 0

    # The strongest point of the peaks lies at 2θ_B of (100)
    y_peaks = intensity_vs_angle(two_θ, indices, multiplicities, λ, a, peak_width)
    i = argmin(abs.(two_θ .- 2asin(λ / (2a))))
    @test y_peaks[i] ≈ maximum(y_peaks[max(1, i-20):i+20])

    Random.seed!(42)
    y1 = compute_xrd_pattern(two_θ, indices, multiplicities, λ, a, peak_width; noise_level=0.1)
    Random.seed!(43)
    y2 = compute_xrd_pattern(two_θ, indices, multiplicities, λ, a, peak_width; noise_level=0.1)
    @test y1 != y2

    Random.seed!(42)
    y3 = compute_xrd_pattern(two_θ, indices, multiplicities, λ, a, peak_width; noise_level=0.1)
    @test y1 == y3
end
