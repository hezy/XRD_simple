using Test
using Random

const DATA_TOML = joinpath(@__DIR__, "..", "data.toml")

@testset "sum_peaks" begin
    θ = collect(LinRange(deg2rad(5.0), deg2rad(60.0), 1000))
    θ_list = [0.3, 0.5, 0.8]
    mult = [1, 1, 1]
    w_L = fill(0.01, length(θ))
    w_G = fill(0.005, length(θ))

    result = sum_peaks(θ, θ_list, mult, w_L, w_G)
    @test length(result) == length(θ)
    @test all(result .>= 0)
    @test sum(result) > 0

    single = sum_peaks(θ, [θ_list[1]], [1], w_L, w_G)
    @test sum(result) > sum(single)

    # Doubling the multiplicity doubles the contribution of that peak
    double = sum_peaks(θ, [θ_list[1]], [2], w_L, w_G)
    @test sum(double) ≈ 2 * sum(single)

    @test_throws DimensionMismatch sum_peaks(θ, θ_list, [1, 1], w_L, w_G)
end

@testset "compute_peak_widths" begin
    instrument, peak_width, lattice = read_xrd_config(DATA_TOML)
    θ = collect(LinRange(instrument["two_theta_min"]/2, instrument["two_theta_max"]/2, 100))

    w_L, w_G = compute_peak_widths(θ, peak_width, instrument)
    @test length(w_L) == length(θ)
    @test length(w_G) == length(θ)
    @test all(w_L .> 0)
    @test all(w_G .> 0)
end

@testset "compute_xrd_pattern" begin
    instrument, peak_width, lattice = read_xrd_config(DATA_TOML)
    θ = collect(LinRange(instrument["two_theta_min"]/2, instrument["two_theta_max"]/2, 100))
    λ = instrument["lambda"]
    a = lattice["SC"][2]
    max_hkl_sq = bragg_max_hkl_sq(a, λ)
    indices, multiplicities = Miller_indices("SC", max_hkl_sq)

    w_L, w_G = compute_peak_widths(θ, peak_width, instrument)
    y = compute_xrd_pattern(θ, indices, multiplicities, λ, a, w_L, w_G)
    @test length(y) == length(θ)
    @test all(y .>= 0)
    @test sum(y) > 0

    Random.seed!(42)
    y1 = compute_xrd_pattern(θ, indices, multiplicities, λ, a, w_L, w_G; noise_level=0.1)
    Random.seed!(43)
    y2 = compute_xrd_pattern(θ, indices, multiplicities, λ, a, w_L, w_G; noise_level=0.1)
    @test y1 != y2

    Random.seed!(42)
    y3 = compute_xrd_pattern(θ, indices, multiplicities, λ, a, w_L, w_G; noise_level=0.1)
    @test y1 == y3
end
