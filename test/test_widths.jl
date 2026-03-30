using Test

include(joinpath(@__DIR__, "..", "functions.jl"))

@testset "Gaussian_peaks_width" begin
    θ = collect(LinRange(deg2rad(5.0), deg2rad(60.0), 100))
    U, V, W = 0.0001, -0.00005, 0.00001

    result = Gaussian_peaks_width(θ, U, V, W)
    @test length(result) == length(θ)
    @test all(result .> 0)
    @test result[end] > result[1]
end

@testset "Lorentzian_peaks_width" begin
    θ = collect(LinRange(deg2rad(5.0), deg2rad(60.0), 100))
    K, ϵ, λ, D = 0.9, 0.001, 1.5418, 500.0

    result = Lorentzian_peaks_width(θ, K, ϵ, λ, D)
    @test length(result) == length(θ)
    @test all(result .> 0)
    @test result[end] > result[1]
end
