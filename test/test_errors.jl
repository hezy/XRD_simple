using Test

include(joinpath(@__DIR__, "..", "functions.jl"))

@testset "Miller_indices errors" begin
    @test_throws ArgumentError Miller_indices("HCP", -1, 1)
    @test_throws ArgumentError Miller_indices("SC", 5, 2)
end

@testset "bragg_angles errors" begin
    @test_throws ArgumentError bragg_angles(-1.0, [1.0])
    @test_throws ArgumentError bragg_angles(1.0, [-1.0])
end

@testset "d_list errors" begin
    @test_throws DomainError d_list([[1,1,1]], -1.0)
    @test_throws DomainError d_list([[1,1,1]], 0.0)
end

@testset "peak function errors" begin
    θ = collect(LinRange(0.0, 2.0, 100))
    @test_throws ArgumentError pseudo_Voigt_peak(θ, 1.0, -1.0, 0.01, 0.005)
    @test_throws ArgumentError Voigt_peak(θ, 1.0, 1.0, -0.01, 0.005)
end
