using Test

include(joinpath(@__DIR__, "..", "functions.jl"))

@testset "background" begin
    θ = collect(LinRange(deg2rad(5.0), deg2rad(60.0), 1000))

    bg = background(θ)
    @test all(bg .>= 0)
    @test bg[1] > bg[end]

    bg2 = background(θ)
    @test bg == bg2

    @test_throws DomainError background(θ; noise_level=-0.1)
    @test_throws DomainError background(θ; noise_level=1.1)
end
