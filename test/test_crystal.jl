using Test

include(joinpath(@__DIR__, "..", "functions.jl"))

@testset "Miller_indices" begin
    sc = Miller_indices("SC", -1, 1)
    @test length(sc) == 26

    bcc = Miller_indices("BCC", -1, 1)
    @test all(h -> iseven(sum(h)), bcc)

    fcc = Miller_indices("FCC", -1, 1)
    @test all(h -> (all(isodd, h) || all(iseven, h)), fcc)

    @test [0,0,0] ∉ sc
    @test [0,0,0] ∉ bcc
    @test [0,0,0] ∉ fcc
end

@testset "d_list" begin
    a_cu = 3.594
    @test d_list([[1,1,1]], a_cu)[1] ≈ 2.075 atol=1e-3
    @test d_list([[2,0,0]], a_cu)[1] ≈ 1.797 atol=1e-3
    @test d_list([[2,2,0]], a_cu)[1] ≈ 1.271 atol=1e-3
    @test d_list([[3,1,1]], a_cu)[1] ≈ 1.084 atol=1e-3
    @test d_list([[2,2,2]], a_cu)[1] ≈ 1.038 atol=1e-3

    a_fe = 2.866
    @test d_list([[1,1,0]], a_fe)[1] ≈ 2.027 atol=1e-3
    @test d_list([[2,0,0]], a_fe)[1] ≈ 1.433 atol=1e-3
    @test d_list([[2,1,1]], a_fe)[1] ≈ 1.170 atol=1e-3
    @test d_list([[2,2,0]], a_fe)[1] ≈ 1.013 atol=1e-3

    @test_throws DomainError d_list([[1,1,1]], -1.0)
end

@testset "bragg_angles" begin
    λ = 1.5418

    d_spacings = [2.075, 1.797, 1.271]
    angles, idx = bragg_angles(λ, d_spacings)
    @test length(angles) == 3
    @test all(a -> 0 < a < π/2, angles)
    @test idx == [1, 2, 3]

    d_impossible = [0.5, 0.7]
    angles2, idx2 = bragg_angles(λ, d_impossible)
    @test length(angles2) == 0

    d_mixed = [2.075, 0.5, 1.797]
    angles3, idx3 = bragg_angles(λ, d_mixed)
    @test length(angles3) == 2
    @test idx3 == [1, 3]

    @test_throws ArgumentError bragg_angles(-1.0, [1.0])
end
