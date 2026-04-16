using Test

@testset "cubic_multiplicity" begin
    # Canonical h ≥ k ≥ l ≥ 0, standard cubic multiplicities
    @test cubic_multiplicity(1, 0, 0) == 6    # {100}
    @test cubic_multiplicity(1, 1, 0) == 12   # {110}
    @test cubic_multiplicity(1, 1, 1) == 8    # {111}
    @test cubic_multiplicity(2, 1, 0) == 24   # {210}
    @test cubic_multiplicity(2, 1, 1) == 24   # {211}
    @test cubic_multiplicity(2, 2, 1) == 24   # {221}
    @test cubic_multiplicity(3, 2, 1) == 48   # {321} — all distinct, all nonzero
    @test cubic_multiplicity(2, 2, 2) == 8    # {222}
end

@testset "bragg_max_hkl_sq" begin
    # Cu Kα on a=3Å → h²+k²+l² ≤ floor((2·3/1.5418)²) = 15
    @test bragg_max_hkl_sq(3.0, 1.5418) == 15

    # Doubling the lattice parameter quadruples the reachable range
    @test bragg_max_hkl_sq(6.0, 1.5418) == 60

    @test_throws ArgumentError bragg_max_hkl_sq(-1.0, 1.5)
    @test_throws ArgumentError bragg_max_hkl_sq(3.0, -1.0)
end

@testset "Miller_indices" begin
    # SC up to h²+k²+l² ≤ 3 yields {100}, {110}, {111}
    sc_idx, sc_mult = Miller_indices("SC", 3)
    @test sc_idx == [[1,0,0], [1,1,0], [1,1,1]]
    @test sc_mult == [6, 12, 8]

    # BCC up to 4: only h+k+l even — {110}, {200}
    bcc_idx, bcc_mult = Miller_indices("BCC", 4)
    @test bcc_idx == [[1,1,0], [2,0,0]]
    @test bcc_mult == [12, 6]
    @test all(h -> iseven(sum(h)), bcc_idx)

    # FCC up to 4: all-parity — {111}, {200}
    fcc_idx, fcc_mult = Miller_indices("FCC", 4)
    @test fcc_idx == [[1,1,1], [2,0,0]]
    @test fcc_mult == [8, 6]
    @test all(h -> (all(isodd, h) || all(iseven, h)), fcc_idx)

    # Canonical order: h ≥ k ≥ l ≥ 0; no [0,0,0]
    sc_big, _ = Miller_indices("SC", 15)
    @test all(hkl -> hkl[1] ≥ hkl[2] ≥ hkl[3] ≥ 0, sc_big)
    @test [0,0,0] ∉ sc_big

    # Multiplicity sum over canonical set equals count of full [-n,n]³ enumeration
    # (for any cutoff that includes all variants of every family it contains).
    # For SC with max_hkl_sq=3, all variants of {100},{110},{111} fit in [-1,1]³.
    @test sum(sc_mult) == 6 + 12 + 8  # == 26, matching old [-1,1]³ count
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
