using Test

const SF_XRAY_TOML = joinpath(@__DIR__, "reference", "xray.toml")
const SF_ELECTRON_TOML = joinpath(@__DIR__, "reference", "electron.toml")

# The configuration of `file` with the reflection method `model` and no noise
function with_model(file, model)
    c = read_xrd_config(file)
    return XRDConfig(c.mode, model, c.N, 0.0, c.K, c.Epsilon, c.D, c.samples)
end

origin(element, B=0.0) = Atom(element, (0.0, 0.0, 0.0), B)
cell(name, centering, a, atoms...) = Sample(name, centering, a, collect(atoms))

# Scattering weights of the families `hkls` at s (X-ray form factors)
function sf_weights(sample, hkls, s)
    mode = read_xrd_config(SF_XRAY_TOML).mode
    return scattering_weights(StructureFactor(), mode, sample, hkls, fill(s, length(hkls)))
end

has(indices, hkl) = hkl in indices

@testset "family_members" begin
    for hkl in ([1, 0, 0], [1, 1, 0], [1, 1, 1], [2, 1, 0], [2, 2, 1], [3, 2, 1], [4, 0, 0])
        members = family_members(hkl)
        @test length(members) == cubic_multiplicity(hkl...)
        @test all(sum(abs2, m) == sum(abs2, hkl) for m in members)
    end
end

@testset "unit_cell" begin
    @test length(unit_cell(lattice_sample("SC", "Po", 3.352))) == 1
    @test [a.xyz for a in unit_cell(lattice_sample("BCC", "Fe", 2.866))] ==
          [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]
    nacl = cell("NaCl", "FCC", 5.64, origin("Na"), Atom("Cl", (0.5, 0.0, 0.0), 0.0))
    @test length(unit_cell(nacl)) == 8
    @test all(0 ≤ x < 1 for a in unit_cell(nacl) for x in a.xyz)

    # CsCl is not BCC: the body-centre translation puts Cl on the Cs site
    @test_throws ArgumentError unit_cell(cell("CsCl", "BCC", 4.12, origin("Cs"),
                                              Atom("Cl", (0.5, 0.5, 0.5), 0.0)))
    # A centering copy listed in the basis falls on its own site
    @test_throws ArgumentError unit_cell(cell("Cu", "FCC", 3.594, origin("Cu"),
                                              Atom("Cu", (0.0, 0.5, 0.5), 0.0)))
end

@testset "structure factor equals the absence rules for monatomic cells" begin
    for structure in ("SC", "BCC", "FCC"), max_hkl_sq in (3, 20, 60)
        sample = lattice_sample(structure, "Cu", 3.6)
        @test reflections(StructureFactor(), sample, max_hkl_sq) ==
              reflections(AbsenceRules(), sample, max_hkl_sq)
    end

    for file in (SF_XRAY_TOML, SF_ELECTRON_TOML)
        rules, sf = with_model(file, AbsenceRules()), with_model(file, StructureFactor())
        for (structure, element, a, B) in (("SC", "Po", 3.352, 0.0), ("BCC", "Fe", 2.866, 0.325),
                                           ("BCC", "Li", 3.491, 4.81), ("FCC", "Cu", 3.594, 0.55),
                                           ("FCC", "Au", 4.065, 0.0))
            x, y_rules = simulate(rules, structure, element, a, B)
            _, y_sf = simulate(sf, structure, element, a, B)
            @test y_sf ≈ y_rules rtol=1e-10
        end
    end
end

@testset "structure factor of multi-atom cells" begin
    # Diamond (Si): FCC with a second atom at (¼,¼,¼). Beyond the FCC rule,
    # h+k+l = 4n+2 is absent: 200 and 222 vanish, 400 stays.
    si = cell("Si", "FCC", 5.431, origin("Si"), Atom("Si", (0.25, 0.25, 0.25), 0.0))
    indices, _ = reflections(StructureFactor(), si, 20)
    @test indices == [[1, 1, 1], [2, 2, 0], [3, 1, 1], [3, 3, 1], [4, 0, 0]]
    @test !has(indices, [2, 0, 0]) && !has(indices, [2, 2, 2])

    # NaCl: FCC reflections; F(111) = 4(f_Na − f_Cl), F(200) = 4(f_Na + f_Cl)
    nacl = cell("NaCl", "FCC", 5.64, origin("Na"), Atom("Cl", (0.5, 0.0, 0.0), 0.0))
    @test reflections(StructureFactor(), nacl, 20) == Miller_indices("FCC", 20)
    s = 0.2
    fNa, fCl = atomic_form_factor("Na", s), atomic_form_factor("Cl", s)
    F₀ = 4 * (atomic_form_factor("Na", 0.0) + atomic_form_factor("Cl", 0.0))
    @test sf_weights(nacl, [[1, 1, 1], [2, 0, 0]], s) ≈ [(4(fNa - fCl) / F₀)^2, (4(fNa + fCl) / F₀)^2]

    # KCl: K⁺ and Cl⁻ have nearly equal f, so the all-odd reflections are weak, not absent
    kcl = cell("KCl", "FCC", 6.29, origin("K"), Atom("Cl", (0.5, 0.0, 0.0), 0.0))
    @test has(reflections(StructureFactor(), kcl, 20)[1], [1, 1, 1])
    w = sf_weights(kcl, [[1, 1, 1], [2, 0, 0]], 0.1)
    @test w[1] < 0.01 * w[2]

    # CsCl: SC with Cs at the corner and Cl at the body centre; 100 is
    # present, F(100) = f_Cs − f_Cl, F(110) = f_Cs + f_Cl
    cscl = cell("CsCl", "SC", 4.12, origin("Cs"), Atom("Cl", (0.5, 0.5, 0.5), 0.0))
    @test has(reflections(StructureFactor(), cscl, 3)[1], [1, 0, 0])
    fCs, fCl = atomic_form_factor("Cs", s), atomic_form_factor("Cl", s)
    F₀ = atomic_form_factor("Cs", 0.0) + atomic_form_factor("Cl", 0.0)
    @test sf_weights(cscl, [[1, 0, 0], [1, 1, 0]], s) ≈ [((fCs - fCl) / F₀)^2, ((fCs + fCl) / F₀)^2]

    # Ordered Cu₃Au (L1₂): SC; the superlattice reflections 100 and 110 appear
    cu3au = cell("Cu3Au", "SC", 3.75, origin("Au"), Atom("Cu", (0.0, 0.5, 0.5), 0.0),
                 Atom("Cu", (0.5, 0.0, 0.5), 0.0), Atom("Cu", (0.5, 0.5, 0.0), 0.0))
    @test reflections(StructureFactor(), cu3au, 3)[1] == [[1, 0, 0], [1, 1, 0], [1, 1, 1]]

    # B per atom: damping only the Cl atoms changes F(111) and F(200) differently
    nacl_B = cell("NaCl", "FCC", 5.64, origin("Na"), Atom("Cl", (0.5, 0.0, 0.0), 1.0))
    fCl_B = atomic_form_factor("Cl", s) * exp(-s^2)
    F₀ = 4 * (atomic_form_factor("Na", 0.0) + atomic_form_factor("Cl", 0.0))
    @test sf_weights(nacl_B, [[1, 1, 1], [2, 0, 0]], s) ≈
          [(4(fNa - fCl_B) / F₀)^2, (4(fNa + fCl_B) / F₀)^2]

    # Weights are of order 1 and positive in both modes
    for file in (SF_XRAY_TOML, SF_ELECTRON_TOML)
        x, y = simulate(with_model(file, StructureFactor()), nacl)
        @test all(isfinite, y) && all(≥(0), y)
    end
end
