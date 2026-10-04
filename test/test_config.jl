using Test
using TOML

const DATA_TOML = joinpath(@__DIR__, "..", "data.toml")

# A minimal valid configuration as a parsed-TOML Dict; each test edits a copy.
function base_config(radiation="xray")
    Dict{String,Any}(
        "instrument" => Dict{String,Any}(
            "radiation" => radiation,
            "two_theta_min" => 10.0, "two_theta_max" => 120.0, "lambda" => 1.5418,
            "g_max" => 1.2, "N" => 1000),
        "peak_width" => Dict{String,Any}(
            "U" => 1e-4, "V" => 5e-5, "W" => 1e-5, "K" => 0.9, "Epsilon" => 0.001, "D" => 500),
        "lattice" => Dict{String,Any}(
            "SC" => Dict{String,Any}("Po" => 3.352),
            "FCC" => Dict{String,Any}("Cu" => 3.594, "Ag" => 4.079)))
end

function with(cfg, section, key, value)
    c = deepcopy(cfg)
    value === nothing ? delete!(c[section], key) : (c[section][key] = value)
    return c
end

@testset "read_xrd_config" begin
    cfg = read_xrd_config(DATA_TOML)
    @test cfg isa XRDConfig
    @test cfg.mode isa Radiation
    @test cfg.model isa ReflectionModel
    @test cfg.samples isa Vector{Sample}
    @test all(s.centering in ("SC", "BCC", "FCC") for s in cfg.samples)
    @test all(s.a > 0 for s in cfg.samples)
    @test all(atom.B ≥ 0 for s in cfg.samples for atom in s.atoms)

    # The file and the Dict methods agree
    @test read_xrd_config(DATA_TOML).samples == read_xrd_config(TOML.parsefile(DATA_TOML)).samples

    cfg = read_xrd_config(base_config())
    @test cfg.mode isa XRay
    @test cfg.mode.two_theta_min ≈ deg2rad(10.0)
    @test cfg.mode.two_theta_max ≈ deg2rad(120.0)
    @test cfg.mode.lambda == 1.5418
    @test cfg.D === 500.0                          # integer in TOML, Float64 here
    @test cfg.model isa AbsenceRules
    @test cfg.samples == [lattice_sample("FCC", "Ag", 4.079), lattice_sample("FCC", "Cu", 3.594),
                          lattice_sample("SC", "Po", 3.352)]
    @test cfg.samples[1] == Sample("Ag-FCC", "FCC", 4.079, [Atom("Ag", (0.0, 0.0, 0.0), 0.0)])

    # [debye_waller]: an entry of the sample's structure, else the default;
    # others are ignored (Ag under BCC does not apply to FCC Ag)
    c = deepcopy(base_config())
    c["debye_waller"] = Dict{String,Any}("default" => 0.5,
        "FCC" => Dict{String,Any}("Cu" => 0.55, "Fe" => 0.56),
        "BCC" => Dict{String,Any}("Ag" => 0.9, "Fe" => 0.33))
    @test [only(s.atoms).B for s in read_xrd_config(c).samples] == [0.5, 0.55, 0.5]
end

# A rock-salt cell as parsed TOML: FCC lattice, Na at the origin, Cl at (½,0,0)
nacl_cell() = Dict{String,Any}("lattice" => "FCC", "a" => 5.64,
    "basis" => Any[Dict{String,Any}("element" => "Na", "xyz" => Any[0, 0, 0]),
                   Dict{String,Any}("element" => "Cl", "xyz" => Any[0.5, 0, 0], "B" => 1.2)])

function with_cell(cfg, cell=nacl_cell(); name="NaCl", reflections="structure_factor")
    c = deepcopy(cfg)
    c["model"] = Dict{String,Any}("reflections" => reflections)
    c["cell"] = Dict{String,Any}(name => cell)
    return c
end

@testset "read_xrd_config [model] and [cell]" begin
    @test read_xrd_config(base_config()).model isa AbsenceRules
    c = deepcopy(base_config()); c["model"] = Dict{String,Any}("reflections" => "structure_factor")
    @test read_xrd_config(c).model isa StructureFactor

    c = with_cell(base_config())
    c["debye_waller"] = Dict{String,Any}("default" => 0.4)
    cfg = read_xrd_config(c)
    @test [s.name for s in cfg.samples] == ["Ag-FCC", "Cu-FCC", "NaCl", "Po-SC"]
    @test cfg.samples[3] == Sample("NaCl", "FCC", 5.64,
        [Atom("Na", (0.0, 0.0, 0.0), 0.4), Atom("Cl", (0.5, 0.0, 0.0), 1.2)])
end

@testset "read_xrd_config defaults" begin
    c = with(base_config(), "instrument", "radiation", nothing)
    cfg = read_xrd_config(c)
    @test cfg.mode isa XRay
    @test cfg.noise_level == 0.0

    e = read_xrd_config(base_config("electron")).mode
    @test e isa Electron
    @test e.voltage_kV == 200.0
    @test e.g_min == 0.0
    @test e.g_max == 1.2
    @test e.G_inst == 0.005
    @test e.camera_constant == 50.0
    @test e.image_px == 800
    @test e.beam_stop_mm == 2.5
    @test e.ring_phosphor == true
    @test e.ring_gamma == 0.5
    @test e.ring_noise == 0.0
    @test (e.ring_ellipticity, e.ring_axis, e.ring_centre_x_mm, e.ring_centre_y_mm,
           e.ring_radial_distortion) == (0.0, 0.0, 0.0, 0.0, 0.0)
    @test read_xrd_config(with(base_config("electron"), "instrument", "ring_axis_deg", 90)).mode.ring_axis ≈ π / 2

    # Zero samples is valid, with or without a [lattice] section
    c = deepcopy(base_config()); delete!(c, "lattice")
    @test isempty(read_xrd_config(c).samples)
    c = deepcopy(base_config()); c["lattice"] = Dict{String,Any}("BCC" => Dict{String,Any}())
    @test isempty(read_xrd_config(c).samples)
end

@testset "read_xrd_config ignores the unused mode" begin
    # Electron mode needs no X-ray keys, and the reverse
    c = base_config("electron")
    for key in ("two_theta_min", "two_theta_max", "lambda")
        c = with(c, "instrument", key, nothing)
    end
    for key in ("U", "V", "W")
        c = with(c, "peak_width", key, nothing)
    end
    @test read_xrd_config(c).mode isa Electron

    @test read_xrd_config(with(base_config(), "instrument", "g_max", nothing)).mode isa XRay

    # Keys of the unused mode are not checked, even when invalid
    @test read_xrd_config(with(base_config(), "peak_width", "G_inst", -1.0)).mode isa XRay
    @test read_xrd_config(with(base_config("electron"), "instrument", "lambda", "x")).mode isa Electron
end

@testset "read_xrd_config errors" begin
    b = base_config()
    e = base_config("electron")

    c = deepcopy(b); delete!(c, "peak_width")
    @test_throws ArgumentError read_xrd_config(c)
    @test_throws ArgumentError read_xrd_config(with(b, "instrument", "radiation", "neutron"))
    @test_throws ArgumentError read_xrd_config(with(b, "instrument", "N", nothing))
    @test_throws ArgumentError read_xrd_config(with(b, "instrument", "lambda", nothing))
    @test_throws ArgumentError read_xrd_config(with(e, "instrument", "g_max", nothing))
    @test_throws ArgumentError read_xrd_config(with(b, "peak_width", "D", nothing))

    # Wrong types
    @test_throws ArgumentError read_xrd_config(with(b, "instrument", "lambda", "1.54"))
    @test_throws ArgumentError read_xrd_config(with(b, "instrument", "N", 10.5))
    @test_throws ArgumentError read_xrd_config(with(e, "instrument", "ring_phosphor", 1))

    # Out-of-range values
    @test_throws ArgumentError read_xrd_config(with(b, "instrument", "N", 1))
    @test_throws ArgumentError read_xrd_config(with(b, "instrument", "noise_level", 1.5))
    @test_throws ArgumentError read_xrd_config(with(b, "instrument", "lambda", -1.0))
    @test_throws ArgumentError read_xrd_config(with(b, "instrument", "two_theta_max", 5.0))
    @test_throws ArgumentError read_xrd_config(with(b, "peak_width", "D", -500.0))
    @test_throws ArgumentError read_xrd_config(with(b, "peak_width", "K", 0.0))
    @test_throws ArgumentError read_xrd_config(with(b, "peak_width", "Epsilon", -0.001))
    @test_throws ArgumentError read_xrd_config(with(b, "peak_width", "W", -1e-5))
    @test_throws ArgumentError read_xrd_config(with(b, "peak_width", "V", -1e-3))
    @test_throws ArgumentError read_xrd_config(with(e, "peak_width", "G_inst", -0.005))
    @test_throws ArgumentError read_xrd_config(with(e, "instrument", "g_min", 2.0))
    @test_throws ArgumentError read_xrd_config(with(e, "instrument", "image_px", 1))
    @test_throws ArgumentError read_xrd_config(with(e, "instrument", "ring_ellipticity", -0.01))
    @test_throws ArgumentError read_xrd_config(with(e, "instrument", "ring_ellipticity", 1.0))
    @test_throws ArgumentError read_xrd_config(with(e, "instrument", "ring_centre_x_mm", 61.0))   # r_max = 60 mm
    @test_throws ArgumentError read_xrd_config(with(e, "instrument", "ring_centre_y_mm", -61.0))
    @test_throws ArgumentError read_xrd_config(with(e, "instrument", "ring_radial_distortion", 0.1))

    # A negative V is valid when the Caglioti FWHM² stays positive
    @test read_xrd_config(with(b, "peak_width", "V", -5e-5)).mode.V == -5e-5

    # Samples
    c = deepcopy(b); c["lattice"]["HCP"] = Dict{String,Any}("Mg" => 3.21)
    @test_throws ArgumentError read_xrd_config(c)
    c = deepcopy(b); c["lattice"]["SC"]["Po"] = -3.352
    @test_throws ArgumentError read_xrd_config(c)
    c = deepcopy(b); c["lattice"]["SC"]["Xx"] = 3.0
    @test_throws ArgumentError read_xrd_config(c)

    # Debye–Waller
    c = deepcopy(b); c["debye_waller"] = Dict{String,Any}("FCC" => Dict{String,Any}("Cu" => -0.5))
    @test_throws ArgumentError read_xrd_config(c)
    c = deepcopy(b); c["debye_waller"] = Dict{String,Any}("default" => -0.5)
    @test_throws ArgumentError read_xrd_config(c)
    c = deepcopy(b); c["debye_waller"] = Dict{String,Any}("FCC" => Dict{String,Any}("Cu" => "0.5"))
    @test_throws ArgumentError read_xrd_config(c)
    c = deepcopy(b); c["debye_waller"] = Dict{String,Any}("Cu" => 0.5)      # element without structure
    @test_throws ArgumentError read_xrd_config(c)
    c = deepcopy(b); c["debye_waller"] = Dict{String,Any}("HCP" => Dict{String,Any}("Mg" => 1.8))
    @test_throws ArgumentError read_xrd_config(c)

    # [model]
    c = deepcopy(b); c["model"] = Dict{String,Any}("reflections" => "kinematic")
    @test_throws ArgumentError read_xrd_config(c)

    # [cell.*]: only with the structure factor, and well formed
    @test_throws ArgumentError read_xrd_config(with_cell(b; reflections="rules"))
    for (key, value) in (("lattice", "HCP"), ("lattice", nothing), ("a", -5.64), ("a", nothing),
                         ("basis", Any[]), ("basis", nothing), ("basis", "Na"))
        @test_throws ArgumentError read_xrd_config(with_cell(b, with(Dict("x" => nacl_cell()), "x", key, value)["x"]))
    end
    for (key, value) in (("element", "Xx"), ("element", nothing), ("xyz", Any[0, 0]),
                         ("xyz", Any[0, "0", 0]), ("xyz", nothing), ("B", -1.0))
        cell = nacl_cell()
        value === nothing ? delete!(cell["basis"][2], key) : (cell["basis"][2][key] = value)
        @test_throws ArgumentError read_xrd_config(with_cell(b, cell))
    end
    @test_throws ArgumentError read_xrd_config(with_cell(b; name="Cu-FCC"))     # name taken
    cell = nacl_cell(); cell["lattice"] = "BCC"; cell["basis"][2]["xyz"] = Any[0.5, 0.5, 0.5]
    @test_throws ArgumentError read_xrd_config(with_cell(b, cell))               # Cl on the Na site
end
