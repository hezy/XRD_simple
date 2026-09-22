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
    @test cfg.radiation in ("xray", "electron")
    @test cfg.samples isa Vector{Tuple{String,String,Float64}}
    @test all(s[1] in ("SC", "BCC", "FCC") for s in cfg.samples)
    @test all(s[3] > 0 for s in cfg.samples)

    # The file and the Dict methods agree
    @test read_xrd_config(DATA_TOML).samples == read_xrd_config(TOML.parsefile(DATA_TOML)).samples

    cfg = read_xrd_config(base_config())
    @test cfg.two_theta_min ≈ deg2rad(10.0)
    @test cfg.two_theta_max ≈ deg2rad(120.0)
    @test cfg.D === 500.0                          # integer in TOML, Float64 here
    @test cfg.samples == [("FCC", "Ag", 4.079), ("FCC", "Cu", 3.594), ("SC", "Po", 3.352)]
end

@testset "read_xrd_config defaults" begin
    c = with(base_config(), "instrument", "radiation", nothing)
    cfg = read_xrd_config(c)
    @test cfg.radiation == "xray"
    @test cfg.noise_level == 0.0
    @test cfg.voltage_kV == 200.0
    @test cfg.g_min == 0.0
    @test cfg.G_inst == 0.005
    @test cfg.camera_constant == 50.0
    @test cfg.image_px == 800
    @test cfg.beam_stop_mm == 2.5
    @test cfg.ring_phosphor == true
    @test cfg.ring_gamma == 0.5
    @test cfg.ring_noise == 0.0

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
    cfg = read_xrd_config(c)
    @test cfg.radiation == "electron"
    @test isnan(cfg.lambda) && isnan(cfg.U)

    cfg = read_xrd_config(with(base_config(), "instrument", "g_max", nothing))
    @test isnan(cfg.g_max)
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

    # A negative V is valid when the Caglioti FWHM² stays positive
    @test read_xrd_config(with(b, "peak_width", "V", -5e-5)).V == -5e-5

    # Samples
    c = deepcopy(b); c["lattice"]["HCP"] = Dict{String,Any}("Mg" => 3.21)
    @test_throws ArgumentError read_xrd_config(c)
    c = deepcopy(b); c["lattice"]["SC"]["Po"] = -3.352
    @test_throws ArgumentError read_xrd_config(c)
end
