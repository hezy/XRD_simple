using Test
using Random

using TOML

# The fixed reference configurations, independent of the user's data.toml
const XRAY_TOML = joinpath(@__DIR__, "reference", "xray.toml")
const ELECTRON_TOML = joinpath(@__DIR__, "reference", "electron.toml")

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

@testset "peak_widths" begin
    cfg = read_xrd_config(XRAY_TOML)
    for two_θ_B in deg2rad.([20.0, 60.0, 100.0])
        w_L, w_G = peak_widths(cfg.mode, two_θ_B, cfg)
        @test w_L > 0 && w_G > 0
        @test w_G == Gaussian_peaks_width(two_θ_B / 2, cfg.mode.U, cfg.mode.V, cfg.mode.W)
    end

    cfg = read_xrd_config(ELECTRON_TOML)
    w_L, w_G = peak_widths(cfg.mode, 0.5, cfg)
    @test w_L == Lorentzian_peaks_width_g(0.5, cfg.K, cfg.Epsilon, cfg.D)
    @test w_G == cfg.mode.G_inst
end

@testset "simulate" begin
    a = 3.352
    for (file, label) in ((XRAY_TOML, "2θ (deg)"), (ELECTRON_TOML, "g (1/Å)"))
        cfg = read_xrd_config(file)
        x, y = simulate(cfg, "SC", a)
        @test length(x) == length(y) == cfg.N
        @test all(y .>= 0)
        @test axis_label(cfg.mode) == label
    end

    # X-ray: x is 2θ in degrees, and the (100) peak lies at 2θ_B
    cfg = read_xrd_config(XRAY_TOML)
    x, y = simulate(cfg, "SC", a)
    @test x[1] ≈ 10.0 && x[end] ≈ 120.0
    i = argmin(abs.(x .- rad2deg(2asin(cfg.mode.lambda / (2a)))))
    @test y[i] ≈ maximum(y[max(1, i-5):i+5])

    # Electron: the (100) peak lies at g = 1/a
    cfg = read_xrd_config(ELECTRON_TOML)
    x, y = simulate(cfg, "SC", a)
    i = argmin(abs.(x .- 1 / a))
    @test y[i] ≈ maximum(y[max(1, i-5):i+5])

    # Noise: reproducible with the same seed, different with another
    toml = TOML.parsefile(XRAY_TOML)
    toml["instrument"]["noise_level"] = 0.1
    cfg = read_xrd_config(toml)
    Random.seed!(42); _, y1 = simulate(cfg, "SC", a)
    Random.seed!(43); _, y2 = simulate(cfg, "SC", a)
    Random.seed!(42); _, y3 = simulate(cfg, "SC", a)
    @test y1 != y2
    @test y1 == y3
end

@testset "ring_image" begin
    # With per-pixel noise, so that only properties robust to noise are tested
    a = 3.352
    toml = TOML.parsefile(ELECTRON_TOML)
    toml["instrument"]["image_px"] = 201    # odd: one pixel at the centre
    toml["instrument"]["ring_noise"] = 0.1
    cfg = read_xrd_config(toml)
    mode = cfg.mode
    Random.seed!(42)
    g, y = simulate(cfg, "SC", a)
    coords, img = ring_image(g, y, mode)

    @test length(coords) == mode.image_px
    @test size(img) == (mode.image_px, mode.image_px)
    @test coords[end] ≈ mode.camera_constant * mode.g_max
    @test extrema(img) == (0.0, 1.0)

    # The (100) ring lies at r = camera_constant · g = camera_constant / a,
    # in all four directions from the centre
    c = (mode.image_px + 1) ÷ 2
    r₀ = mode.camera_constant / a
    for (line, sign) in ((img[:, c], 1), (img[:, c], -1), (img[c, :], 1), (img[c, :], -1))
        i₀ = argmin(abs.(coords .- sign * r₀))
        window = i₀-3:i₀+3
        @test abs(window[argmax(line[window])] - i₀) ≤ 1
    end

    # The beam stop is darker than the (100) ring
    stop = [img[i, j] for i in eachindex(coords), j in eachindex(coords)
            if hypot(coords[i], coords[j]) < mode.beam_stop_mm]
    ring = img[argmin(abs.(coords .- r₀)), c]
    @test !isempty(stop) && maximum(stop) < ring

    @test_throws DimensionMismatch ring_image(g, y[1:end-1], mode)
end
