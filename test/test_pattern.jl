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

@testset "Lorentz_polarization" begin
    @test Lorentz_polarization(π/4) ≈ 1.0
    # Analytic value at 2θ = 30°: (1 + cos²30°) / (sin²15° · cos15°) / (2√2)
    θ = deg2rad(15.0)
    @test Lorentz_polarization(θ) ≈ (1 + cosd(30)^2) / (sind(15)^2 * cosd(15)) / (2√2)
    # Decreasing from low angles to its minimum near 2θ ≈ 100°
    θs = deg2rad.(5.0:5.0:45.0)
    @test issorted(Lorentz_polarization.(θs), rev=true)
    @test_throws ArgumentError Lorentz_polarization(0.0)
    @test_throws ArgumentError Lorentz_polarization(π/2)

end

@testset "Debye_Waller" begin
    @test Debye_Waller(0.5, 0.0) == 1.0
    @test Debye_Waller(0.0, 1.0) == 1.0
    @test Debye_Waller(0.5, 0.4) ≈ exp(-0.2)
    @test Debye_Waller(0.6, 0.4) < Debye_Waller(0.3, 0.4)
    @test_throws ArgumentError Debye_Waller(0.5, -0.1)
end

@testset "atomic_form_factor" begin
    # f(0) = Z for every tabulated neutral atom
    Z = Dict("H" => 1, "Al" => 13, "Fe" => 26, "Cu" => 29, "Ag" => 47, "W" => 74,
             "Au" => 79, "Po" => 84, "Cf" => 98)
    for (el, z) in Z
        @test atomic_form_factor(el, 0.0) ≈ z atol = 0.05
    end
    @test length(FORM_FACTOR_COEFFICIENTS) == 98
    @test all(abs(sum(c[1:6]) - round(sum(c[1:6]))) < 0.05 for c in values(FORM_FACTOR_COEFFICIENTS))
    # Decreasing with s
    @test issorted(atomic_form_factor.("Cu", 0.0:0.1:2.0), rev=true)
    @test_throws ArgumentError atomic_form_factor("Xx", 0.5)
    @test_throws ArgumentError atomic_form_factor("Fe", -0.1)
    @test_throws ArgumentError atomic_form_factor("Fe", 6.5)
end

@testset "electron_form_factor" begin
    # Mott–Bethe: C (Z − f₀(s)) / s², with Z = f₀(0) of the fit
    for el in ("Al", "Fe", "Au"), s in (0.1, 0.5, 2.0)
        Z = atomic_form_factor(el, 0.0)
        @test electron_form_factor(el, s) ≈ 0.023934 * (Z - atomic_form_factor(el, s)) / s^2
    end
    # Finite and continuous at s = 0, and falling with s
    @test electron_form_factor("Fe", 0.0) ≈ electron_form_factor("Fe", 1e-4) rtol = 1e-6
    @test issorted(electron_form_factor.("Fe", 0:0.1:2), rev=true)
    # Falls faster than f₀: f_e(s)/f_e(0) < f₀(s)/Z
    @test electron_form_factor("Fe", 0.5) / electron_form_factor("Fe", 0.0) <
          atomic_form_factor("Fe", 0.5) / atomic_form_factor("Fe", 0.0)
    @test_throws ArgumentError electron_form_factor("Xx", 0.1)
    @test_throws ArgumentError electron_form_factor("Fe", -0.1)
    @test_throws ArgumentError electron_form_factor("Fe", 6.5)
end

@testset "peak_weights" begin
    # X-ray: LP × f² × DW with s = sin θ / λ
    cfg = read_xrd_config(XRAY_TOML)
    λ = cfg.mode.lambda
    two_θ₀ = deg2rad.([20.0, 90.0])
    s = sin.(two_θ₀ ./ 2) ./ λ
    f² = (atomic_form_factor.("Fe", s) ./ atomic_form_factor("Fe", 0.0)) .^ 2
    @test peak_weights(cfg.mode, two_θ₀, "Fe", 0.0) ≈ Lorentz_polarization.(two_θ₀ ./ 2) .* f²
    @test peak_weights(cfg.mode, two_θ₀, "Fe", 0.5) ≈
          Lorentz_polarization.(two_θ₀ ./ 2) .* f² .* exp.(-2 * 0.5 .* s .^ 2)

    # Electron: (f_e/f_e(0))² × DW with s = g/2
    cfg = read_xrd_config(ELECTRON_TOML)
    s = [0.3, 0.6] ./ 2
    fe² = (electron_form_factor.("Fe", s) ./ electron_form_factor("Fe", 0.0)) .^ 2
    @test peak_weights(cfg.mode, [0.3, 0.6], "Fe", 0.0) ≈ fe²
    @test peak_weights(cfg.mode, [0.3, 0.6], "Fe", 0.5) ≈ fe² .* exp.(-2 * 0.5 .* s .^ 2)
    @test 1 > fe²[1] > fe²[2] > 0

    # In a pattern, B lowers a high-angle peak more than a low-angle one
    x, y0 = simulate(read_xrd_config(XRAY_TOML), "SC", "Po", 3.352, 0.0)
    _, y1 = simulate(read_xrd_config(XRAY_TOML), "SC", "Po", 3.352, 1.0)
    peak(y, two_θ) = y[argmin(abs.(x .- two_θ))]
    d(N) = 3.352 / √N
    two_θ(N) = 2asind(λ / (2d(N)))
    @test peak(y1, two_θ(1)) / peak(y0, two_θ(1)) > peak(y1, two_θ(9)) / peak(y0, two_θ(9))
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
        x, y = simulate(cfg, "SC", "Po", a)
        @test length(x) == length(y) == cfg.N
        @test all(y .>= 0)
        @test axis_label(cfg.mode) == label
    end

    # X-ray: x is 2θ in degrees, and the (100) peak lies at 2θ_B
    cfg = read_xrd_config(XRAY_TOML)
    x, y = simulate(cfg, "SC", "Po", a)
    @test x[1] ≈ 10.0 && x[end] ≈ 120.0
    i = argmin(abs.(x .- rad2deg(2asin(cfg.mode.lambda / (2a)))))
    @test y[i] ≈ maximum(y[max(1, i-5):i+5])

    # Electron: the (100) peak lies at g = 1/a
    cfg = read_xrd_config(ELECTRON_TOML)
    x, y = simulate(cfg, "SC", "Po", a)
    i = argmin(abs.(x .- 1 / a))
    @test y[i] ≈ maximum(y[max(1, i-5):i+5])

    # Noise: reproducible with the same seed, different with another
    toml = TOML.parsefile(XRAY_TOML)
    toml["instrument"]["noise_level"] = 0.1
    cfg = read_xrd_config(toml)
    Random.seed!(42); _, y1 = simulate(cfg, "SC", "Po", a)
    Random.seed!(43); _, y2 = simulate(cfg, "SC", "Po", a)
    Random.seed!(42); _, y3 = simulate(cfg, "SC", "Po", a)
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
    g, y = simulate(cfg, "SC", "Po", a)
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

@testset "ring_image distortion" begin
    # Noise-free; the (100) ring is measured along the horizontal and the
    # vertical line through the offset pattern centre (both on the pixel grid)
    a = 3.352
    function ring_radii(settings)
        toml = TOML.parsefile(ELECTRON_TOML)
        toml["instrument"]["image_px"] = 401             # 0.3 mm per pixel
        merge!(toml["instrument"], settings)
        cfg = read_xrd_config(toml)
        mode = cfg.mode
        g, y = simulate(cfg, "SC", "Po", a)
        coords, img = ring_image(g, y, mode)
        i_c = argmin(abs.(coords .- mode.ring_centre_y_mm))
        j_c = argmin(abs.(coords .- mode.ring_centre_x_mm))
        r₀ = mode.camera_constant / a
        # Distance from the centre of the brightest pixel within ±0.15 r₀ of r₀;
        # (110) at √2 r₀ stays outside even when compressed by 10 %
        function measured(line, centre, sign)
            near = [k for k in eachindex(coords) if abs(sign * (coords[k] - coords[centre]) - r₀) < 0.15r₀]
            return abs(coords[near[argmax(line[near])]] - coords[centre])
        end
        h = [measured(img[i_c, :], j_c, s) for s in (1, -1)]   # φ = 0, π
        v = [measured(img[:, j_c], i_c, s) for s in (1, -1)]   # φ = ±π/2
        return r₀, h, v, coords[2] - coords[1]
    end

    # Ellipse along φ₀ = 0: semi-axes r₀(1 + η) horizontal, r₀(1 − η) vertical,
    # about the offset centre; φ₀ = 90° exchanges them
    for (φ₀, fh, fv) in ((0, 1.1, 0.9), (90, 0.9, 1.1))
        r₀, h, v, px = ring_radii(Dict("ring_ellipticity" => 0.1, "ring_axis_deg" => φ₀,
                                       "ring_centre_x_mm" => 3.0, "ring_centre_y_mm" => -2.1))
        @test all(abs.(h .- fh * r₀) .≤ px)
        @test all(abs.(v .- fv * r₀) .≤ px)
    end

    # Pincushion: ρ = r₀ (1 + κ (ρ/R)²), solved by iteration
    κ = 0.05
    r₀, h, v, px = ring_radii(Dict("ring_radial_distortion" => κ))
    R = 50.0 * 1.2
    ρ = r₀
    for _ in 1:50
        ρ = r₀ * (1 + κ * (ρ / R)^2)
    end
    @test all(abs.(vcat(h, v) .- ρ) .≤ px)

    # Without distortion the result does not change
    toml = TOML.parsefile(ELECTRON_TOML)
    cfg = read_xrd_config(toml)
    g, y = simulate(cfg, "SC", "Po", a)
    toml["instrument"]["ring_axis_deg"] = 30.0          # an axis without ellipticity
    @test ring_image(g, y, cfg.mode) == ring_image(g, y, read_xrd_config(toml).mode)
end
