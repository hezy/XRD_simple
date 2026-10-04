# Electron diffraction (1D)
#
# Powder electron diffraction expressed over the scattering vector g = 1/d (1/Å).
# At electron wavelengths (~0.025 Å) every Bragg angle is a fraction of a degree,
# so a 2θ axis is useless; reflections map linearly to g = √(h²+k²+l²)/a and
# Bragg's law drops out. The crystallography (Miller indices, multiplicities) and
# peak profiles are shared with the X-ray path; only the geometry, the reflection
# cutoff, the broadening and the background differ.
#
# Kinematical approximation only — valid for thin specimens; real SAED is dynamical.

# Background model parameters (electron) — central-beam tail + inelastic floor,
# expressed over the scattering-vector axis g = 1/d (1/Å)
const CENTRAL_BEAM_AMPLITUDE = 80.0
const CENTRAL_BEAM_DECAY = 8.0
const INELASTIC_LEVEL = 5.0


"""
    electron_wavelength(V::Real)::Float64

Relativistic de Broglie wavelength (Å) of an electron accelerated through `V`
volts. E.g. 200 kV → 0.0251 Å. Not needed to place the g-axis peaks (which are
purely geometric); used for labelling and as a hook for future camera-length /
ring-radius extensions.
"""
function electron_wavelength(V::Real)::Float64
    V > 0 || throw(ArgumentError("Accelerating voltage must be positive"))
    return 12.2643 / sqrt(V * (1 + 0.978476e-6 * V))
end


# Mott–Bethe constant m e² / (8π ε₀ h²) in 1/Å, for s = sin θ / λ
const MOTT_BETHE = 0.023934


"""
    electron_form_factor(element::String, s::Real)::Float64

Electron scattering factor f_e of a neutral atom, in Å, from the X-ray form
factor by the Mott–Bethe relation:

    f_e(s) = C (Z − f₀(s)) / s²,   C = 0.023934 1/Å,   s = sin θ / λ = g / 2

Z is taken as f₀(0) of the Waasmaier–Kirfel fit, so that

    Z − f₀(s) = Σᵢ aᵢ (1 − exp(−bᵢ s²))

and the relation has the finite limit f_e(0) = C Σᵢ aᵢ bᵢ instead of a 0/0 at
s = 0. f_e falls much faster with s than f₀, so low-g reflections dominate.
Non-relativistic: the factor γ at the accelerating voltage scales every f_e
equally and is omitted.

# Arguments
- `element::String`: Element symbol, e.g. "Fe"
- `s::Real`: sin θ / λ in 1/Å, 0 ≤ s ≤ 6

# Returns
- `Float64`: f_e(s) in Å

# Throws
* ArgumentError: If the element is not in `FORM_FACTOR_COEFFICIENTS`, or s is
  outside [0, 6]

# Examples
```julia
electron_form_factor("Fe", 0.0)    # ≈ 7.0 Å
```
"""
function electron_form_factor(element::String, s::Real)::Float64
    haskey(FORM_FACTOR_COEFFICIENTS, element) ||
        throw(ArgumentError("no atomic form factor for element \"$element\""))
    0 ≤ s ≤ 6 || throw(ArgumentError("sin θ/λ must be in [0, 6] 1/Å, got $s"))
    a₁, a₂, a₃, a₄, a₅, _, b₁, b₂, b₃, b₄, b₅ = FORM_FACTOR_COEFFICIENTS[element]
    s² = s^2
    # (1 − exp(−b s²)) / s², with the limit b at s = 0
    term(a, b) = s² == 0 ? a * b : -a * expm1(-b * s²) / s²
    return MOTT_BETHE * (term(a₁, b₁) + term(a₂, b₂) + term(a₃, b₃) +
                         term(a₄, b₄) + term(a₅, b₅))
end


"""
    ed_max_hkl_sq(a::Real, g_max::Real)::Int

Largest h²+k²+l² with g = √(h²+k²+l²)/a ≤ g_max — i.e. reflections that fall
within the plotted detector range. The electron analogue of `bragg_max_hkl_sq`;
at electron wavelengths the Bragg `sinθ ≤ 1` bound is never binding, so the
detector range sets the cutoff instead.
"""
function ed_max_hkl_sq(a::Real,
                       g_max::Real
                       )::Int
    a > 0 || throw(ArgumentError("Lattice parameter must be positive"))
    g_max > 0 || throw(ArgumentError("g_max must be positive"))

    return max(1, floor(Int, (g_max * a)^2))
end


"""
    Lorentzian_peaks_width_g(g, K, E, D)

Lorentzian FWHM in g-space (1/Å). Size broadening is the Scherrer width in
reciprocal units — constant K/D — and strain broadening is Δg/g = 2E (from
Δd/d = E). Replaces the angle-space `Lorentzian_peaks_width` for electrons.

# Arguments
- `g::Real`: Scattering vector (1/Å) at which to evaluate the width
- `K::Real`: Scherrer constant (≈ 0.9)
- `E::Real`: Microstrain (dimensionless)
- `D::Real`: Crystallite size in nanometres (converted to Å internally)
"""
function Lorentzian_peaks_width_g(g::Real,
                                  K::Real,
                                  E::Real,
                                  D::Real
                                  )::Float64
    D > 0 || throw(ArgumentError("Crystallite size D must be positive"))
    D_Å = D * 10.0                      # nm → Å
    return K / D_Å + 2 * E * g
end


"""
    background(mode::Electron, g::AbstractVector{<:Real})::Vector{Float64}

Simplified powder-ED background over g: an exponential central-beam tail plus a
constant inelastic (plasmon) floor. Non-negative. Noise is applied by
`simulate`, not here.
"""
function background(::Electron, g::AbstractVector{<:Real})::Vector{Float64}
    return @. CENTRAL_BEAM_AMPLITUDE * exp(-CENTRAL_BEAM_DECAY * g) + INELASTIC_LEVEL
end


# Electron steps of `simulate`. The grid and the peak centres are g (1/Å). The
# Lorentzian width is Scherrer size + strain in g-space; the Gaussian width is
# the constant instrumental point-spread G_inst, which replaces the Caglioti
# U/V/W terms (degenerate at θ ≈ 0).

grid(m::Electron, N::Int) = collect(LinRange(m.g_min, m.g_max, N))

max_hkl_sq(m::Electron, a::Real) = ed_max_hkl_sq(a, m.g_max)

peak_centres(::Electron, indices::AbstractVector{<:AbstractVector{<:Integer}},
             multiplicities::AbstractVector{<:Integer}, a::Real) =
    g_list(indices, a), indices, multiplicities

# The Lorentz factor is constant in g.
angular_factor(::Electron, g₀::AbstractVector{<:Real}) = ones(length(g₀))

# s = sin θ / λ = g/2
scattering_s(::Electron, g₀::AbstractVector{<:Real}) = g₀ ./ 2

form_factor(::Electron, element::String, s::Real) = electron_form_factor(element, s)

peak_widths(m::Electron, g₀::Real, cfg::XRDConfig) =
    (Lorentzian_peaks_width_g(g₀, cfg.K, cfg.Epsilon, cfg.D), m.G_inst)

display_axis(::Electron, g::AbstractVector{<:Real}) = g
axis_label(::Electron) = "g (1/Å)"


"""
    reflection_table(model::ReflectionModel, mode::Electron, sample::Sample)

Discrete answer key for the ring pattern. Returns every reflection family of
`sample` present under the reflection method `model` (see `reflections`), with
g = √(h²+k²+l²)/a ≤ g_max, as a NamedTuple of equal-length vectors sorted by g:

- `indices`      : canonical `[h,k,l]` representatives
- `N`            : N = h²+k²+l² (the ring's squared-index; ring r² ∝ N)
- `g`            : scattering vector g = √N/a (1/Å) — ring radius is `camera_constant·g`
- `multiplicity` : reflection multiplicity (relative ring brightness, geometric)
- `weight`       : scattering weight at s = g/2 (see `scattering_weights`); with
  the structure factor, |F|²/F(000)² with the Debye–Waller factor

This is the hidden key students reconstruct from measured ring radii (r² ratios →
N-sequence → SC/BCC/FCC selection rule → lattice constant a).
"""
function reflection_table(model::ReflectionModel, mode::Electron, sample::Sample)::NamedTuple
    indices, multiplicities = reflections(model, sample, ed_max_hkl_sq(sample.a, mode.g_max))
    g = g_list(indices, sample.a)
    N = [h^2 + k^2 + l^2 for (h, k, l) in indices]
    weight = scattering_weights(model, mode, sample, indices, scattering_s(mode, g))

    perm = sortperm(g)
    return (indices = indices[perm],
            N = N[perm],
            g = g[perm],
            multiplicity = multiplicities[perm],
            weight = weight[perm])
end


"""
    radial_profile_value(y, g_min, g_max, gg)

Linear interpolation of the uniform g-grid profile `y` (over `[g_min, g_max]`) at
scattering vector `gg`. Clamps to the end samples outside the grid. Internal
helper for `ring_image`.
"""
@inline function radial_profile_value(y::Vector{Float64},
                                      g_min::Float64,
                                      g_max::Float64,
                                      gg::Float64)::Float64
    n = length(y)
    gg ≤ g_min && return @inbounds y[1]
    gg ≥ g_max && return @inbounds y[n]
    t = (gg - g_min) / (g_max - g_min) * (n - 1)   # 0-based fractional index
    i = floor(Int, t)
    f = t - i
    @inbounds return (1 - f) * y[i + 1] + f * y[i + 2]
end


"""
    ring_true_radius(ρ, φ, mode::Electron)

Undistorted ring radius r (mm) of the point at distance `ρ` (mm) from the
pattern centre and azimuth `φ` (radians, from the +x axis). It inverts the
distortion model of `ring_image`,

    ρ = r · (1 + η cos 2(φ − φ₀)) · (1 + κ (ρ/R)²),   R = camera_constant · g_max

which is written in ρ on the right so that the inverse is explicit. Internal
helper for `ring_image`.
"""
@inline function ring_true_radius(ρ::Float64, φ::Float64, mode::Electron)::Float64
    R = mode.camera_constant * mode.g_max
    return ρ / ((1 + mode.ring_ellipticity * cos(2 * (φ - mode.ring_axis))) *
                (1 + mode.ring_radial_distortion * (ρ / R)^2))
end


"""
    ring_image(g, y, mode::Electron) -> (coords, img)

Map the 1D powder electron-diffraction profile `y(g)` to a 2D Debye–Scherrer
ring pattern. A powder pattern is rotationally symmetric, so the image is a pure
radial lookup: a pixel at distance r (mm) from the centre maps to g = r /
`camera_constant` (1/Å) and takes intensity `y(g)`. This is the SAED-style ring
image students measure — ring radius r = `camera_constant`·g, so r² ∝ N.

Optional geometric distortion, as in a real microscope, moves a ring of radius
r to the distance ρ from the pattern centre at azimuth φ (see
`ring_true_radius`):

    ρ = r · (1 + η cos 2(φ − φ₀)) · (1 + κ (ρ/R)²),   R = camera_constant · g_max

- η, φ₀ (`ring_ellipticity`, `ring_axis`): elliptical distortion from
  projector-lens astigmatism; to first order an ellipse with semi-axes r(1 ± η)
  along φ₀ and φ₀ + 90°. Typical η is 0.005–0.02.
- κ (`ring_radial_distortion`): barrel (κ < 0) or pincushion (κ > 0)
  distortion; the radius error grows as r³.
- (`ring_centre_x_mm`, `ring_centre_y_mm`): pattern centre, offset from the
  image centre. The beam stop is centred on it.
All are zero by default, which gives exact circles about the image centre.

# Arguments
- `g::Vector{Float64}`: profile g-grid (1/Å), uniform, from `simulate`
- `y::Vector{Float64}`: profile intensity at each g
- `mode::Electron`: ring settings — `camera_constant` (λL, mm·Å), `image_px`
  (side length, pixels), `beam_stop_mm` (central radius blanked to the floor),
  `ring_gamma` (display gamma, <1 lifts faint outer rings), `ring_noise`
  (per-pixel multiplicative noise, seeded upstream), and the distortion
  settings above

# Returns
- `coords::Vector{Float64}`: pixel positions (mm) along both axes, centred on 0
- `img::Matrix{Float64}`: `image_px × image_px` intensities, normalised to
  [0, 1] and gamma-compressed; `img[i, j]` is at x = `coords[j]`,
  y = `coords[i]`, the convention of `heatmap`

`plot_ring_image` in `plotting.jl` draws the result.
"""
function ring_image(g::Vector{Float64},
                    y::Vector{Float64},
                    mode::Electron
                    )::Tuple{Vector{Float64}, Matrix{Float64}}
    length(g) == length(y) || throw(DimensionMismatch("g and y must have equal length"))
    camera_constant, image_px = mode.camera_constant, mode.image_px
    beam_stop_mm, gamma, noise_level = mode.beam_stop_mm, mode.ring_gamma, mode.ring_noise

    g_min, g_max = first(g), last(g)
    floor_val = last(y)                      # dark background outside the ring field
    r_max = camera_constant * g_max          # mm, half-frame (image edge midpoint)
    coords = collect(LinRange(-r_max, r_max, image_px))   # mm, both axes

    img = Matrix{Float64}(undef, image_px, image_px)
    x_c, y_c = mode.ring_centre_x_mm, mode.ring_centre_y_mm
    @inbounds for j in 1:image_px
        dx = coords[j] - x_c
        for i in 1:image_px
            dy = coords[i] - y_c
            ρ = hypot(dx, dy)                # mm from the pattern centre
            img[i, j] = ρ < beam_stop_mm ? floor_val :
                        radial_profile_value(y, g_min, g_max,
                                             ring_true_radius(ρ, atan(dy, dx), mode) / camera_constant)
        end
    end

    if noise_level > 0
        img .*= rand(Normal(1, noise_level), size(img))
    end

    # Normalise to [0,1] then gamma-compress so faint high-g rings stay visible.
    lo, hi = extrema(img)
    disp = hi > lo ? @.(((img - lo) / (hi - lo))^gamma) : zero(img)

    return coords, disp
end
