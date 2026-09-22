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
    electron_wavelength(V::Float64)::Float64

Relativistic de Broglie wavelength (Å) of an electron accelerated through `V`
volts. E.g. 200 kV → 0.0251 Å. Not needed to place the g-axis peaks (which are
purely geometric); used for labelling and as a hook for future camera-length /
ring-radius extensions.
"""
function electron_wavelength(V::Float64)::Float64
    V > 0 || throw(ArgumentError("Accelerating voltage must be positive"))
    return 12.2643 / sqrt(V * (1 + 0.978476e-6 * V))
end


"""
    ed_max_hkl_sq(a::Float64, g_max::Float64)::Int

Largest h²+k²+l² with g = √(h²+k²+l²)/a ≤ g_max — i.e. reflections that fall
within the plotted detector range. The electron analogue of `bragg_max_hkl_sq`;
at electron wavelengths the Bragg `sinθ ≤ 1` bound is never binding, so the
detector range sets the cutoff instead.
"""
function ed_max_hkl_sq(a::Float64,
                       g_max::Float64
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
- `g::Float64`: Scattering vector (1/Å) at which to evaluate the width
- `K::Float64`: Scherrer constant (≈ 0.9)
- `E::Float64`: Microstrain (dimensionless)
- `D::Float64`: Crystallite size in nanometres (converted to Å internally)
"""
function Lorentzian_peaks_width_g(g::Float64,
                                  K::Float64,
                                  E::Float64,
                                  D::Float64
                                  )::Float64
    D > 0 || throw(ArgumentError("Crystallite size D must be positive"))
    D_Å = D * 10.0                      # nm → Å
    return K / D_Å + 2 * E * g
end


"""
    background(mode::Electron, g::Vector{Float64})::Vector{Float64}

Simplified powder-ED background over g: an exponential central-beam tail plus a
constant inelastic (plasmon) floor. Non-negative. Noise is applied by
`simulate`, not here.
"""
function background(::Electron, g::Vector{Float64})::Vector{Float64}
    return @. CENTRAL_BEAM_AMPLITUDE * exp(-CENTRAL_BEAM_DECAY * g) + INELASTIC_LEVEL
end


# Electron steps of `simulate`. The grid and the peak centres are g (1/Å). The
# Lorentzian width is Scherrer size + strain in g-space; the Gaussian width is
# the constant instrumental point-spread G_inst, which replaces the Caglioti
# U/V/W terms (degenerate at θ ≈ 0).

grid(m::Electron, N::Int) = collect(LinRange(m.g_min, m.g_max, N))

max_hkl_sq(m::Electron, a::Float64) = ed_max_hkl_sq(a, m.g_max)

peak_centres(::Electron, indices::Vector{Vector{Int}},
             multiplicities::Vector{Int}, a::Float64) = g_list(indices, a), multiplicities

peak_widths(m::Electron, g₀::Float64, cfg::XRDConfig) =
    (Lorentzian_peaks_width_g(g₀, cfg.K, cfg.Epsilon, cfg.D), m.G_inst)

display_axis(::Electron, g::Vector{Float64}) = g
axis_label(::Electron) = "g (1/Å)"


"""
    reflection_table(structure, a, g_max)

Discrete answer key for the ring pattern. Returns every allowed reflection family
with g = √(h²+k²+l²)/a ≤ g_max, as a NamedTuple of equal-length vectors sorted by
g:

- `indices`      : canonical `[h,k,l]` representatives
- `N`            : N = h²+k²+l² (the ring's squared-index; ring r² ∝ N)
- `g`            : scattering vector g = √N/a (1/Å) — ring radius is `camera_constant·g`
- `multiplicity` : reflection multiplicity (relative ring brightness, geometric)

This is the hidden key students reconstruct from measured ring radii (r² ratios →
N-sequence → SC/BCC/FCC selection rule → lattice constant a).
"""
function reflection_table(structure::String,
                          a::Float64,
                          g_max::Float64
                          )::NamedTuple
    a > 0 || throw(ArgumentError("Lattice parameter must be positive"))
    g_max > 0 || throw(ArgumentError("g_max must be positive"))

    max_hkl_sq = ed_max_hkl_sq(a, g_max)
    indices, multiplicities = Miller_indices(structure, max_hkl_sq)
    g = g_list(indices, a)
    N = [h^2 + k^2 + l^2 for (h, k, l) in indices]

    perm = sortperm(g)
    return (indices = indices[perm],
            N = N[perm],
            g = g[perm],
            multiplicity = multiplicities[perm])
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
    ring_image(g, y, mode::Electron) -> (coords, img)

Map the 1D powder electron-diffraction profile `y(g)` to a 2D Debye–Scherrer
ring pattern. A powder pattern is rotationally symmetric, so the image is a pure
radial lookup: a pixel at distance r (mm) from the centre maps to g = r /
`camera_constant` (1/Å) and takes intensity `y(g)`. This is the SAED-style ring
image students measure — ring radius r = `camera_constant`·g, so r² ∝ N.

# Arguments
- `g::Vector{Float64}`: profile g-grid (1/Å), uniform, from `simulate`
- `y::Vector{Float64}`: profile intensity at each g
- `mode::Electron`: ring settings — `camera_constant` (λL, mm·Å), `image_px`
  (side length, pixels), `beam_stop_mm` (central radius blanked to the floor),
  `ring_gamma` (display gamma, <1 lifts faint outer rings), `ring_noise`
  (per-pixel multiplicative noise, seeded upstream)

# Returns
- `coords::Vector{Float64}`: pixel positions (mm) along both axes, centred on 0
- `img::Matrix{Float64}`: `image_px × image_px` intensities, normalised to
  [0, 1] and gamma-compressed; `img[i, j]` is at `(coords[i], coords[j])`

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
    @inbounds for j in 1:image_px
        yj = coords[j]
        for i in 1:image_px
            r = hypot(coords[i], yj)         # mm from centre
            img[i, j] = r < beam_stop_mm ? floor_val :
                        radial_profile_value(y, g_min, g_max, r / camera_constant)
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
