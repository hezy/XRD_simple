# The generic simulation of one sample; the mode-specific steps are methods on
# `XRay` (xray.jl) and `Electron` (electron.jl).

"""
    simulate(cfg::XRDConfig, structure::String, a::Real)

Compute the powder pattern of one sample in the radiation mode `cfg.mode`.

The steps are the same for every mode; each step is a method on the mode type:
grid (`grid`), reflections up to the cutoff (`max_hkl_sq`, `Miller_indices`),
their centres and multiplicities (`peak_centres`), their widths at each centre
(`peak_widths`), the sum of multiplicity-weighted pseudo-Voigt peaks
(`sum_peaks`) on the `background`, and multiplicative noise of standard
deviation `cfg.noise_level`.

# Arguments
- `cfg::XRDConfig`: Configuration from `read_xrd_config`
- `structure::String`: Crystal structure ("SC", "BCC", or "FCC")
- `a::Real`: Lattice parameter in Angstroms

# Returns
- `(x, y)`: the x axis in display units (2θ in degrees, or g in 1/Å; see
  `axis_label`) and the intensity at each x
"""
function simulate(cfg::XRDConfig,
                  structure::String,
                  a::Real
                  )::Tuple{Vector{Float64}, Vector{Float64}}
    mode = cfg.mode
    x = grid(mode, cfg.N)

    indices, multiplicities = Miller_indices(structure, max_hkl_sq(mode, a))
    x₀, m = peak_centres(mode, indices, multiplicities, a)
    widths = [peak_widths(mode, xᵢ, cfg) for xᵢ in x₀]
    w_L, w_G = first.(widths), last.(widths)

    y = background(mode, x) .+ sum_peaks(x, x₀, m, w_L, w_G)

    if cfg.noise_level > 0
        y .*= rand(Normal(1, cfg.noise_level), length(x))
        y = max.(y, 0)
    end

    return display_axis(mode, x), y
end
