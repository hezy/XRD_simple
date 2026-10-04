# The generic simulation of one sample; the mode-specific steps are methods on
# `XRay` (xray.jl) and `Electron` (electron.jl).

"""
    simulate(cfg::XRDConfig, sample::Sample)
    simulate(cfg::XRDConfig, structure::String, element::String, a::Real, B::Real=0.0)

Compute the powder pattern of one sample in the radiation mode `cfg.mode`.
The second method simulates the monatomic sample `lattice_sample(structure,
element, a, B)`.

The steps are the same for every mode and model. The mode gives the grid
(`grid`), the cutoff (`max_hkl_sq`), the centres of the reflections
(`peak_centres`), the Lorentz–polarization factor (`angular_factor`), s = sin θ / λ
(`scattering_s`) and the widths at each centre (`peak_widths`). The model
`cfg.model` gives the reflections (`reflections`) and their scattering weights
(`scattering_weights`: the form factors and the Debye–Waller factor). The pattern
is the sum of pseudo-Voigt peaks of area multiplicity × angular factor ×
scattering weight (`sum_peaks`) on the `background`, with multiplicative noise
of standard deviation `cfg.noise_level`.

# Arguments
- `cfg::XRDConfig`: Configuration from `read_xrd_config`
- `sample::Sample`: the unit cell; with `AbsenceRules`, a one-atom basis
- `structure::String`: Crystal structure ("SC", "BCC", or "FCC")
- `element::String`: Chemical symbol of the (single) element, e.g. "Fe"
- `a::Real`: Lattice parameter in Angstroms
- `B::Real`: Debye–Waller parameter in Å² (0: no thermal damping)

# Returns
- `(x, y)`: the x axis in display units (2θ in degrees, or g in 1/Å; see
  `axis_label`) and the intensity at each x
"""
function simulate(cfg::XRDConfig, sample::Sample)::Tuple{Vector{Float64}, Vector{Float64}}
    cfg.model isa AbsenceRules ||
        throw(ArgumentError("reflections = \"structure_factor\" is not implemented yet"))
    mode, model, a = cfg.mode, cfg.model, sample.a
    x = grid(mode, cfg.N)

    indices, multiplicities = reflections(model, sample, max_hkl_sq(mode, a))
    x₀, indices, m = peak_centres(mode, indices, multiplicities, a)
    widths = [peak_widths(mode, xᵢ, cfg) for xᵢ in x₀]
    w_L, w_G = first.(widths), last.(widths)

    A = m .* angular_factor(mode, x₀) .*
        scattering_weights(model, mode, sample, indices, scattering_s(mode, x₀))

    y = background(mode, x) .+ sum_peaks(x, x₀, A, w_L, w_G)

    if cfg.noise_level > 0
        y .*= rand(Normal(1, cfg.noise_level), length(x))
        y = max.(y, 0)
    end

    return display_axis(mode, x), y
end

simulate(cfg::XRDConfig, structure::String, element::String, a::Real, B::Real=0.0) =
    simulate(cfg, lattice_sample(structure, element, a, B))
