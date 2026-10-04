# The generic simulation of one sample; the mode-specific steps are methods on
# `XRay` (xray.jl) and `Electron` (electron.jl).

"""
    simulate(cfg::XRDConfig, sample::Sample)
    simulate(cfg::XRDConfig, structure::String, element::String, a::Real, B::Real=0.0)

Compute the powder pattern of one sample in the radiation mode `cfg.mode`.
The second method simulates the monatomic sample `lattice_sample(structure,
element, a, B)`.

The steps are the same for every mode; each step is a method on the mode type:
grid (`grid`), reflections up to the cutoff (`max_hkl_sq`, `Miller_indices`),
their centres and multiplicities (`peak_centres`), an angle-dependent weight of
each (`peak_weights`: the Debye–Waller factor, and for X-rays also the
Lorentz–polarization factor and the squared atomic form factor f²/Z²), their widths
at each centre (`peak_widths`), the sum of pseudo-Voigt peaks of area
multiplicity × weight (`sum_peaks`) on the `background`, and multiplicative noise of standard
deviation `cfg.noise_level`.

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
    length(sample.atoms) == 1 ||
        throw(ArgumentError("sample $(sample.name): the absence rules need a one-atom basis"))
    atom = only(sample.atoms)
    mode = cfg.mode
    a = sample.a
    x = grid(mode, cfg.N)

    indices, multiplicities = Miller_indices(sample.centering, max_hkl_sq(mode, a))
    x₀, m = peak_centres(mode, indices, multiplicities, a)
    widths = [peak_widths(mode, xᵢ, cfg) for xᵢ in x₀]
    w_L, w_G = first.(widths), last.(widths)

    A = m .* peak_weights(mode, x₀, atom.element, atom.B)

    y = background(mode, x) .+ sum_peaks(x, x₀, A, w_L, w_G)

    if cfg.noise_level > 0
        y .*= rand(Normal(1, cfg.noise_level), length(x))
        y = max.(y, 0)
    end

    return display_axis(mode, x), y
end

simulate(cfg::XRDConfig, structure::String, element::String, a::Real, B::Real=0.0) =
    simulate(cfg, lattice_sample(structure, element, a, B))
