# Configuration: radiation modes, the `XRDConfig` struct, and `read_xrd_config`.

"""
    Radiation

Radiation mode of a simulation: `XRay` or `Electron`. Each subtype holds the
instrument parameters of its mode, and the mode-specific steps of `simulate`
(`grid`, `max_hkl_sq`, `peak_centres`, `peak_weights`, `peak_widths`, `background`,
`display_axis`, `axis_label`) are methods on it. The plot title of each mode
(`plot_title`) is in `plotting.jl`.
"""
abstract type Radiation end


"""
    XRay <: Radiation

X-ray instrument: Bragg geometry, pattern over 2θ.

# Fields
- `lambda`: wavelength (Å)
- `two_theta_min`, `two_theta_max`: 2θ range (radians)
- `U`, `V`, `W`: Caglioti parameters (FWHM² in radians² of 2θ)
"""
struct XRay <: Radiation
    lambda::Float64
    two_theta_min::Float64
    two_theta_max::Float64
    U::Float64
    V::Float64
    W::Float64
end


"""
    Electron <: Radiation

Electron instrument: reciprocal-space geometry, pattern over g = 1/d.

# Fields
- `voltage_kV`: accelerating voltage (kV)
- `g_min`, `g_max`: scattering-vector range (1/Å)
- `G_inst`: instrumental Gaussian FWHM (1/Å)
- `camera_constant`, `image_px`, `beam_stop_mm`, `ring_phosphor`, `ring_gamma`,
  `ring_noise`: ring-image settings
- `ring_ellipticity`, `ring_axis` (radians), `ring_centre_x_mm`,
  `ring_centre_y_mm`, `ring_radial_distortion`: geometric distortion of the
  ring image (see `ring_image`)
"""
struct Electron <: Radiation
    voltage_kV::Float64
    g_min::Float64
    g_max::Float64
    G_inst::Float64
    camera_constant::Float64
    image_px::Int
    beam_stop_mm::Float64
    ring_phosphor::Bool
    ring_gamma::Float64
    ring_noise::Float64
    ring_ellipticity::Float64
    ring_axis::Float64
    ring_centre_x_mm::Float64
    ring_centre_y_mm::Float64
    ring_radial_distortion::Float64
end


"""
    ReflectionModel

Method that decides which reflections a sample has and their scattering weights:
`AbsenceRules` or `StructureFactor`, chosen once per run by `reflections` in
`[model]`.
"""
abstract type ReflectionModel end

"""
    AbsenceRules <: ReflectionModel

The fixed systematic-absence rules of the centering (SC, BCC, FCC) for a
monatomic cell. Cannot describe a `[cell.*]` sample.
"""
struct AbsenceRules <: ReflectionModel end

"""
    StructureFactor <: ReflectionModel

The full structure factor F(hkl) = Σⱼ fⱼ exp(2πi (h xⱼ + k yⱼ + l zⱼ)) over the
atoms of the unit cell. The absences follow from F = 0.
"""
struct StructureFactor <: ReflectionModel end


"""
    Atom

One atom of the basis of a unit cell.

# Fields
- `element::String`: chemical symbol
- `xyz::NTuple{3,Float64}`: fractional position in the cubic cell
- `B::Float64`: Debye–Waller parameter (Å²)
"""
struct Atom
    element::String
    xyz::NTuple{3,Float64}
    B::Float64
end


"""
    Sample

One simulated sample: a cubic lattice of centering `centering` (SC, BCC or
FCC) with lattice parameter `a` (Å), and a basis of atoms. The full cell is the
basis repeated by the centering translations. An `element = a` entry of
`[lattice.STRUCTURE]` is a one-atom basis at the origin, named "element-STRUCTURE";
a `[cell.NAME]` section is a sample named NAME.

# Fields
- `name::String`: title of the sample in the plots and the output files
- `centering::String`: "SC", "BCC" or "FCC"
- `a::Float64`: lattice parameter (Å)
- `atoms::Vector{Atom}`: the basis
"""
struct Sample
    name::String
    centering::String
    a::Float64
    atoms::Vector{Atom}
end

Base.:(==)(s::Sample, t::Sample) =
    (s.name, s.centering, s.a, s.atoms) == (t.name, t.centering, t.a, t.atoms)

"""
    lattice_sample(structure, element, a, B=0.0) -> Sample

The monatomic sample of an `element = a` entry of `[lattice.STRUCTURE]`: one
atom at the origin, named "element-structure".
"""
lattice_sample(structure::String, element::String, a::Real, B::Real=0.0) =
    Sample("$element-$structure", structure, a, [Atom(element, (0.0, 0.0, 0.0), B)])


"""
    XRDConfig

Every parameter of one simulation run, read from `data.toml` by
`read_xrd_config`. All defaults are applied and all values are validated there;
no other function reads the TOML file or supplies a default.

# Fields
- `mode::Radiation`: `XRay` or `Electron`, with its instrument parameters
- `model::ReflectionModel`: `AbsenceRules` or `StructureFactor`
- `N::Int`: number of grid points
- `noise_level::Float64`: multiplicative noise standard deviation (0–1)
- `K`: Scherrer constant; `Epsilon`: microstrain; `D`: crystallite size (nm)
- `samples::Vector{Sample}`: sorted by centering, then name; may be empty
"""
struct XRDConfig
    mode::Radiation
    model::ReflectionModel
    N::Int
    noise_level::Float64
    K::Float64
    Epsilon::Float64
    D::Float64
    samples::Vector{Sample}
end


# Read `key` from the TOML table `section` as type T (Float64, Int or Bool).
# A missing key gives `default`; with `default === nothing` it is an error.
function config_value(table::Dict, section::String, key::String, ::Type{T},
                      default=nothing) where {T}
    if !haskey(table, key)
        default === nothing &&
            throw(ArgumentError("[$section] is missing the key \"$key\""))
        return T(default)
    end
    v = table[key]
    if T === Bool
        v isa Bool || throw(ArgumentError("[$section] $key must be true or false, got $(repr(v))"))
    else
        (v isa Real && !(v isa Bool)) ||
            throw(ArgumentError("[$section] $key must be a number, got $(repr(v))"))
        T === Int && !isinteger(v) &&
            throw(ArgumentError("[$section] $key must be an integer, got $v"))
    end
    return T(v)
end

config_check(ok::Bool, msg::String) = ok || throw(ArgumentError(msg))


"""
    read_xrd_config(filename::String) -> XRDConfig
    read_xrd_config(config::Dict) -> XRDConfig

Read the TOML configuration once, apply every default, validate every value,
and return one `XRDConfig`. The `Dict` method takes an already parsed file.

`radiation` selects the mode: `XRay` or `Electron` is constructed from the keys
of that mode. Keys of the other mode are ignored, not read or checked: both
X-ray and electron parameters can stay in one file. The 2θ limits are
converted from degrees to radians. Each uncommented `element = a` entry under
`[lattice.STRUCTURE]` becomes one sample; zero samples is valid. Each element
must be a chemical symbol with an atomic form factor (H to Cf). The optional
section `[debye_waller]` gives the Debye–Waller parameter B (Å²) of each sample:
an `element = B` entry under `[debye_waller.STRUCTURE]`, or else the key
`default` of `[debye_waller]`. B is keyed by structure because it differs
between phases of one element (BCC and FCC Fe). Entries that are not samples
are ignored.

The optional section `[model]` selects the reflection method by the key
`reflections`: `"rules"` (`AbsenceRules`) or `"structure_factor"`
(`StructureFactor`). Each `[cell.NAME]` section is one more sample, with the
keys `lattice` (SC, BCC or FCC), `a` (Å) and `basis`, a non-empty array of
`{ element = "Na", xyz = [0, 0, 0] }` tables with an optional `B` (Å², else the
`[debye_waller]` default). Cells need `reflections = "structure_factor"`.
Sample names must be unique.

# Defaults
`radiation = "xray"`, `noise_level = 0`, `voltage_kV = 200`, `g_min = 0`,
`G_inst = 0.005`, `camera_constant = 50`, `image_px = 800`,
`beam_stop_mm = 2.5`, `ring_phosphor = true`, `ring_gamma = 0.5`,
`ring_noise = 0`, `ring_ellipticity = 0`, `ring_axis_deg = 0`,
`ring_centre_x_mm = 0`, `ring_centre_y_mm = 0`, `ring_radial_distortion = 0`
(no distortion), `[debye_waller] default = 0` (no thermal damping),
`[model] reflections = "rules"`. All other keys of the selected mode are required.

# Throws
- `ArgumentError`: missing section or key, value of the wrong type, unknown
  `radiation` or structure, or a value out of range (for example a negative
  width, size or lattice parameter)
"""
read_xrd_config(filename::String) = read_xrd_config(TOML.parsefile(filename))

function read_xrd_config(config::Dict)
    for section in ("instrument", "peak_width")
        haskey(config, section) && config[section] isa Dict ||
            throw(ArgumentError("the configuration has no [$section] section"))
    end
    inst, pw = config["instrument"], config["peak_width"]

    radiation = get(inst, "radiation", "xray")
    radiation in ("xray", "electron") ||
        throw(ArgumentError("[instrument] radiation must be \"xray\" or \"electron\", got $(repr(radiation))"))

    N           = config_value(inst, "instrument", "N", Int)
    noise_level = config_value(inst, "instrument", "noise_level", Float64, 0.0)
    K           = config_value(pw, "peak_width", "K", Float64)
    Epsilon     = config_value(pw, "peak_width", "Epsilon", Float64)
    D           = config_value(pw, "peak_width", "D", Float64)

    config_check(N ≥ 2, "[instrument] N must be at least 2, got $N")
    config_check(0 ≤ noise_level ≤ 1, "[instrument] noise_level must be between 0 and 1, got $noise_level")
    config_check(K > 0, "[peak_width] K must be positive, got $K")
    config_check(Epsilon ≥ 0, "[peak_width] Epsilon must not be negative, got $Epsilon")
    config_check(D > 0, "[peak_width] D must be positive, got $D")

    mode = radiation == "xray" ? read_xray(inst, pw) : read_electron(inst, pw)

    dw = get(config, "debye_waller", Dict{String,Any}())
    dw isa Dict || throw(ArgumentError("[debye_waller] must be a table"))
    B_default = config_value(dw, "debye_waller", "default", Float64, 0.0)
    config_check(B_default ≥ 0, "[debye_waller] default must not be negative, got $B_default")
    for (key, value) in dw
        key == "default" || (key in ("SC", "BCC", "FCC") && value isa Dict) ||
            throw(ArgumentError("[debye_waller] $key: expected default or a [debye_waller.SC], [debye_waller.BCC] or [debye_waller.FCC] table of element = B entries"))
    end

    model_table = get(config, "model", Dict{String,Any}())
    model_table isa Dict || throw(ArgumentError("[model] must be a table"))
    reflections = get(model_table, "reflections", "rules")
    reflections in ("rules", "structure_factor") ||
        throw(ArgumentError("[model] reflections must be \"rules\" or \"structure_factor\", got $(repr(reflections))"))
    model = reflections == "rules" ? AbsenceRules() : StructureFactor()

    samples = Sample[]
    for (structure, elements) in get(config, "lattice", Dict{String,Any}())
        section = "lattice.$structure"
        structure in ("SC", "BCC", "FCC") ||
            throw(ArgumentError("[$section]: unknown structure; use SC, BCC or FCC"))
        elements isa Dict || throw(ArgumentError("[$section] must be a table of element = a entries"))
        for element in keys(elements)
            a = config_value(elements, section, element, Float64)
            config_check(a > 0, "[$section] $element: lattice parameter must be positive, got $a")
            config_check(haskey(FORM_FACTOR_COEFFICIENTS, element),
                "[$section] $element: unknown element; use a symbol from H to Cf, e.g. Fe")
            B = config_value(get(dw, structure, Dict{String,Any}()), "debye_waller.$structure",
                             element, Float64, B_default)
            config_check(B ≥ 0, "[debye_waller.$structure] $element: B must not be negative, got $B")
            push!(samples, lattice_sample(structure, element, a, B))
        end
    end

    cells = get(config, "cell", Dict{String,Any}())
    cells isa Dict || throw(ArgumentError("[cell] must be a table of [cell.NAME] sections"))
    for (name, cell) in cells
        model isa StructureFactor ||
            throw(ArgumentError("[cell.$name] needs reflections = \"structure_factor\" in [model]; the absence rules describe only one atom per lattice point"))
        push!(samples, read_cell(name, cell, B_default))
    end

    names = [s.name for s in samples]
    allunique(names) ||
        throw(ArgumentError("two samples have the same name: $(first(n for n in names if count(==(n), names) > 1))"))
    sort!(samples, by = s -> (s.centering, s.name))

    return XRDConfig(mode, model, N, noise_level, K, Epsilon, D, samples)
end


# One [cell.NAME] section: lattice, a, and a non-empty basis of atoms.
function read_cell(name::String, cell, B_default::Float64)::Sample
    section = "cell.$name"
    cell isa Dict || throw(ArgumentError("[$section] must be a table with lattice, a and basis"))
    lattice = get(cell, "lattice", nothing)
    lattice in ("SC", "BCC", "FCC") ||
        throw(ArgumentError("[$section] lattice must be \"SC\", \"BCC\" or \"FCC\", got $(repr(lattice))"))
    a = config_value(cell, section, "a", Float64)
    config_check(a > 0, "[$section] a must be positive, got $a")

    basis = get(cell, "basis", nothing)
    (basis isa AbstractVector && !isempty(basis) && all(b -> b isa Dict, basis)) ||
        throw(ArgumentError("[$section] basis must be a non-empty array of { element = ..., xyz = [x, y, z] } tables"))
    atoms = Atom[]
    for (i, b) in enumerate(basis)
        where_ = "[$section] basis atom $i"
        element = get(b, "element", nothing)
        (element isa String && haskey(FORM_FACTOR_COEFFICIENTS, element)) ||
            throw(ArgumentError("$where_: unknown element $(repr(element)); use a symbol from H to Cf, e.g. Fe"))
        xyz = get(b, "xyz", nothing)
        (xyz isa AbstractVector && length(xyz) == 3 && all(v -> v isa Real && !(v isa Bool), xyz)) ||
            throw(ArgumentError("$where_: xyz must be three numbers (fractional coordinates), got $(repr(xyz))"))
        B = config_value(b, "$section.basis", "B", Float64, B_default)
        config_check(B ≥ 0, "$where_: B must not be negative, got $B")
        push!(atoms, Atom(element, Tuple(Float64.(xyz)), B))
    end
    return Sample(name, lattice, a, atoms)
end


# The X-ray keys of [instrument] and [peak_width]; all are required.
function read_xray(inst::Dict, pw::Dict)::XRay
    two_theta_min = deg2rad(config_value(inst, "instrument", "two_theta_min", Float64))
    two_theta_max = deg2rad(config_value(inst, "instrument", "two_theta_max", Float64))
    lambda        = config_value(inst, "instrument", "lambda", Float64)
    U = config_value(pw, "peak_width", "U", Float64)
    V = config_value(pw, "peak_width", "V", Float64)
    W = config_value(pw, "peak_width", "W", Float64)

    config_check(0 ≤ two_theta_min < two_theta_max ≤ π,
        "[instrument] need 0 ≤ two_theta_min < two_theta_max ≤ 180 (degrees)")
    config_check(lambda > 0, "[instrument] lambda must be positive, got $lambda")
    # FWHM² = U tan²θ + V tanθ + W is positive for every tanθ ≥ 0 exactly when:
    config_check(W > 0 && U ≥ 0 && (V ≥ 0 || V^2 < 4U * W),
        "[peak_width] U, V, W give a negative or zero Caglioti FWHM² at some angle")

    return XRay(lambda, two_theta_min, two_theta_max, U, V, W)
end


# The electron keys of [instrument] and [peak_width]; only g_max is required.
function read_electron(inst::Dict, pw::Dict)::Electron
    voltage_kV      = config_value(inst, "instrument", "voltage_kV", Float64, 200.0)
    g_min           = config_value(inst, "instrument", "g_min", Float64, 0.0)
    g_max           = config_value(inst, "instrument", "g_max", Float64)
    G_inst          = config_value(pw, "peak_width", "G_inst", Float64, 0.005)
    camera_constant = config_value(inst, "instrument", "camera_constant", Float64, 50.0)
    image_px        = config_value(inst, "instrument", "image_px", Int, 800)
    beam_stop_mm    = config_value(inst, "instrument", "beam_stop_mm", Float64, 2.5)
    ring_phosphor   = config_value(inst, "instrument", "ring_phosphor", Bool, true)
    ring_gamma      = config_value(inst, "instrument", "ring_gamma", Float64, 0.5)
    ring_noise      = config_value(inst, "instrument", "ring_noise", Float64, 0.0)
    η               = config_value(inst, "instrument", "ring_ellipticity", Float64, 0.0)
    ring_axis       = deg2rad(config_value(inst, "instrument", "ring_axis_deg", Float64, 0.0))
    x_c             = config_value(inst, "instrument", "ring_centre_x_mm", Float64, 0.0)
    y_c             = config_value(inst, "instrument", "ring_centre_y_mm", Float64, 0.0)
    κ               = config_value(inst, "instrument", "ring_radial_distortion", Float64, 0.0)

    config_check(voltage_kV > 0, "[instrument] voltage_kV must be positive, got $voltage_kV")
    config_check(0 ≤ g_min < g_max, "[instrument] need 0 ≤ g_min < g_max")
    config_check(G_inst > 0, "[peak_width] G_inst must be positive, got $G_inst")
    config_check(camera_constant > 0, "[instrument] camera_constant must be positive, got $camera_constant")
    config_check(image_px ≥ 2, "[instrument] image_px must be at least 2, got $image_px")
    config_check(beam_stop_mm ≥ 0, "[instrument] beam_stop_mm must not be negative, got $beam_stop_mm")
    config_check(ring_gamma > 0, "[instrument] ring_gamma must be positive, got $ring_gamma")
    config_check(0 ≤ ring_noise ≤ 1, "[instrument] ring_noise must be between 0 and 1, got $ring_noise")
    config_check(0 ≤ η < 1, "[instrument] ring_ellipticity must be in [0, 1), got $η")
    # The pattern centre stays inside the frame; with |κ| < 0.1 the radial
    # correction is then monotonic up to the farthest corner.
    r_max = camera_constant * g_max
    config_check(max(abs(x_c), abs(y_c)) ≤ r_max,
        "[instrument] ring_centre_x_mm and ring_centre_y_mm must lie within ±$r_max mm (camera_constant · g_max)")
    config_check(abs(κ) < 0.1, "[instrument] ring_radial_distortion must be between -0.1 and 0.1, got $κ")

    return Electron(voltage_kV, g_min, g_max, G_inst, camera_constant, image_px,
                    beam_stop_mm, ring_phosphor, ring_gamma, ring_noise,
                    η, ring_axis, x_c, y_c, κ)
end
