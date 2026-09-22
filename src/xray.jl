# X-ray diffraction: Bragg geometry over 2θ, Caglioti and Scherrer widths,
# background, and the `XRay` methods of `simulate`.

# Background model parameters (X-ray)
const AIR_SCATTER_AMPLITUDE = 50.0
const AIR_SCATTER_DECAY = 5.0
const FLUORESCENCE_LEVEL = 10.0
const AMORPHOUS_AMPLITUDE = 25.0
const AMORPHOUS_CENTER = 0.18
const AMORPHOUS_WIDTH = 0.08


"""
    Gaussian_peaks_width(θ, U, V, W)

Calculate the Gaussian peak width using the Caglioti formula.

The Caglioti formula models instrumental resolution as a function of
Bragg angle: FWHM² = U·tan²(θ) + V·tan(θ) + W

# Arguments
- `θ::Real`: Bragg angle θ (half of 2θ) in radians
- `U::Real`: Caglioti parameter, rad² (typically positive)
- `V::Real`: Caglioti parameter, rad² (typically negative)
- `W::Real`: Caglioti parameter, rad² (typically positive)

# Returns
- `Float64`: Gaussian FWHM at θ, in radians of 2θ
"""
function Gaussian_peaks_width(θ::Real,
                              U::Real,
                              V::Real,
                              W::Real
                              )::Float64

        return √(U * tan(θ)^2 + V * tan(θ) + W)
end


"""
    Lorentzian_peaks_width(θ, K, E, λ, D)

Calculate the Lorentzian peak width from sample broadening effects.

Combines crystallite size broadening (Scherrer equation) and
microstrain broadening (Stokes-Wilson equation).

# Arguments
- `θ::Real`: Bragg angle θ (half of 2θ) in radians
- `K::Real`: Scherrer constant (typically ≈ 0.9)
- `E::Real`: Microstrain (dimensionless)
- `λ::Real`: X-ray wavelength in Angstroms
- `D::Real`: Crystallite size in nanometers (converted to Å internally)

# Returns
- `Float64`: Lorentzian FWHM at θ, in radians of 2θ
"""
function Lorentzian_peaks_width(θ::Real,
                                K::Real,
                                E::Real,
                                λ::Real,
                                D::Real,
                                )::Float64

    D > 0 || throw(ArgumentError("Crystallite size D must be positive"))
    D_Å = D * 10.0                      # nm → Å

    # Strain broadening (Stokes-Wilson)
    w_L_strain = 4 * E * tan(θ)
    # ε is microstrain

    # Size broadening (Scherrer)
    w_L_size = K * λ / (D_Å * cos(θ))
    # K is the Scherrer constant (typically ≈ 0.9)
    # λ is wavelength
    # D is crystallite size

    # Combined broadening
    return w_L_strain + w_L_size
end


"""
    bragg_angles(wavelength::Real, d_spacings::AbstractVector{<:Real})::Tuple{Vector{Float64}, Vector{Int}}

Calculate the Bragg diffraction angles (θ) for a given X-ray wavelength and set of crystal plane d-spacings.

Uses Bragg's law: nλ = 2d·sin(θ), where n=1, λ is the wavelength, and d is the interplanar spacing.

# Arguments
- `wavelength::Real`: X-ray wavelength in Angstroms (Å)
- `d_spacings::AbstractVector{<:Real}`: Vector of interplanar spacings in Angstroms (Å)

# Returns
- `Tuple{Vector{Float64}, Vector{Int}}`: 
   - First element: Vector of Bragg angles in radians where |sin(θ)| ≤ 1
   - Second element: Vector of indices corresponding to the valid angles in the original d_spacings

# Examples
```julia
λ = 1.54  # Cu Kα radiation
d = [2.814, 2.024, 1.431]  # d-spacings in Å
angles, valid_indices = bragg_angles(λ, d)
```

# Throws
* ArgumentError: If wavelength ≤ 0 or any d-spacing ≤ 0
"""
function bragg_angles(wavelength::Real,
                      d_spacings::AbstractVector{<:Real}
                      )::Tuple{Vector{Float64}, Vector{Int}}
    wavelength <= 0 && throw(ArgumentError("Wavelength must be positive"))
    any(d_spacings .<= 0) && throw(ArgumentError("d-spacings must be positive"))

    sinθ = wavelength ./ (2 * d_spacings)
    valid_idx = findall(x -> abs(x) <= 1, sinθ)
    angles = asin.(sinθ[valid_idx])
    return angles, valid_idx
end


"""
    bragg_max_hkl_sq(a::Real, λ::Real)::Int

Largest h²+k²+l² for a cubic lattice that admits a real Bragg angle.

Derived from Bragg's law: sin(θ) = λ/(2d) ≤ 1, combined with the cubic
d-spacing 1/d² = (h²+k²+l²)/a². Reflections above this bound have no real
Bragg angle and cannot appear in any scan.

Note: peaks whose *centers* lie outside the scan window can still contribute
through their wings, so the cutoff deliberately does not depend on
`two_theta_max`.
"""
function bragg_max_hkl_sq(a::Real,
                          λ::Real
                          )::Int
    a > 0 || throw(ArgumentError("Lattice parameter must be positive"))
    λ > 0 || throw(ArgumentError("Wavelength must be positive"))

    return floor(Int, (2 * a / λ)^2)
end


"""
    background(mode::XRay, two_θ::AbstractVector{<:Real})::Vector{Float64}

Simplified X-ray background over the 2θ grid (radians), for educational
simulation purposes. Includes common physical effects seen in XRD patterns:
- Air scattering (exponential decay at low angles)
- Fluorescence (constant background)
- Amorphous scattering (a broad Gaussian hump)

The model is written in θ = 2θ/2. Noise is applied by `simulate`, not here.

# Returns
- `Vector{Float64}`: Background intensity at each grid point (non-negative)
"""
function background(::XRay, two_θ::AbstractVector{<:Real})::Vector{Float64}
    θ = two_θ ./ 2

    air_scatter = @. AIR_SCATTER_AMPLITUDE * exp(-AIR_SCATTER_DECAY * θ)
    fluorescence = FLUORESCENCE_LEVEL
    amorphous = @. AMORPHOUS_AMPLITUDE * exp(-((θ - AMORPHOUS_CENTER)^2) / (2 * AMORPHOUS_WIDTH^2))
    return air_scatter .+ fluorescence .+ amorphous
end


# X-ray steps of `simulate`. The grid and the peak centres are 2θ in radians;
# the widths are FWHM in radians of 2θ, evaluated at θ_B.

grid(m::XRay, N::Int) = collect(LinRange(m.two_theta_min, m.two_theta_max, N))

max_hkl_sq(m::XRay, a::Real) = bragg_max_hkl_sq(a, m.lambda)

# Reflections without a real Bragg angle are dropped with their multiplicities.
function peak_centres(m::XRay, indices::AbstractVector{<:AbstractVector{<:Integer}},
                      multiplicities::AbstractVector{<:Integer}, a::Real)
    θ_B, visible = bragg_angles(m.lambda, d_list(indices, a))
    return 2 .* θ_B, multiplicities[visible]
end

function peak_widths(m::XRay, two_θ₀::Real, cfg::XRDConfig)
    θ_B = two_θ₀ / 2
    return (Lorentzian_peaks_width(θ_B, cfg.K, cfg.Epsilon, m.lambda, cfg.D),
            Gaussian_peaks_width(θ_B, m.U, m.V, m.W))
end

display_axis(::XRay, two_θ::AbstractVector{<:Real}) = rad2deg.(two_θ)
axis_label(::XRay) = "2θ (deg)"
