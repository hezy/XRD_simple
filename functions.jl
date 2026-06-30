"""
XRD sim
by Hezy Amiel
2023--2026
Julia > 1.8
"""


using Plots; gr()
using SpecialFunctions
#using Random
using Distributions
#using DataFrames
#using CSV
using TOML



""" 
=========
Constants
=========
"""

# Background model parameters (X-ray)
const AIR_SCATTER_AMPLITUDE = 50.0
const AIR_SCATTER_DECAY = 5.0
const FLUORESCENCE_LEVEL = 10.0
const AMORPHOUS_AMPLITUDE = 25.0
const AMORPHOUS_CENTER = 0.18
const AMORPHOUS_WIDTH = 0.08

# Background model parameters (electron) — central-beam tail + inelastic floor,
# expressed over the scattering-vector axis g = 1/d (1/Å)
const CENTRAL_BEAM_AMPLITUDE = 80.0
const CENTRAL_BEAM_DECAY = 8.0
const INELASTIC_LEVEL = 5.0


""" 
=========
Functions
=========
"""


"""
    abstract_peak(θ, θ₀, A, w_L, w_G; cutoff_sigma=5.0, normalize=false)

Template for peak profile functions (Voigt and Pseudo-Voigt).

Arguments:
- `θ::Vector{Float64}`: Position values where to evaluate the peak
- `θ₀::Float64`: Center position of the peak
- `A::Float64`: Peak amplitude (must be positive)
- `w_L::Union{Float64, Vector{Float64}}`: Lorentzian full width at half maximum (FWHM) (must be positive)
- `w_G::Union{Float64, Vector{Float64}}`: Gaussian full width at half maximum (FWHM) (must be positive)

Keyword Arguments:
- `cutoff_sigma::Float64=5.0`: Number of standard deviations beyond which to set peak to zero
- `normalize::Bool=false`: If true, normalize peak height to 1.0

Returns:
- `Vector{Float64}`: Peak intensity at each θ position
"""


"""
    Voigt_peak(θ, θ₀, A, w_L, w_G; cutoff_sigma=5.0, normalize=false)

Computes Voigt peak profile as the convolution of Gaussian and Lorentzian functions
using the complex error function.

See `abstract_peak` for parameter descriptions.

Notes:
- Uses the scaled complementary error function (erfcx) for numerical stability
- More computationally expensive but more accurate than pseudo-Voigt approximation
- Implements bounds checking to improve performance for large datasets
- The cutoff region is based on both Gaussian and Lorentzian widths
"""
# Scalar w's version
function Voigt_peak(θ::Vector{Float64},
                    θ₀::Float64,
                    A::Float64,
                    w_L::Float64,
                    w_G::Float64;
                    cutoff_sigma::Float64=5.0,
                    normalize::Bool=false
                    )::Vector{Float64}

    # Validate parameters
    A > 0 || throw(ArgumentError("Amplitude A must be positive"))
    w_L > 0 || throw(ArgumentError("Lorentzian width w_L must be positive"))
    w_G > 0 || throw(ArgumentError("Gaussian width w_G must be positive"))
    cutoff_sigma > 0 || throw(ArgumentError("cutoff_sigma must be positive"))

    # Initialize output array
    result = zeros(Float64, length(θ))

    # Calculate width parameters
    γ = w_L / 2                                    # Lorentzian HWHM
    σ = w_G / (2√(2log(2)))                        # Gaussian standard deviation

    # Calculate effective width
    w_eff = peak_fwhm(w_L, w_G)

    # Calculate profile only for points within the cutoff region
    for i in eachindex(θ)
        # Check if point is within cutoff region
        if abs(θ[i] - θ₀) ≤ cutoff_sigma * w_eff
            z = -im * (θ[i] - θ₀ + im * γ) / (√2 * σ)    # Complex argument for erfcx
            result[i] = A * real(erfcx(z)) / (√(2π) * σ)
        end
    end

    # Normalize if needed
    if normalize
        maxval = maximum(result)
        if maxval > 0
            result ./= maxval
        end
    end
    
    return result
end

# Vector w's version
function Voigt_peak(θ::Vector{Float64},
                    θ₀::Float64,
                    A::Float64,
                    w_L::Vector{Float64},
                    w_G::Vector{Float64};
                    cutoff_sigma::Float64=5.0,
                    normalize::Bool=false
                    )::Vector{Float64}

    # Validate parameters
    A > 0 || throw(ArgumentError("Amplitude A must be positive"))
    all(w_L .> 0) || throw(ArgumentError("Lorentzian width w_L must be positive"))
    all(w_G .> 0) || throw(ArgumentError("Gaussian width w_G must be positive"))
    cutoff_sigma > 0 || throw(ArgumentError("cutoff_sigma must be positive"))    
    
    # Initialize output array
    result = zeros(Float64, length(θ))                   

    # Calculate width parameters
    γ = w_L / 2                                    # Lorentzian HWHM    
    σ = w_G / (2√(2log(2)))                        # Gaussian standard deviation

    # Calculate effective width 
    w_eff = peak_fwhm(w_L, w_G)    
    
    # Calculate profile only for points within the cutoff region
    for i in eachindex(θ)
        # Check if point is within cutoff region
        if abs(θ[i] - θ₀) ≤ cutoff_sigma * w_eff[i]
            z = -im * (θ[i] - θ₀ + im * γ[i]) / (√2 * σ[i])    # Complex argument for erfcx
            result[i] = A * real(erfcx(z)) / (√(2π) * σ[i])
        end
    end

    # Normalize if needed
    if normalize
        maxval = maximum(result)
        if maxval > 0
            result ./= maxval
        end
    end
    
    return result
end



"""
    pseudo_Voigt_peak(θ, θ₀, A, w_L, w_G; cutoff_sigma=5.0, normalize=false)

Computes pseudo-Voigt peak profile as a linear combination of Gaussian and Lorentzian functions.
The mixing factor is calculated based on the relative widths of the components.

See `abstract_peak` for parameter descriptions.

Notes:
- Mixing factor is computed using the Humps2 approximation
- Implements bounds checking to improve performance for large datasets
"""
# Scalar w's version
function pseudo_Voigt_peak(θ::Vector{Float64},
                           θ₀::Float64,
                           A::Float64,
                           w_L::Float64,
                           w_G::Float64;
                           cutoff_sigma::Float64=5.0,
                           normalize::Bool=false
                           )::Vector{Float64}

    # Validate parameters
    A > 0 || throw(ArgumentError("Amplitude A must be positive"))
    w_L > 0 || throw(ArgumentError("Lorentzian width w_L must be positive"))
    w_G > 0 || throw(ArgumentError("Gaussian width w_G must be positive"))
    cutoff_sigma > 0 || throw(ArgumentError("cutoff_sigma must be positive"))
    
    # Calculate width parameters
    γ = w_L / 2
    σ = w_G / (2√(2log(2)))
    
    # Calculate effective width 
    w_eff = peak_fwhm(w_L, w_G)    
    
    # Calculate mixing factor using Humps2 approximation
    η = 1.36603 * (w_L/w_eff) - 0.47719 * (w_L/w_eff)^2 + 0.11116 * (w_L/w_eff)^3
    
    # Initialize output array
    result = zeros(Float64, length(θ))
 
        # Calculate profile only for points within the cutoff region
    for i in eachindex(θ)
        if abs(θ[i] - θ₀) ≤ cutoff_sigma * w_eff
            # Lorentzian component
            L = γ / (π * ((θ[i] - θ₀)^2 + γ^2))
            # Gaussian component
            G = exp(-(θ[i] - θ₀)^2 / (2σ^2)) / (σ * √(2π))
            # Combined profile
            result[i] = A * (η * L + (1 - η) * G)
        end
    end
    
    if normalize
        maxval = maximum(result)
        if maxval > 0
            result ./= maxval
        end
    end
    
    return result
end

# Vector w's version
function pseudo_Voigt_peak(θ::Vector{Float64},
                         θ₀::Float64,
                         A::Float64,
                         w_L::Vector{Float64},
                         w_G::Vector{Float64};
                         cutoff_sigma::Float64=5.0,
                         normalize::Bool=false
                         )::Vector{Float64}
    
    # Validate parameters
    A > 0 || throw(ArgumentError("Amplitude A must be positive"))
    all(w_L .> 0) || throw(ArgumentError("Lorentzian width w_L must be positive"))
    all(w_G .> 0) || throw(ArgumentError("Gaussian width w_G must be positive"))
    cutoff_sigma > 0 || throw(ArgumentError("cutoff_sigma must be positive"))
    
    # Calculate width parameters
    γ = w_L / 2
    σ = w_G / (2√(2log(2)))
    
    # Calculate effective width 
    w_eff = peak_fwhm(w_L, w_G)    
    
    # Calculate mixing factor using Humps2 approximation
    η = @. 1.36603 * (w_L/w_eff) - 0.47719 * (w_L/w_eff)^2 + 0.11116 * (w_L/w_eff)^3
    
    # Initialize output array
    result = zeros(Float64, length(θ))
 
    
    # Calculate profile only for points within the cutoff region
    for i in eachindex(θ)
        if abs(θ[i] - θ₀) ≤ cutoff_sigma * w_eff[i]
            # Lorentzian component
            L = γ[i] / (π * ((θ[i] - θ₀)^2 + γ[i]^2))
            # Gaussian component
            G = exp(-(θ[i] - θ₀)^2 / (2σ[i]^2)) / (σ[i] * √(2π))
            # Combined profile
            result[i] = A * (η[i] * L + (1 - η[i]) * G)
        end
    end
    
    if normalize
        maxval = maximum(result)
        if maxval > 0
            result ./= maxval
        end
    end
    
    return result
end


# Utility functions that work with both Voigt and pseudo Voigt


"""
    peak_fwhm(w_L::Float64, w_G::Float64)
    peak_fwhm(w_L::Vector{Float64}, w_G::Vector{Float64})

Calculates the full width at half maximum for either Voigt or pseudo-Voigt profile.
Handles both scalar and vector inputs.
"""
# Scalar version
function peak_fwhm(w_L::Float64,
                   w_G::Float64
                   )::Float64
                   
    return 0.5346 * w_L + √(0.2166 * w_L^2 + w_G^2)
end

# Vector version
function peak_fwhm(w_L::Vector{Float64},
                   w_G::Vector{Float64}
                   )::Vector{Float64}
                   
    length(w_L) == length(w_G) || throw(DimensionMismatch("w_L and w_G must have same length"))
    
    return @. 0.5346 * w_L + √(0.2166 * w_L^2 + w_G^2)
end


"""
    Gaussian_peaks_width(θ, U, V, W)

Calculate Gaussian peak widths using the Caglioti formula.

The Caglioti formula models instrumental resolution as a function of
scattering angle: HWHM² = U·tan²(θ) + V·tan(θ) + W

# Arguments
- `θ::Vector{Float64}`: Scattering angles in radians
- `U::Float64`: Caglioti parameter (typically positive)
- `V::Float64`: Caglioti parameter (typically negative)
- `W::Float64`: Caglioti parameter (typically positive)

# Returns
- `Vector{Float64}`: Gaussian FWHM at each angle
"""
function Gaussian_peaks_width(θ::Vector{Float64},
                              U::Float64,
                              V::Float64,
                              W::Float64
                              )::Vector{Float64}

        return @. √(U * tan(θ)^2 + V * tan(θ) + W)
end


"""
    Lorentzian_peaks_width(θ, K, E, λ, D)

Calculate Lorentzian peak widths from sample broadening effects.

Combines crystallite size broadening (Scherrer equation) and
microstrain broadening (Stokes-Wilson equation).

# Arguments
- `θ::Vector{Float64}`: Scattering angles in radians
- `K::Float64`: Scherrer constant (typically ≈ 0.9)
- `E::Float64`: Microstrain (dimensionless)
- `λ::Float64`: X-ray wavelength in Angstroms
- `D::Float64`: Crystallite size in nanometers

# Returns
- `Vector{Float64}`: Lorentzian FWHM at each angle
"""
function Lorentzian_peaks_width(θ::Vector{Float64},
                                K::Float64,
                                E::Float64,
                                λ::Float64,
                                D::Float64,
                                )::Vector{Float64}

    # Strain broadening (Stokes-Wilson)
    w_L_strain = @. 4 * E * tan(θ)
    # ε is microstrain

    # Size broadening (Scherrer)
    w_L_size = @. K * λ / (D * cos(θ))
    # K is the Scherrer constant (typically ≈ 0.9)
    # λ is wavelength
    # D is crystallite size

    # Combined broadening
    return @. w_L_strain + w_L_size        
end



"""
   bragg_angles(wavelength::Float64, d_spacings::Vector{Float64})::Tuple{Vector{Float64}, Vector{Int}}

Calculate the Bragg diffraction angles (θ) for a given X-ray wavelength and set of crystal plane d-spacings.

Uses Bragg's law: nλ = 2d·sin(θ), where n=1, λ is the wavelength, and d is the interplanar spacing.

# Arguments
- `wavelength::Float64`: X-ray wavelength in Angstroms (Å)
- `d_spacings::Vector{Float64}`: Vector of interplanar spacings in Angstroms (Å)

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
function bragg_angles(wavelength::Float64,
                      d_spacings::Vector{Float64}
                      )::Tuple{Vector{Float64}, Vector{Int}}
    wavelength <= 0 && throw(ArgumentError("Wavelength must be positive"))
    any(d_spacings .<= 0) && throw(ArgumentError("d-spacings must be positive"))

    sinθ = wavelength ./ (2 * d_spacings)
    valid_idx = findall(x -> abs(x) <= 1, sinθ)
    angles = asin.(sinθ[valid_idx])
    return angles, valid_idx
end



"""
    d_list(indices::Vector{Vector{Int}}, a::Float64)::Vector{Float64}

Calculate the interplanar distances (d-spacing) for a cubic crystal structure given Miller indices
and lattice parameter.

# Arguments
- `indices::Vector{Vector{Int}}`: Array of Miller indices, where each index is a vector of three 
   integers [h,k,l] representing crystallographic planes
- `a::Float64`: Lattice parameter (unit cell edge length) in appropriate units

# Returns
- `Vector{Float64}`: Array of interplanar distances corresponding to each set of Miller indices

# Throws
- `DimensionMismatch`: If any Miller index vector doesn't contain exactly 3 components
- `DomainError`: If lattice parameter is not positive
"""
function d_list(indices::Vector{Vector{Int}}, a::Float64)::Vector{Float64}
    # Validate lattice parameter
    a > 0 || throw(DomainError(a, "Lattice parameter must be positive"))
    
    # Validate indices structure and dimensions
    for (idx, hkl) in enumerate(indices)
        length(hkl) == 3 || throw(DimensionMismatch(
            "Miller index at position $idx must have exactly 3 components"))
    end
    
    # Pre-allocate output array for better performance
    result = Vector{Float64}(undef, length(indices))
    
    # Calculate d-spacings using direct iteration instead of array comprehension
    # This avoids creating temporary arrays and is more memory efficient
    @inbounds for (i, (h, k, l)) in enumerate(indices)
        result[i] = a / sqrt(h^2 + k^2 + l^2)
    end
    
    return result
end


"""
    sum_peaks(θ, θ_list, multiplicities, w_L, w_G)

Sum pseudo-Voigt peak profiles at given Bragg angles, weighted by multiplicity.

Each entry in `θ_list` is one canonical reflection; its amplitude is the
multiplicity of that family. Summing one weighted peak per family is
mathematically identical to summing every sign+permutation variant at unit
amplitude, and far cheaper.

# Arguments
- `θ::Vector{Float64}`: Scattering angles in radians
- `θ_list::Vector{Float64}`: Peak center positions (Bragg angles) in radians
- `multiplicities::Vector{Int}`: Multiplicity of each reflection family
- `w_L::Vector{Float64}`: Lorentzian FWHM parameters
- `w_G::Vector{Float64}`: Gaussian FWHM parameters

# Returns
- `Vector{Float64}`: Combined peak intensities at each θ
"""
function sum_peaks(θ::Vector{Float64},
                   θ_list::Vector{Float64},
                   multiplicities::Vector{Int},
                   w_L::Vector{Float64},
                   w_G::Vector{Float64},
                   )::Vector{Float64}

    length(θ_list) == length(multiplicities) || throw(DimensionMismatch(
        "θ_list and multiplicities must have same length"))

    y = zeros(length(θ))
    for (θ₀, m) in zip(θ_list, multiplicities)
        y .+= pseudo_Voigt_peak(θ, θ₀, Float64(m), w_L, w_G)
    end
    return y
end


"""
   intensity_vs_angle(θ, indices, multiplicities, λ, a, w_L, w_G)

Calculate X-ray diffraction pattern by summing peak profiles at allowed Bragg angles.

# Arguments
- `θ::Vector{Float64}`: Scattering angles for intensity calculation (radians)
- `indices::Vector{Vector{Int}}`: Canonical Miller indices
- `multiplicities::Vector{Int}`: Multiplicity of each reflection family
- `λ::Float64`: X-ray wavelength (Å)
- `a::Float64`: Lattice parameter (Å)
- `w_L::Vector{Float64}`: Lorentzian width parameters
- `w_G::Vector{Float64}`: Gaussian width parameters

# Returns
- `Vector{Float64}`: XRD intensities at each θ angle

# Throws
- `ArgumentError`: If λ ≤ 0, a ≤ 0, any width ≤ 0, or w_L and w_G have different lengths
- `DimensionMismatch`: If `indices` and `multiplicities` have different lengths
"""
function intensity_vs_angle(θ::Vector{Float64},
                         indices::Vector{Vector{Int}},
                         multiplicities::Vector{Int},
                         λ::Float64,
                         a::Float64,
                         w_L::Vector{Float64},
                         w_G::Vector{Float64}
                         )::Vector{Float64}

   λ <= 0 && throw(ArgumentError("Wavelength must be positive"))
   a <= 0 && throw(ArgumentError("Lattice parameter must be positive"))
   length(indices) == length(multiplicities) || throw(DimensionMismatch(
       "indices and multiplicities must have same length"))
   length(w_L) != length(w_G) && throw(ArgumentError("Width parameter vectors must have same length"))
   any(w_L .<= 0) && throw(ArgumentError("Lorentzian widths must be positive"))
   any(w_G .<= 0) && throw(ArgumentError("Gaussian widths must be positive"))

   θ_list, valid_idx = bragg_angles(λ, d_list(indices, a))
   y = sum_peaks(θ, θ_list, multiplicities[valid_idx], w_L, w_G)
   return y
end


"""
    cubic_multiplicity(h::Int, k::Int, l::Int)::Int

Compute reflection multiplicity for a cubic crystal from a canonical Miller index.

Assumes the input is in canonical form `h ≥ k ≥ l ≥ 0` and not `[0,0,0]`.
Returns the number of (sign, permutation) variants that share the same |G|² =
h²+k²+l², which equals the multiplicity of the {hkl} family under cubic (m-3m)
point-group symmetry.

# Examples
- {100} → 6, {110} → 12, {111} → 8, {210} → 24, {211} → 24, {321} → 48
"""
function cubic_multiplicity(h::Int, k::Int, l::Int)::Int
    nonzero = (h != 0) + (k != 0) + (l != 0)
    sign_variants = 2^nonzero

    # Given h ≥ k ≥ l ≥ 0, repeats appear only as h==k or k==l
    perms = if h == k == l
        1                 # {hhh}
    elseif h == k || k == l
        3                 # {hhl}, {hh0}, {h00}
    else
        6                 # all distinct
    end

    return perms * sign_variants
end


"""
    bragg_max_hkl_sq(a::Float64, λ::Float64)::Int

Largest h²+k²+l² for a cubic lattice that admits a real Bragg angle.

Derived from Bragg's law: sin(θ) = λ/(2d) ≤ 1, combined with the cubic
d-spacing 1/d² = (h²+k²+l²)/a². Reflections above this bound have no real
Bragg angle and cannot appear in any scan.

Note: peaks whose *centers* lie outside the scan window can still contribute
through their wings, so the cutoff deliberately does not depend on
`two_theta_max`.
"""
function bragg_max_hkl_sq(a::Float64,
                          λ::Float64
                          )::Int
    a > 0 || throw(ArgumentError("Lattice parameter must be positive"))
    λ > 0 || throw(ArgumentError("Wavelength must be positive"))

    return floor(Int, (2 * a / λ)^2)
end


"""
    Miller_indices(cell_type::String, max_hkl_sq::Int)::Tuple{Vector{Vector{Int}}, Vector{Int}}

Generate canonical Miller indices and reflection multiplicities for cubic crystals.

Enumerates representatives `h ≥ k ≥ l ≥ 0` (excluding `[0,0,0]`) with
`h² + k² + l² ≤ max_hkl_sq`, applies the systematic absence rule for the given
centering, and returns each allowed reflection together with its multiplicity.
Callers sum one peak per representative, weighted by multiplicity — equivalent
to summing over every sign and permutation, at a fraction of the cost.

# Arguments
- `cell_type::String`: "SC", "BCC", or "FCC"
- `max_hkl_sq::Int`: Upper bound on h²+k²+l² (see `bragg_max_hkl_sq`)

# Returns
- `Vector{Vector{Int}}`: Canonical [h,k,l] representatives
- `Vector{Int}`: Multiplicity of each reflection family

# Throws
- `ArgumentError`: If `cell_type` is not "SC", "BCC", or "FCC"
- `ArgumentError`: If `max_hkl_sq < 1`
"""
function Miller_indices(cell_type::String,
                        max_hkl_sq::Int
                        )::Tuple{Vector{Vector{Int}}, Vector{Int}}

    cell_type in ("SC", "BCC", "FCC") || throw(ArgumentError("cell_type must be 'SC', 'BCC', or 'FCC', got '$cell_type'"))
    max_hkl_sq ≥ 1 || throw(ArgumentError("max_hkl_sq must be ≥ 1, got $max_hkl_sq"))

    max_idx = floor(Int, sqrt(max_hkl_sq))
    indices = Vector{Vector{Int}}()
    multiplicities = Vector{Int}()

    for h in 0:max_idx, k in 0:h, l in 0:k
        (h == 0 && k == 0 && l == 0) && continue
        h^2 + k^2 + l^2 > max_hkl_sq && continue

        allowed = if cell_type == "SC"
            true
        elseif cell_type == "BCC"
            iseven(h + k + l)
        else  # FCC
            (iseven(h) && iseven(k) && iseven(l)) || (isodd(h) && isodd(k) && isodd(l))
        end
        allowed || continue

        push!(indices, [h, k, l])
        push!(multiplicities, cubic_multiplicity(h, k, l))
    end

    return indices, multiplicities
end


"""
    background(θ::Vector{Float64}; noise_level::Float64=0.0)::Vector{Float64}

Generate a simplified XRD background for educational simulation purposes.
Includes common physical effects seen in XRD patterns:
- Air scattering (exponential decay at low angles)
- Fluorescence (constant background)
- Optional random noise

# Arguments
- `θ::Vector{Float64}`: Scattering angles in radians
- `noise_level::Float64=0.0`: Amount of random noise to add (0.0 to 1.0)

# Returns
- `Vector{Float64}`: Background intensity at each angle
"""
function background(θ::Vector{Float64}; 
                   noise_level::Float64=0.0)::Vector{Float64}
    
    # Validate inputs
    0 ≤ noise_level ≤ 1 || throw(DomainError(noise_level, "noise_level must be between 0 and 1"))
    
    # Basic background components
    air_scatter = @. AIR_SCATTER_AMPLITUDE * exp(-AIR_SCATTER_DECAY * θ)
    fluorescence = FLUORESCENCE_LEVEL
    amorphous = @. AMORPHOUS_AMPLITUDE * exp(-((θ - AMORPHOUS_CENTER)^2) / (2 * AMORPHOUS_WIDTH^2))
    base = air_scatter .+ fluorescence .+ amorphous
    
    # Add optional noise
    if noise_level > 0
        noise = noise_level * randn(length(θ))
        return max.(base .+ noise, 0)  # Ensure non-negative intensity
    else
        return base
    end
end


"""
        read_xrd_config(filename::String) ->
        (Dict{String,Any}, Dict{String,Float64}, Vector{Tuple{String,String,Float64}})

Read XRD configuration from TOML file, returning instrument, peak width, and a
flat list of samples to simulate.

# Arguments
- `filename`: Path to TOML configuration file

# Returns
Tuple with:
- "instrument": Dict of instrument parameters (two_theta_min, two_theta_max, N, lambda)
- "peak_width": Dict of peak width parameters (U, V, W, K, Epsilon, D)
- "samples": Vector of `(structure, element, a)` triples, one per uncommented
  lattice entry. Empty blocks contribute nothing; zero samples is a valid result.

Note: Angular parameters (two_theta_min, two_theta_max) are automatically converted to radians.
"""
function read_xrd_config(filename::String)
    config = TOML.parsefile(filename)

    instrument = Dict{String,Any}(
        k => (k in ["two_theta_min", "two_theta_max"] ? deg2rad(v) : v)
        for (k,v) in config["instrument"]
    )

    peak_width = Dict{String,Float64}(config["peak_width"])

    samples = Tuple{String,String,Float64}[]
    for (structure, elements) in config["lattice"]
        for (element, a) in elements
            push!(samples, (structure, element, Float64(a)))
        end
    end
    sort!(samples)

    return instrument, peak_width, samples
end



"""
    do_it_zero(file_name)

Read XRD configuration and return the scattering angle array.

A helper function that extracts instrument parameters from the TOML
config file and constructs the 2θ angle grid. Used to initialize
a DataFrame in the main workflow.

# Arguments
- `file_name::String`: Path to TOML configuration file

# Returns
- `Vector{Float64}`: Scattering angles in radians (half of 2θ range)
"""
function do_it_zero(file_name::String
                    )::Vector{Float64}

    instrument, _, _ = read_xrd_config(file_name)

    if get(instrument, "radiation", "xray") == "electron"
        return collect(LinRange(get(instrument, "g_min", 0.0),
                                instrument["g_max"],
                                instrument["N"]))
    end

    θ = collect(LinRange((instrument["two_theta_min"]/2),
                         (instrument["two_theta_max"]/2),
                         instrument["N"]))
    return θ
end


"""
    compute_peak_widths(θ, peak_width, instrument)

Calculate Lorentzian and Gaussian peak widths from configuration parameters.

Combines Scherrer size + Stokes-Wilson strain broadening for the Lorentzian
component, and Caglioti instrumental resolution for the Gaussian component.

# Arguments
- `θ::Vector{Float64}`: Scattering angles in radians
- `peak_width::Dict{String,Float64}`: Peak width parameters (U, V, W, K, Epsilon, D)
- `instrument::Dict{String,Any}`: Instrument parameters (lambda)

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: (w_L, w_G) peak widths at each angle
"""
function compute_peak_widths(θ::Vector{Float64},
                             peak_width::Dict{String,Float64},
                             instrument::Dict{String,Any}
                             )::Tuple{Vector{Float64}, Vector{Float64}}
    K, ϵ, D = peak_width["K"], peak_width["Epsilon"], peak_width["D"]
    U, V, W = peak_width["U"], peak_width["V"], peak_width["W"]
    λ = instrument["lambda"]

    w_L = Lorentzian_peaks_width(θ, K, ϵ, λ, D)
    w_G = Gaussian_peaks_width(θ, U, V, W)
    return w_L, w_G
end


"""
    compute_xrd_pattern(θ, indices, multiplicities, λ, a, w_L, w_G; noise_level=0.0)

Compute XRD intensity pattern from peak parameters.

Sums pseudo-Voigt peak profiles at Bragg angles, adds background,
and optionally applies multiplicative noise.

# Arguments
- `θ::Vector{Float64}`: Scattering angles in radians
- `indices::Vector{Vector{Int}}`: Canonical Miller indices
- `multiplicities::Vector{Int}`: Multiplicity of each reflection family
- `λ::Float64`: X-ray wavelength in Angstroms
- `a::Float64`: Lattice parameter in Angstroms
- `w_L::Vector{Float64}`: Lorentzian FWHM at each angle
- `w_G::Vector{Float64}`: Gaussian FWHM at each angle

# Keyword Arguments
- `noise_level::Float64=0.0`: Multiplicative noise standard deviation

# Returns
- `Vector{Float64}`: XRD intensities at each angle
"""
function compute_xrd_pattern(θ::Vector{Float64},
                             indices::Vector{Vector{Int}},
                             multiplicities::Vector{Int},
                             λ::Float64,
                             a::Float64,
                             w_L::Vector{Float64},
                             w_G::Vector{Float64};
                             noise_level::Float64=0.0
                             )::Vector{Float64}
    y = background(θ) .+ intensity_vs_angle(θ, indices, multiplicities, λ, a, w_L, w_G)

    if noise_level > 0
        y .*= rand(Normal(1, noise_level), length(θ))
        y = max.(y, 0)
    end

    return y
end


"""
=========================
Electron diffraction (1D)
=========================

Powder electron diffraction expressed over the scattering vector g = 1/d (1/Å).
At electron wavelengths (~0.025 Å) every Bragg angle is a fraction of a degree,
so a 2θ axis is useless; reflections map linearly to g = √(h²+k²+l²)/a and
Bragg's law drops out. The crystallography (Miller indices, multiplicities) and
peak profiles are shared with the X-ray path; only the geometry, the reflection
cutoff, the broadening and the background differ.

Kinematical approximation only — valid for thin specimens; real SAED is dynamical.
"""


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
    g_list(indices::Vector{Vector{Int}}, a::Float64)::Vector{Float64}

Scattering-vector magnitudes g = |G| = √(h²+k²+l²)/a (1/Å) for cubic Miller
indices. The reciprocal-space analogue of `d_list` (g = 1/d).
"""
function g_list(indices::Vector{Vector{Int}}, a::Float64)::Vector{Float64}
    a > 0 || throw(DomainError(a, "Lattice parameter must be positive"))

    result = Vector{Float64}(undef, length(indices))
    @inbounds for (i, hkl) in enumerate(indices)
        length(hkl) == 3 || throw(DimensionMismatch(
            "Miller index at position $i must have exactly 3 components"))
        h, k, l = hkl
        result[i] = sqrt(h^2 + k^2 + l^2) / a
    end
    return result
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
- `g::Vector{Float64}`: Scattering vector grid (1/Å)
- `K::Float64`: Scherrer constant (≈ 0.9)
- `E::Float64`: Microstrain (dimensionless)
- `D::Float64`: Crystallite size in nanometres (converted to Å internally)
"""
function Lorentzian_peaks_width_g(g::Vector{Float64},
                                  K::Float64,
                                  E::Float64,
                                  D::Float64
                                  )::Vector{Float64}
    D > 0 || throw(ArgumentError("Crystallite size D must be positive"))
    D_Å = D * 10.0                      # nm → Å
    return @. K / D_Å + 2 * E * g
end


"""
    compute_peak_widths_g(g, peak_width)

Lorentzian and Gaussian FWHM (1/Å) on the g grid for the electron path.
Lorentzian = Scherrer size + strain (`Lorentzian_peaks_width_g`); Gaussian =
a constant instrumental point-spread `G_inst` (1/Å), replacing the Caglioti
U/V/W terms which are degenerate at θ ≈ 0. `G_inst` defaults to 0.005 1/Å.
"""
function compute_peak_widths_g(g::Vector{Float64},
                               peak_width::Dict{String,Float64}
                               )::Tuple{Vector{Float64}, Vector{Float64}}
    K, ϵ, D = peak_width["K"], peak_width["Epsilon"], peak_width["D"]
    g_inst = get(peak_width, "G_inst", 0.005)

    w_L = Lorentzian_peaks_width_g(g, K, ϵ, D)
    w_G = fill(Float64(g_inst), length(g))
    return w_L, w_G
end


"""
    intensity_vs_g(g, indices, multiplicities, a, w_L, w_G)

Electron-diffraction intensity over the g grid: sum one multiplicity-weighted
pseudo-Voigt per reflection family at its g = √(h²+k²+l²)/a centre. The
g-space counterpart of `intensity_vs_angle`; heights are multiplicity-only,
matching the X-ray path's fidelity.
"""
function intensity_vs_g(g::Vector{Float64},
                        indices::Vector{Vector{Int}},
                        multiplicities::Vector{Int},
                        a::Float64,
                        w_L::Vector{Float64},
                        w_G::Vector{Float64}
                        )::Vector{Float64}
    a <= 0 && throw(ArgumentError("Lattice parameter must be positive"))
    length(indices) == length(multiplicities) || throw(DimensionMismatch(
        "indices and multiplicities must have same length"))
    length(w_L) != length(w_G) && throw(ArgumentError("Width parameter vectors must have same length"))

    g_centers = g_list(indices, a)
    return sum_peaks(g, g_centers, multiplicities, w_L, w_G)
end


"""
    background_electron(g; noise_level=0.0)

Simplified powder-ED background over g: an exponential central-beam tail plus a
constant inelastic (plasmon) floor, with optional additive noise. Non-negative.
"""
function background_electron(g::Vector{Float64};
                             noise_level::Float64=0.0)::Vector{Float64}
    0 ≤ noise_level ≤ 1 || throw(DomainError(noise_level, "noise_level must be between 0 and 1"))

    base = @. CENTRAL_BEAM_AMPLITUDE * exp(-CENTRAL_BEAM_DECAY * g) + INELASTIC_LEVEL

    if noise_level > 0
        noise = noise_level * randn(length(g))
        return max.(base .+ noise, 0)
    else
        return base
    end
end


"""
    compute_ed_pattern(g, indices, multiplicities, a, w_L, w_G; noise_level=0.0)

Full electron-diffraction pattern: background_electron + intensity_vs_g, with
optional multiplicative noise. The g-space analogue of `compute_xrd_pattern`.
"""
function compute_ed_pattern(g::Vector{Float64},
                            indices::Vector{Vector{Int}},
                            multiplicities::Vector{Int},
                            a::Float64,
                            w_L::Vector{Float64},
                            w_G::Vector{Float64};
                            noise_level::Float64=0.0
                            )::Vector{Float64}
    y = background_electron(g) .+ intensity_vs_g(g, indices, multiplicities, a, w_L, w_G)

    if noise_level > 0
        y .*= rand(Normal(1, noise_level), length(g))
        y = max.(y, 0)
    end

    return y
end


"""
    do_it_electron(instrument, peak_width, structure, element, a, plot_theme)

Generate a 1D powder electron-diffraction pattern (intensity vs g = 1/d) for one
sample. Mirrors `do_it`'s X-ray path but in g-space. Returns `(g, intensity,
title, plot)` where `g` is in 1/Å and `title` is `"{element}-{structure}"`.
"""
function do_it_electron(instrument::Dict{String,Any},
                        peak_width::Dict{String,Float64},
                        structure::String,
                        element::String,
                        a::Float64,
                        plot_theme::Symbol
                        )::Tuple{Vector{Float64}, Vector{Float64}, String, Plots.Plot}

    g_min = get(instrument, "g_min", 0.0)
    g_max = instrument["g_max"]
    g = collect(LinRange(g_min, g_max, instrument["N"]))

    max_hkl_sq = ed_max_hkl_sq(a, g_max)
    indices, multiplicities = Miller_indices(structure, max_hkl_sq)
    w_L, w_G = compute_peak_widths_g(g, peak_width)

    noise_level = get(instrument, "noise_level", 0.0)
    y = compute_ed_pattern(g, indices, multiplicities, a, w_L, w_G; noise_level=noise_level)

    title = "$element-$structure"

    V = get(instrument, "voltage_kV", 200.0) * 1000.0
    λe = electron_wavelength(V)
    plot_title = "$title  (e⁻, λ=$(round(λe, digits=4)) Å)"

    theme(plot_theme)
    the_plot = plot(g, y, title=plot_title, xlabel="g = 1/d (1/Å)",
                    ylabel="Intensity (arb.)", show=false)

    return g, y, title, the_plot
end


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


# Phosphor-green colour ramp (black → dark green → bright green → highlight),
# the look of a fluorescent ED viewing screen. `gamma` < 1 lifts faint outer rings.
const PHOSPHOR_RAMP = ["#000000", "#022b06", "#1f9b3a", "#5dff7a", "#e6ffe9"]


"""
    radial_profile_value(y, g_min, g_max, gg)

Linear interpolation of the uniform g-grid profile `y` (over `[g_min, g_max]`) at
scattering vector `gg`. Clamps to the end samples outside the grid. Internal
helper for `render_ring_image`.
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
    render_ring_image(g, y, camera_constant; kwargs...) -> Plots.Plot

Render the 1D powder electron-diffraction profile `y(g)` as a 2D Debye–Scherrer
ring pattern. A powder pattern is rotationally symmetric, so the image is a pure
radial lookup: a pixel at distance r (mm) from the centre maps to g = r /
`camera_constant` (1/Å) and takes intensity `y(g)`. This is the SAED-style ring
image students measure — ring radius r = `camera_constant`·g, so r² ∝ N.

# Arguments
- `g::Vector{Float64}`     : profile g-grid (1/Å), uniform, from `do_it_electron`
- `y::Vector{Float64}`     : profile intensity at each g
- `camera_constant::Float64`: λL (mm·Å); ring radius r = camera_constant·g (mm)

# Keywords
- `image_px::Int=700`       : output image side length (pixels)
- `beam_stop_mm::Float64=2.5`: central beam-stop radius (mm); blanked to floor
- `phosphor::Bool=true`     : phosphor-green colormap (false → grayscale)
- `gamma::Float64=0.5`      : display gamma (<1 lifts faint outer rings)
- `noise_level::Float64=0.0`: per-pixel multiplicative noise (0–1), seeded upstream

Returns a square `Plots.Plot` heatmap (no axes/frame) ready to `savefig`.
"""
function render_ring_image(g::Vector{Float64},
                           y::Vector{Float64},
                           camera_constant::Float64;
                           image_px::Int = 700,
                           beam_stop_mm::Float64 = 2.5,
                           phosphor::Bool = true,
                           gamma::Float64 = 0.5,
                           noise_level::Float64 = 0.0
                           )::Plots.Plot
    length(g) == length(y) || throw(DimensionMismatch("g and y must have equal length"))
    camera_constant > 0 || throw(ArgumentError("camera_constant must be positive"))
    image_px ≥ 2 || throw(ArgumentError("image_px must be ≥ 2"))
    0 ≤ noise_level ≤ 1 || throw(DomainError(noise_level, "noise_level must be between 0 and 1"))

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

    cmap = phosphor ? cgrad(PHOSPHOR_RAMP) : cgrad(:grays)
    return heatmap(coords, coords, disp;
                   c = cmap, aspect_ratio = :equal, colorbar = false,
                   axis = false, ticks = false, framestyle = :none,
                   legend = false, grid = false, widen = false,
                   background_color = :black, margin = 0 * Plots.mm,
                   size = (image_px, image_px), clims = (0, 1), show = false)
end


"""
    do_it(file_name, structure, element, a, plot_theme)

Generate a complete XRD diffraction pattern for one (structure, element) sample.

Reads instrument parameters from a TOML config file, computes peak positions
and widths for the given structure and lattice parameter, adds background and
optional noise, and produces a plot.

# Arguments
- `file_name::String`: Path to TOML configuration file
- `structure::String`: Crystal structure ("SC", "BCC", or "FCC")
- `element::String`: Element label (used in plot title and filename)
- `a::Float64`: Lattice parameter in Angstroms
- `plot_theme::Symbol`: Plots.jl theme (e.g., `:dark`, `:light`)

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}, String, Plots.Plot}`:
  - 2θ angles in degrees
  - Intensity values
  - Plot/file title string (`"{element}-{structure}"`)
  - Plots.jl figure object
"""
function do_it(file_name::String,
               structure::String,
               element::String,
               a::Float64,
               plot_theme::Symbol
               )::Tuple{Vector{Float64}, Vector{Float64}, String, Plots.Plot}

    instrument, peak_width, _ = read_xrd_config(file_name)

    if get(instrument, "radiation", "xray") == "electron"
        return do_it_electron(instrument, peak_width, structure, element, a, plot_theme)
    end

    θ = collect(LinRange(instrument["two_theta_min"]/2,
                         instrument["two_theta_max"]/2,
                         instrument["N"]))
    λ = instrument["lambda"]

    max_hkl_sq = bragg_max_hkl_sq(a, λ)
    indices, multiplicities = Miller_indices(structure, max_hkl_sq)
    w_L, w_G = compute_peak_widths(θ, peak_width, instrument)

    noise_level = get(instrument, "noise_level", 0.0)
    y = compute_xrd_pattern(θ, indices, multiplicities, λ, a, w_L, w_G; noise_level=noise_level)

    title = "$element-$structure"

    theme(plot_theme)

    twoθ_deg = 2 * rad2deg.(θ)
    the_plot = plot(twoθ_deg, y, title=title, xlabel="2θ (deg)", ylabel="Intensity (arb.)", show=false)

    return twoθ_deg, y, title, the_plot
end

