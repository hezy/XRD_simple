"""
simple_XRD_new.jl
by Hezy Amiel
April 2023
Julia 1.8.5
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
- `w_L::Vector{Float64}`: Lorentzian full width at half maximum (FWHM) (must be positive)
- `w_G::Vector{Float64}`: Gaussian full width at half maximum (FWHM) (must be positive)

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
# Scaler w's version
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
    w_L .> 0 || throw(ArgumentError("Lorentzian width w_L must be positive"))
    w_G .> 0 || throw(ArgumentError("Gaussian width w_G must be positive"))
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
# Scaler w's version
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
    w_L .> 0 || throw(ArgumentError("Lorentzian width w_L must be positive"))
    w_G .> 0 || throw(ArgumentError("Gaussian width w_G must be positive"))
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
            L = γ[i] / (π * ((θ[i] - θ₀)^2 + γ^2))
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
    estimate_peak_bounds(w_L::Float64, w_G::Float64, tol::Float64=1e-6)
    estimate_peak_bounds(w_L::Vector{Float64}, w_G::Vector{Float64}, tol::Float64=1e-6)

Estimates the distance from peak center where profile falls below a given tolerance.
Works for both Voigt and pseudo-Voigt profiles. Returns maximum bound for vector inputs.
"""
# Scalar w's version
function estimate_peak_bounds(w_L::Float64,
                              w_G::Float64;
                              tol::Float64=1e-6
                              )::Float64
    γ = w_L / 2
    σ = w_G / (2√(2log(2)))
    
    gaussian_cutoff = √(-2 * log(tol))
    lorentzian_cutoff = γ * √(1/tol - 1) / σ
    
    return max(gaussian_cutoff, lorentzian_cutoff)
end

# Vector w's version
function estimate_peak_bounds(w_L::Vector{Float64},
                              w_G::Vector{Float64};
                              tol::Float64=1e-6
                              )::Float64
    length(w_L) == length(w_G) || throw(DimensionMismatch("w_L and w_G must have same length"))
    
    γ = w_L ./ 2
    σ = w_G ./ (2√(2log(2)))
    
    gaussian_cutoffs = fill(√(-2 * log(tol)), length(w_L))
    lorentzian_cutoffs = @. γ * √(1/tol - 1) / σ
    
    # Return the maximum bound to ensure coverage of entire peak
    return maximum(max.(gaussian_cutoffs, lorentzian_cutoffs))
end



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
Calculates the width of the Gaussian as a function θ with U, V, W parameters -
Caglioti formula
"""
function Gaussian_peaks_width(θ::Vector{Float64},
                              U::Float64,
                              V::Float64,
                              W::Float64
                              )::Vector{Float64}

        return @. √(U * tan(θ)^2 + V * tan(θ) + W)
end


"""
Calculates the width of the Lorntzian as a function θ with K, E, λ, D parameters -
Strain (Stokes-Wilson) and size (Scherrer) broadening
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
    sinθ_cleaned = [item for item in sinθ if abs(item) <= 1]  # removing values outside (-1,1)
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


"sums peak functions to return intensity vs angle"
function sum_peaks(θ::Vector{Float64},
                   θ_list::Vector{Float64},
                   w_L::Vector{Float64},
                   w_G::Vector{Float64},
                   )::Vector{Float64}
    
    y = zeros(size(θ))
    # cutoff = estimate_peak_bounds(w_L, w_G; tol=1e-6)
    for item in θ_list
        y = y + pseudo_Voigt_peak(θ, item, 1.0, w_L, w_G) #; cutoff_sigma=cutoff, normalize=true)
        # y = y + Voigt_peak(θ, item, 1.0, w_L, w_G; cutoff_sigma=cutoff, normalize=true)
    end
    return y
end


"""
   intensity_vs_angle(θ::Vector{Float64}, indices::Vector{Vector{Int}}, λ::Float64, 
                     a::Float64, w_L::Vector{Float64}, w_G::Vector{Float64})::Vector{Float64}

Calculate X-ray diffraction pattern by summing peak profiles at allowed Bragg angles.

# Arguments
- `θ::Vector{Float64}`: Scattering angles for intensity calculation (radians)
- `indices::Vector{Vector{Int}}`: Miller indices of crystal planes
- `λ::Float64`: X-ray wavelength (Å)
- `a::Float64`: Lattice parameter (Å)
- `w_L::Vector{Float64}`: Lorentzian width parameters
- `w_G::Vector{Float64}`: Gaussian width parameters

# Returns
- `Vector{Float64}`: XRD intensities at each θ angle

# Throws
- `ArgumentError`: If λ ≤ 0, a ≤ 0, any width ≤ 0, or w_L and w_G have different lengths
"""
function intensity_vs_angle(θ::Vector{Float64},
                         indices::Vector{Vector{Int}},
                         λ::Float64,
                         a::Float64,
                         w_L::Vector{Float64},
                         w_G::Vector{Float64}
                         )::Vector{Float64}
   
   λ <= 0 && throw(ArgumentError("Wavelength must be positive"))
   a <= 0 && throw(ArgumentError("Lattice parameter must be positive"))
   length(w_L) != length(w_G) && throw(ArgumentError("Width parameter vectors must have same length"))
   any(w_L .<= 0) && throw(ArgumentError("Lorentzian widths must be positive"))
   any(w_G .<= 0) && throw(ArgumentError("Gaussian widths must be positive"))

   θ_list, _ = bragg_angles(λ, d_list(indices, a))
   y = sum_peaks(θ, θ_list, w_L, w_G)
   return y
end


"Returns a list of Miller indices for each one of the cubic symmetries"
function Miller_indices(cell_type::String,
                        min::Int, 
                        max::Int
                        )::Vector{Vector{Int}}
    
    if !(cell_type in ["SC", "BCC", "FCC"])
        error("Invalid cell_type: $cell_type. Expected 'SC', 'BCC', or 'FCC'.")
    end
    if min > max
        error("Minimum value cannot be greater than maximum value.")
    end
    if !(isa(min, Int) && isa(max, Int))
        error("Minimum and maximum values must be integers.")
    end

    if cell_type == "SC"
        # In simple cubic lattice, all Miller indices are allowed
        return [
            [h, k, l] for h = min:max for k = min:max for l = min:max if [h, k, l] != [0, 0, 0]
        ]

    elseif cell_type == "BCC"
        # In body centered cubic lattice, only indices with h+k+l=even are allowed
        return [
            [h, k, l] for h = min:max for k = min:max for l = min:max if iseven(h + k + l) && [h, k, l] != [0, 0, 0]
        ]

    elseif cell_type == "FCC"
        # In face centered cubic lattice, h,k,l must all be either odd or even
        return [
            [h, k, l] for h = min:max for k = min:max for l = min:max if
            ((iseven(h) && iseven(k) && iseven(l)) || (isodd(h) && isodd(k) && isodd(l))) &&
            [h, k, l] != [0, 0, 0]
        ]

    end

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
    air_scatter = @. 50 * exp(-5θ)        # Strong at low angles
    fluorescence = 10.0                    # Constant background
    base = air_scatter .+ fluorescence
    
    # Add optional noise
    if noise_level > 0
        noise = noise_level * randn(length(θ))
        return max.(base .+ noise, 0)  # Ensure non-negative intensity
    else
        return base
    end
end


"Adding some noise to the data"
function make_noisy(θ::Vector{Float64},
                    y::Vector{Float64}
                    )::Vector{Float64}
    
    return (background(θ) + y) .* 0.5 * rand(Normal(1, 0.1), size(θ))
end



"""
        read_xrd_config(filename::String) -> 
        (Dict{String,Any}, Dict{String,Float64}, Dict{String,Tuple{String,Float64}})

Read XRD configuration from TOML file, returning instrument, peak width, and lattice parameters.

# Arguments
- `filename`: Path to TOML configuration file

# Returns
Tuple with:
- "instrument": Dict of instrument parameters (two_theta_min, two_theta_max, N, lambda)
- "peak_width": Dict of peak width parameters (U, V, W, K, Epsilon, D)
- "lattice": Dict of lattice parameters by structure type

Note: Angular parameters (two_theta_min, two_theta_max) are automatically converted to radians.
"""
function read_xrd_config(filename::String)
    config = TOML.parsefile(filename)
    
    instrument = Dict{String,Any}(
        k => (k in ["two_theta_min", "two_theta_max"] ? deg2rad(v) : v)
        for (k,v) in config["instrument"]
    )
    
    peak_width = Dict{String,Float64}(config["peak_width"])
    
    lattice = Dict{String,Tuple{String,Float64}}()
    for (structure, elements) in config["lattice"]
        for (element, value) in elements
            lattice[structure] = (element, value)
        end
    end
    
    return instrument, peak_width, lattice
end



"colecting input data, building the XRD pattern with background and noise, plotting it"
function do_it_zero(file_name::String
                    )::Vector{Float64}
    
    instrument_data, peak_width, lattice_params = read_xrd_config(file_name)
    θ = collect(LinRange((instrument["two_theta_min"]/2),
                         (instrument["two_theta_max"]/2),
                         instrument["N"]))
    return θ
end


"colecting input data, building the XRD pattern with background and noise, plotting it"
function do_it(file_name::String,
               lattice_type::String,
               plot_theme::Symbol
               )::Tuple{Vector{Float64}, Vector{Float64}, String, Plots.Plot}
    
    instrument, peak_width, lattice = read_xrd_config(file_name)

    N = instrument["N"]
    θ = collect(LinRange((instrument["two_theta_min"]/2),
                         (instrument["two_theta_max"]/2),
                         instrument["N"]))
    y = zeros(instrument["N"])
    λ = instrument["lambda"]
    U, V, W = peak_width["U"], peak_width["V"], peak_width["W"]
    K, ϵ, D = peak_width["K"], peak_width["Epsilon"], peak_width["D"]
    a = lattice[lattice_type][2]

    index_min::Int = -5
    index_max::Int = 5
    indices = Miller_indices(lattice_type, index_min, index_max)

    w_L = Lorentzian_peaks_width(θ, K, ϵ, λ, D)
    w_G = Gaussian_peaks_width(θ, U, V, W)
    y = (background(θ) +
         intensity_vs_angle(θ, indices, λ, a, w_L, w_G)) # .* 0.01 .* rand(Normal(1, 0.1), N)

    the_title = "XRD - " * lattice_type

    theme(plot_theme)

    twoθ_deg = 2 * rad2deg.(θ)
    the_plot = plot(twoθ_deg, y, title=the_title, xlabel="2θ (deg)", ylabel="Intensity (arb.)")

    return twoθ_deg, y, the_title, the_plot
end

