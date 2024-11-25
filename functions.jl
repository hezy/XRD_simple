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



# Utility functions that work with both peak types


"""
    estimate_peak_bounds(w_L::Float64, w_G::Float64, tol::Float64=1e-6)
    estimate_peak_bounds(w_L::Vector{Float64}, w_G::Vector{Float64}, tol::Float64=1e-6)

Estimates the distance from peak center where profile falls below a given tolerance.
Works for both Voigt and pseudo-Voigt profiles. Returns maximum bound for vector inputs.
"""
# Scalar version
function estimate_peak_bounds(w_L::Float64, w_G::Float64, tol::Float64=1e-6)::Float64
    γ = w_L / 2
    σ = w_G / (2√(2log(2)))
    
    gaussian_cutoff = √(-2 * log(tol))
    lorentzian_cutoff = γ * √(1/tol - 1) / σ
    
    return max(gaussian_cutoff, lorentzian_cutoff)
end

# Vector version
function estimate_peak_bounds(w_L::Vector{Float64},
                              w_G::Vector{Float64};
                              tol::Float64=1e-5
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
Calculates the width of the Lorntzian as a function θ with K, ϵ, λ, D parameters -
Strain (Stokes-Wilson) and size (Scherrer) broadening
"""
function Lorentzian_peaks_width(θ::Vector{Float64},
                                K::Float64,
                                ϵ::Float64,
                                λ::Float64,
                                D::Float64,
                                )::Vector{Float64}

    # Strain broadening (Stokes-Wilson)
    w_L_strain = @. 4 * ε * tan(θ)
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



"Returnes the inter-layers distances as a function of Miller_indices"
function d_list(indices::Vector{Vector{Int}},
                a::Float64
                )::Vector{Float64}

    return a ./ [√(i^2 + j^2 + k^2) for (i, j, k) in indices]
end


"sums peak functions to return intensity vs angle"
function sum_peaks(θ::Vector{Float64},
                   θ_list::Vector{Float64},
                   w_L::Vector{Float64},
                   w_G::Vector{Float64},
                   )::Vector{Float64}
    
    y = zeros(size(θ))
    cutoff = estimate_peak_bounds(w_L, w_G; tol=1e-6)
    for item in θ_list
        # y = y + pseudo_Voigt_peak(θ, item, 1.0, w_L, w_G; cutoff_sigma=cutoff, normalize=true)
        y = y + Voigt_peak(θ, item, 1.0, w_L, w_G; cutoff_sigma=cutoff, normalize=true)
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


"background function for the XRD pattern"
function background(θ::Vector{Float64}
                    )::Vector{Float64}

    return @. 2 + 0.5 * θ * (π - θ)       
end


"Adding some noise to the data"
function make_noisy(θ::Vector{Float64},
                    y::Vector{Float64}
                    )::Vector{Float64}
    
    return (background(θ) + y) .* 0.5 * rand(Normal(1, 0.1), size(θ))
end


"Reading a text file with instrument data, and lattice parameters"
function read_file(filename::String
                   )::Tuple{Dict,Dict}
    
    instrument_data = Dict{AbstractString,Any}()
    lattice_params = Dict{AbstractString,Float64}()

    # read file line by line
    for line in eachline(filename)
        # split the line by whitespace and remove empty strings
        tokens = filter(x -> x ≠ "", split(line))

        if length(tokens) > 0 && tokens[1] ≠ "#"
            if tokens[1] in ["2θ_min", "2θ_max"]
                instrument_data[tokens[1]] = deg2rad(parse(Float64, tokens[2]))
            elseif tokens[1] == "N"
                instrument_data[tokens[1]] = parse(Int64, tokens[2])    
            elseif tokens[1] == "λ"
                instrument_data[tokens[1]] = parse(Float64, tokens[2])
            elseif tokens[1] in ["U", "V", "W"]
                instrument_data[tokens[1]] = parse(Float64, tokens[2])
            elseif tokens[1] in ["K", "ϵ", "D"]
                instrument_data[tokens[1]] = parse(Float64, tokens[2])
            elseif tokens[1] in ["SC", "BCC", "FCC"]
                lattice_params[tokens[1]] = parse(Float64, tokens[3])
            end
        elseif length(tokens) > 1 && tokens[1] == "BCC" && tokens[2] ≠ "#"
            lattice_params[tokens[1]] = parse(Float64, tokens[3])
        elseif length(tokens) > 1 && tokens[1] == "FCC" && tokens[2] ≠ "#"
            lattice_params[tokens[1]] = parse(Float64, tokens[3])
        elseif length(tokens) > 1 && tokens[1] == "SC" && tokens[2] ≠ "#"
            lattice_params[tokens[1]] = parse(Float64, tokens[3])
        end
    end

    return instrument_data, lattice_params
end


"colecting input data, building the XRD pattern with background and noise, plotting it"
function do_it_zero(file_name::String
                    )::Vector{Float64}
    
    instrument_data, lattice_params = read_file(file_name)
    θ = collect(LinRange((instrument_data["2θ_min"]/2),
                         (instrument_data["2θ_max"]/2),
                         instrument_data["N"]))
    return θ
end


"colecting input data, building the XRD pattern with background and noise, plotting it"
function do_it(file_name::String,
               lattice_type::String,
               plot_theme::Symbol
               )::Tuple{Vector{Float64}, Vector{Float64}, String, Plots.Plot}
    
    instrument_data, lattice_params = read_file(file_name)

    N = instrument_data["N"]
    θ = collect(LinRange((instrument_data["2θ_min"]/2),
                         (instrument_data["2θ_max"]/2),
                         instrument_data["N"]))
    y = zeros(instrument_data["N"])
    λ = instrument_data["λ"]
    U, V, W = instrument_data["U"], instrument_data["V"], instrument_data["W"]
    K, ϵ, D = instrument_data["K"], instrument_data["ϵ"], instrument_data["D"]
    a = lattice_params[lattice_type]

    index_min::Int = -5
    index_max::Int = 5
    indices = Miller_indices(lattice_type, index_min, index_max)

    w_L = Lorentzian_peaks_width(θ, K, ϵ, λ, D)
    w_G = Gaussian_peaks_width(θ, U, V, W)
    y = (background(θ) +
         intensity_vs_angle(θ, indices, λ, a, w_L, w_G)) .*
        0.1 .* rand(Normal(1, 0.1), N)

    the_title = "XRD - " * lattice_type

    theme(plot_theme)

    twoθ_deg = 2 * rad2deg.(θ)
    the_plot = plot(twoθ_deg, y, title=the_title, xlabel="2θ (deg)", ylabel="Intensity (arb.)")

    return twoθ_deg, y, the_title, the_plot
end

