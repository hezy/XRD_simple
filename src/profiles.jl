# Peak profiles (Voigt, pseudo-Voigt), their combined FWHM, and the sum of peaks.

"""
    abstract_peak(θ, θ₀, A, w_L, w_G; cutoff_sigma=5.0, normalize=false)

Template for peak profile functions (Voigt and Pseudo-Voigt).

Arguments:
- `θ::Vector{Float64}`: Position values where to evaluate the peak
- `θ₀::Float64`: Center position of the peak
- `A::Float64`: Peak amplitude (must be positive)
- `w_L::Float64`: Lorentzian full width at half maximum (FWHM) (must be positive)
- `w_G::Float64`: Gaussian full width at half maximum (FWHM) (must be positive)

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



"""
    pseudo_Voigt_peak(θ, θ₀, A, w_L, w_G; cutoff_sigma=5.0, normalize=false)

Computes the Thompson–Cox–Hastings pseudo-Voigt approximation of a Voigt peak:
a linear combination η·L + (1-η)·G of a Lorentzian and a Gaussian that both have
the combined FWHM f = `peak_fwhm(w_L, w_G)`.

See `abstract_peak` for parameter descriptions.

Notes:
- Mixing factor η is a cubic in w_L/f (Thompson, Cox & Hastings 1987)
- Implements bounds checking to improve performance for large datasets
"""
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
    
    # Both components share the combined FWHM f
    f = peak_fwhm(w_L, w_G)
    γ = f / 2                                      # Lorentzian HWHM
    σ = f / (2√(2log(2)))                          # Gaussian standard deviation
    
    # Mixing factor (Thompson-Cox-Hastings)
    q = w_L / f
    η = 1.36603 * q - 0.47719 * q^2 + 0.11116 * q^3
    
    # Initialize output array
    result = zeros(Float64, length(θ))
 
    # Calculate profile only for points within the cutoff region
    for i in eachindex(θ)
        if abs(θ[i] - θ₀) ≤ cutoff_sigma * f
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


# Utility functions that work with both Voigt and pseudo Voigt


"""
    peak_fwhm(w_L::Float64, w_G::Float64)

Calculates the full width at half maximum for either Voigt or pseudo-Voigt profile
(Olivero–Longbothum approximation, accurate to about 0.02 %).
"""
function peak_fwhm(w_L::Float64,
                   w_G::Float64
                   )::Float64
                   
    return 0.5346 * w_L + √(0.2166 * w_L^2 + w_G^2)
end


"""
    sum_peaks(x, x_list, multiplicities, w_L, w_G)

Sum pseudo-Voigt peak profiles at given peak centres, weighted by multiplicity.

Each entry in `θ_list` is one canonical reflection; its amplitude is the
multiplicity of that family. Summing one weighted peak per family is
mathematically identical to summing every sign+permutation variant at unit
amplitude, and far cheaper.

# Arguments
- `x::Vector{Float64}`: Grid on which the pattern is evaluated (2θ in radians, or g)
- `x_list::Vector{Float64}`: Peak centre positions, same unit as `x`
- `multiplicities::Vector{Int}`: Multiplicity of each reflection family
- `w_L::Vector{Float64}`: Lorentzian FWHM of each peak, evaluated at its centre
- `w_G::Vector{Float64}`: Gaussian FWHM of each peak, evaluated at its centre

# Returns
- `Vector{Float64}`: Combined peak intensities at each x
"""
function sum_peaks(x::Vector{Float64},
                   x_list::Vector{Float64},
                   multiplicities::Vector{Int},
                   w_L::Vector{Float64},
                   w_G::Vector{Float64},
                   )::Vector{Float64}

    length(x_list) == length(multiplicities) == length(w_L) == length(w_G) ||
        throw(DimensionMismatch(
            "x_list, multiplicities, w_L and w_G must have same length"))

    y = zeros(length(x))
    for i in eachindex(x_list)
        y .+= pseudo_Voigt_peak(x, x_list[i], Float64(multiplicities[i]), w_L[i], w_G[i])
    end
    return y
end
