# Peak profiles (Voigt, pseudo-Voigt), their combined FWHM, the Debye–Waller
# factor shared by both modes, and the sum of peaks.

"""
    Voigt_peak(θ, θ₀, A, w_L, w_G; cutoff_sigma=5.0, normalize=false)

Computes Voigt peak profile as the convolution of Gaussian and Lorentzian functions
using the complex error function.

# Arguments
- `θ::AbstractVector{<:Real}`: Position values where to evaluate the peak
- `θ₀::Real`: Center position of the peak
- `A::Real`: Peak area (must be positive)
- `w_L::Real`: Lorentzian full width at half maximum (FWHM) (must be positive)
- `w_G::Real`: Gaussian full width at half maximum (FWHM) (must be positive)

# Keyword Arguments
- `cutoff_sigma::Real=5.0`: The profile is zero farther than `cutoff_sigma`
  combined FWHMs from θ₀
- `normalize::Bool=false`: If true, normalize peak height to 1.0

# Returns
- `Vector{Float64}`: Peak intensity at each θ position

Notes:
- Uses the scaled complementary error function (erfcx) for numerical stability
- More computationally expensive but more accurate than pseudo-Voigt approximation
- Implements bounds checking to improve performance for large datasets
- The cutoff region is based on both Gaussian and Lorentzian widths
"""
function Voigt_peak(θ::AbstractVector{<:Real},
                    θ₀::Real,
                    A::Real,
                    w_L::Real,
                    w_G::Real;
                    cutoff_sigma::Real=5.0,
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

Arguments and return value as for `Voigt_peak`.

Notes:
- Mixing factor η is a cubic in w_L/f (Thompson, Cox & Hastings 1987)
- Implements bounds checking to improve performance for large datasets
"""
function pseudo_Voigt_peak(θ::AbstractVector{<:Real},
                           θ₀::Real,
                           A::Real,
                           w_L::Real,
                           w_G::Real;
                           cutoff_sigma::Real=5.0,
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
    peak_fwhm(w_L::Real, w_G::Real)

Calculates the full width at half maximum for either Voigt or pseudo-Voigt profile
(Olivero–Longbothum approximation, accurate to about 0.02 %).
"""
function peak_fwhm(w_L::Real,
                   w_G::Real
                   )::Float64
    return 0.5346 * w_L + √(0.2166 * w_L^2 + w_G^2)
end


"""
    Debye_Waller(s::Real, B::Real)::Float64

Debye–Waller (thermal) factor of the intensity of a reflection:

    exp(−2B s²),   s = sin θ / λ = 1 / (2d) = g / 2

where B = 8π²⟨u²⟩ is the parameter of the atomic temperature factor
exp(−B s²) on the scattering amplitude. It damps high-angle reflections.

# Arguments
- `s::Real`: sin θ / λ in 1/Å
- `B::Real`: Debye–Waller parameter in Å² (must not be negative)

# Returns
- `Float64`: the factor, in (0, 1]

# Throws
* ArgumentError: If B < 0
"""
function Debye_Waller(s::Real, B::Real)::Float64
    B ≥ 0 || throw(ArgumentError("Debye–Waller parameter B must not be negative"))
    return exp(-2B * s^2)
end


"""
    sum_peaks(x, x_list, amplitudes, w_L, w_G)

Sum pseudo-Voigt peak profiles at given peak centres, each with its own area.

Each entry in `x_list` is one canonical reflection; its area is the
multiplicity of that family times any intensity weight (see `peak_weights`).
Summing one weighted peak per family is mathematically identical to summing
every sign+permutation variant separately, and far cheaper.

# Arguments
- `x::AbstractVector{<:Real}`: Grid on which the pattern is evaluated (2θ in radians, or g)
- `x_list::AbstractVector{<:Real}`: Peak centre positions, same unit as `x`
- `amplitudes::AbstractVector{<:Real}`: Area of each peak
- `w_L::AbstractVector{<:Real}`: Lorentzian FWHM of each peak, evaluated at its centre
- `w_G::AbstractVector{<:Real}`: Gaussian FWHM of each peak, evaluated at its centre

# Returns
- `Vector{Float64}`: Combined peak intensities at each x
"""
function sum_peaks(x::AbstractVector{<:Real},
                   x_list::AbstractVector{<:Real},
                   amplitudes::AbstractVector{<:Real},
                   w_L::AbstractVector{<:Real},
                   w_G::AbstractVector{<:Real},
                   )::Vector{Float64}

    length(x_list) == length(amplitudes) == length(w_L) == length(w_G) ||
        throw(DimensionMismatch(
            "x_list, amplitudes, w_L and w_G must have same length"))

    y = zeros(length(x))
    for i in eachindex(x_list)
        y .+= pseudo_Voigt_peak(x, x_list[i], amplitudes[i], w_L[i], w_G[i])
    end
    return y
end
