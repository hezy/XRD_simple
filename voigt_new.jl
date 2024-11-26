"""
Refactored code with repeated elements isolated into separate helper functions:
"""


function peak_fwhm(w_L, w_G)
    return 0.5346 * w_L + √(0.2166 * w_L^2 + w_G^2)
end


function validate_parameters(A, w_L, w_G, cutoff_sigma)
    A > 0 || throw(ArgumentError("Amplitude A must be positive"))
    all(w_L .> 0) || throw(ArgumentError("Lorentzian width w_L must be positive"))
    all(w_G .> 0) || throw(ArgumentError("Gaussian width w_G must be positive"))
    cutoff_sigma > 0 || throw(ArgumentError("cutoff_sigma must be positive"))
end


function normalize_result(result)
    maxval = maximum(result)
    if maxval > 0
        result ./= maxval
    end
    return result
end


function Voigt_peak(θ::Vector{Float64},
                    θ₀::Float64,
                    A::Float64,
                    w_L::Union{Float64,Vector{Float64}},
                    w_G::Union{Float64,Vector{Float64}};
                    cutoff_sigma::Float64=5.0,
                    normalize::Bool=false
                    )::Vector{Float64}

    validate_parameters(A, w_L, w_G, cutoff_sigma)
    
    result = zeros(Float64, length(θ))                   
    γ = w_L / 2                                    
    σ = w_G / (2√(2log(2)))                        
    w_eff = peak_fwhm(w_L, w_G)    
    
    for i in eachindex(θ)
        if abs(θ[i] - θ₀) ≤ cutoff_sigma * (isa(w_eff, Number) ? w_eff : w_eff[i])
            z = -im * (θ[i] - θ₀ + im * (isa(γ, Number) ? γ : γ[i])) / (√2 * (isa(σ, Number) ? σ : σ[i]))   
            result[i] = A * real(erfcx(z)) / (√(2π) * (isa(σ, Number) ? σ : σ[i]))
        end
    end

    if normalize
        result = normalize_result(result)
    end
    
    return result
end


function pseudo_Voigt_peak(θ::Vector{Float64},
                         θ₀::Float64,
                         A::Float64,
                         w_L::Union{Float64,Vector{Float64}},
                         w_G::Union{Float64,Vector{Float64}};
                         cutoff_sigma::Float64=5.0,
                         normalize::Bool=false
                         )::Vector{Float64}
    
    validate_parameters(A, w_L, w_G, cutoff_sigma)
    
    γ = w_L / 2
    σ = w_G / (2√(2log(2)))
    w_eff = peak_fwhm(w_L, w_G)    
    η = @. 1.36603 * (w_L/w_eff) - 0.47719 * (w_L/w_eff)^2 + 0.11116 * (w_L/w_eff)^3
    
    result = zeros(Float64, length(θ))
 
    for i in eachindex(θ)
        if abs(θ[i] - θ₀) ≤ cutoff_sigma * (isa(w_eff, Number) ? w_eff : w_eff[i])
            L = (isa(γ, Number) ? γ : γ[i]) / (π * ((θ[i] - θ₀)^2 + (isa(γ, Number) ? γ : γ[i])^2))
            G = exp(-(θ[i] - θ₀)^2 / (2(isa(σ, Number) ? σ : σ[i])^2)) / ((isa(σ, Number) ? σ : σ[i]) * √(2π))
            result[i] = A * ((isa(η, Number) ? η : η[i]) * L + (1 - (isa(η, Number) ? η : η[i])) * G)
        end
    end
    
    if normalize
        result = normalize_result(result)
    end
    
    return result
end
