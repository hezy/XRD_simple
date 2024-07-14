function Voigt_peak(θ::Vector, θ₀, A, w_L, w_G, n)
    """Returns a Voigt peak centered around θ₀, with amplitude A, width w, and mixing factor n """
    """untested"""
    γ = w_L / 2
    σ = w_G / (2√(2log(2)))
    z = @. -im * ((θ - θ₀) + im * γ) / (√2 * σ)
    return @. real(erfcx(z)) / (√(2pi) * σ)
end

