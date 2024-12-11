using SpecialFunctions
using Plots

function voigt(x, amplitude, center, sigma, gamma)
    # Proper normalization for Voigt profile
    z = (x .- center .+ im*gamma) ./ (sigma * sqrt(2))
    return amplitude * real.(erfcx.(-im*z)) / (sqrt(π))
end

function pseudo_voigt(x, amplitude, center, w_G, w_L, eta=nothing)
    sigma = w_G / (2 * sqrt(2 * log(2)))
    
    if isnothing(eta)
        f_G = w_G^5
        f_L = w_L^5
        f = (f_G + 2.69269 * f_G * f_L + 2.42843 * f_L * f_G +
             4.47163 * f_L * f_G^2 + 0.07842 * f_L^2 * f_G +
             f_L^2)^(1/5)
        eta = 1.36603 * (w_L/f) - 0.47719 * (w_L/f)^2 + 0.11116 * (w_L/f)^3
    end
    
    gaussian = amplitude * exp.(-log(2) * ((x .- center)/(w_G/2)).^2)
    lorentzian = amplitude ./ (1 .+ ((x .- center)/(w_L/2)).^2)
    
    return eta * lorentzian + (1 - eta) * gaussian
end

function compare_profiles(w_G, w_L; x_range=(-10, 10), points=1000)
    x = range(x_range[1], x_range[2], length=points)
    
    sigma = w_G / (2 * sqrt(2 * log(2)))
    gamma = w_L / 2
    
    v = voigt(x, 1, 0, sigma, gamma)
    pv = pseudo_voigt(x, 1, 0, w_G, w_L)
    
    p = plot(x, v, label="True Voigt", linewidth=2)
    plot!(x, pv, label="Pseudo-Voigt", linestyle=:dash, linewidth=2)
    xlabel!("x")
    ylabel!("Intensity")
    title!("Comparison (w_G=$(w_G), w_L=$(w_L))")
    
    max_diff = maximum(abs.(v .- pv))
    println("Maximum absolute difference: ", max_diff)
    
    return x, v, pv, p
end


# example usage:
x, v, pv, p = compare_profiles(2.0, 2.0)
display(p)
