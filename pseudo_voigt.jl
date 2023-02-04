"""
pseudo-voigt.jl
by Hezy Amiel
January 2023
Julia 1.8.5

https://en.m.wikipedia.org/wiki/Voigt_profile
http://journals.iucr.org/j/issues/1997/04/00/gl0484/gl0484.pdf
http://journals.iucr.org/j/issues/2000/06/00/nt0146/nt0146.pdf
https://www.onlinelibrary.wiley.com/doi/epdf/10.1002/sia.5521
"""

#using PyPlot
using Plots
using SpecialFunctions


function Lorentzian(x, fwhm)
    γ = fwhm / 2
    return @. (γ / pi) / (x^2 + γ^2)
end


function Gaussian(x, fwhm)
    σ = fwhm / (2√(2log(2)))
    return @. 1 / √(2π) / σ * exp(-x^2 / 2σ^2)
end


mix_functions(f1, f2, n) = n * f1 + (1 - n) * f2


pseudo_Voigt(x, fwhm, n) =  mix_functions(Lorentzian(x, fwhm), Gaussian(x, fwhm), n)


function Voigt(x, fwhm_L, fwhm_G)
    γ = fwhm_L / 2
    σ = fwhm_G / (2√(2log(2)))
    z = @. -im * (x + im * γ) / (√2 * σ)
    return @. real(erfcx(z)) / (√(2pi) * σ)
end


x = range(start = -3, stop = 3, step = 0.01)

y1 = Lorentzian(x, 1)
y2 = Gaussian(x, 1)
y3 = pseudo_Voigt(x, 1, 0.5)
y4 = Voigt(x, 0.45, 0.72)

plot1 = plot(x, [y1 y2 y3 y4], label = ["Lorentzian" "Gaussian" "Pseudo Voigt" "Voigt"])
title!("peak functions")
xlabel!(raw"x")
ylabel!(raw"y")

display(plot1)
savefig(plot1, "psedo_voigt")
println("The End")
