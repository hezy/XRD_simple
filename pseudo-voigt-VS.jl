### A Pluto.jl notebook ###
# v0.19.19

using Plots
using SpecialFunctions

x = range(start=-2.5, stop=2.5, step=0.01)

function Gaussian(x, fwhm)
	g1 = 2.0 * sqrt(log(2)/π)
	g2 = 4.0 * log(2)
    return @. (g1 / fwhm) * exp((-g2 * x^2) / (fwhm^2)) 
end

function Lorentzian(x, fwhm)
	l1 = 2/π
	l2 = 4.0
    return @. (l1 / fwhm) / ((1 + l2 * x^2) / fwhm^2 ) 
end

function pseudo_Voigt(x, fwhm, n)
	return n * Gaussian(x, fwhm) + (1 - n) * Lorentzian(x, fwhm)
end

function Voigt(x, gamma, sigma)
    z = @. -1im * (x + 1im * gamma) / (sqrt(2)*sigma)
    return @. real(erfcx(z))/(sqrt(2pi)*sigma)
end
#erfcx

y1 = Lorentzian(x, 1)
y2 = Gaussian(x, 1)
y3 = pseudo_Voigt(x, 1, 0.5)
y4 = Voigt(x, 0.25, 0.25)

p = plot(x,[y1 y2 y3 y4], label=["Lorentzian" "Gaussian" "Pseudo Voigt" "Voigt"])
title!("peak functions")
xlabel!("x")
ylabel!("y")

display(p)