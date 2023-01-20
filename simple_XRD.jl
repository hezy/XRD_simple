"""
simple_XRD.jl
by Hezy Amiel
January 2023
Julia 1.8.5
"""

using Plots
using SpecialFunctions

x = range(start=-3, stop=3, step=0.01)

function Gaussian(x, fwhm)
    σ = fwhm/(2√(2log(2)))
    return @. 1/√(2π)/σ * exp(-x^2/2σ^2)
end


function Lorentzian(x, fwhm)
    γ = fwhm / 2
    return @. (γ/pi) / (x^2 + γ^2)
end


function Pseudo_Voigt(x, fwhm, n)
	return n * Lorentzian(x, fwhm) + (1 - n) * Gaussian(x, fwhm)
end


function Voigt(x, fwhm_L, fwhm_G)
    γ = fwhm_L/2
    σ = fwhm_G/(2√(2log(2)))
    z = @. -im * (x + im * γ) / (√2 * σ)
    return @. real(erfcx(z)) / (√(2pi) * σ)
end


function peak(x, x0, A, w, n)
    return A * Pseodo_Voigt(x-x0, w, n)
end


function intensity(theta_space, peaks_positions, peaks_width)
    y = []
    for peaks_positions
        #print(n, peaks_positions[n], peaks_width[n])
        y = y + peak(theta_space, peaks_positions[n], 1, peaks_width[n], 0.5)
    return y


y1 = Lorentzian(x, 1)
y2 = Gaussian(x, 1)
y3 = pseudo_Voigt(x, 1, 0.5)
y4 = Voigt(x, 0.45, 0.72)

p = plot(x,[y1 y2 y3 y4], label=["Lorentzian" "Gaussian" "Pseudo Voigt" "Voigt"])
title!("peak functions")
xlabel!(raw"x")
ylabel!(raw"y")



indices = [[h,k,l] for h=-2:2 for k=-2:2 for l=-2:2]
deleteat!(indices, findall(x->x==[0,0,0],indices))
print(indices)