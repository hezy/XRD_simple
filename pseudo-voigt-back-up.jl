
using Plots
using SpecialFunctions

x = range(start=-3, stop=3, step=0.01)

function Lorentzian(x, fwhm)
    γ = fwhm / 2
    return @. (γ/pi) / (x^2 + γ^2)
end


function Gaussian(x, fwhm)
    σ = fwhm/(2√(2log(2)))
    return @. 1/√(2π)/σ * exp(-x^2/2σ^2)
end


function pseudo_Voigt(x, fwhm, n)
	return n * Lorentzian(x, fwhm) + (1 - n) * Gaussian(x, fwhm)
end


function Voigt(x, fwhm_l, fwhm_g)
    γ = fwhm_l/2
    σ = fwhm_g/(2√(2log(2)))
    z = @. -im * (x + im * γ) / (√2 * σ)
    return @. real(erfcx(z)) / (√(2pi) * σ)
end

#%%
y1 = Lorentzian(x, 1)
y2 = Gaussian(x, 1)
y3 = pseudo_Voigt(x, 1, 0.5)
y4 = Voigt(x, 0.45, 0.72)
#%%
p = plot(x, [y1 y2 y3 y4],
         label=["Lorentzian" "Gaussian" "Pseudo Voigt" "Voigt"])
title!("peak functions")
xlabel!(raw"x")
ylabel!(raw"y")
##%

#display(p)

#savefig(plot,"file.png")
