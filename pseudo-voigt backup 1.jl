### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils
using Plots
using Ranges

# ╔═╡ ab93661e-92ba-11ed-0db0-a5b6da440669
function Lorentzian(x, fwhm)
	l1 = 2/π
	l2 = 4.0
    return @. (l1 / fwhm) / ((1 + l2 * x^2) / fwhm^2 ) 
end

# ╔═╡ c1159303-2800-460e-a6ff-e957681e4c81
function Gaussian(x, fwhm)
	g1 = 2.0 * sqrt(log(2)/π)
	g2 = 4.0 * log(2)
    return @. (g1 / fwhm) * exp((-g2 * x^2) / (fwhm^2)) 
end

# ╔═╡ f712027f-7f6e-4082-95f3-cae200d7f18a
function pseudo_Voigt(x, fwhm, n)
	return n * Gaussian(x, fwhm) + (1 - n) * Lorentzian(x, fwhm)
end

# ╔═╡ 298367a4-6f27-42dd-a115-e3235371347c
x = range(start=-5, stop=5, step=0.01)

# ╔═╡ 2c83b7e1-5edf-441b-b0c2-210c2927da70
y1 = Lorentzian(x, 1)

# ╔═╡ 8741d4e6-8a32-453b-8e27-f466747628af
y3 = pseudo_Voigt(x, 1, 0.5)


# ╔═╡ 3be58c8e-a791-48fb-ba4e-3913e0247c1a
begin
	plot(x,y1, label="Lorentzian")
	plot!(x, y2, label="Gaussian")
	plot!(x, y3, label="Pseudo Voigt")
end


