### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 2dada36d-12c5-41cc-a941-5da6f45a8e39


# ╔═╡ ab93661e-92ba-11ed-0db0-a5b6da440669
function Lorentzian(x, fwhm)
	l1 = 2/π
	l2 = 4
    return (l/fwhm) / ((1 + l2*x^2) / fwhm^2) 
end

# ╔═╡ c1159303-2800-460e-a6ff-e957681e4c81
function Gaussian(x, position, width)
	g1 = 2*sqrt(log(2)/π)
	g2 = 4*log(2)
    return (g1/fwhm)exp((-g2 * x^2) / (fwhm^2)) 
end

# ╔═╡ f712027f-7f6e-4082-95f3-cae200d7f18a
function pseudo_voigt(x, fwhm)
	return n * Gaussian(x, fwhm) + (1-n) * Lorentzian(x, fwhm)
end

# ╔═╡ 298367a4-6f27-42dd-a115-e3235371347c


# ╔═╡ Cell order:
# ╠═2dada36d-12c5-41cc-a941-5da6f45a8e39
# ╠═ab93661e-92ba-11ed-0db0-a5b6da440669
# ╠═c1159303-2800-460e-a6ff-e957681e4c81
# ╠═f712027f-7f6e-4082-95f3-cae200d7f18a
# ╠═298367a4-6f27-42dd-a115-e3235371347c
