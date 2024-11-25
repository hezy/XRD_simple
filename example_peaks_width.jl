using Plots

include("functions.jl")

θ_min = 1.0
θ_max = 80.0

θ_deg = collect(range(θ_min, θ_max, 1000))
θ_rad = collect(range(deg2rad(θ_min), deg2rad(θ_max), 1000))

w_L = Lorentzian_peaks_width(θ_rad, 0.9, 0.02, 1.5418, 500.0)
w_G = Gaussian_peaks_width(θ_rad, 0.001, 0.005, 0.01)
w_eff = peak_fwhm(w_L, w_G)

plot(θ_deg, w_L)
plot!(θ_deg, w_G)
plot!(θ_deg, w_eff)

