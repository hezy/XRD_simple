using Plots

include("functions.jl")

# Set up angle range and parameters
θ_deg = range(10, 80, length=1000)  # degrees
θ = deg2rad.(θ_deg)                 # convert to radians

# Parameters for broadening
λ = 0.154056   # Cu Kα wavelength in nm
D = 100.0      # crystallite size in nm
K = 0.9        # Scherrer constant
Y = K * λ / D  # size parameter

ε = 0.001      # 0.1% microstrain
X = 4 * ε      # strain parameter

# Calculate angle-dependent widths
w_L = @. X * tan(θ) + Y / cos(θ)    # total Lorentzian width

# Gaussian instrumental broadening (Caglioti formula)
U, V, W = 0.001, 0.0005, 0.00001
w_G = @. √(U * tan(θ)^2 + V * tan(θ) + W)

# Calculate peaks at different angles
angles = [15, 30, 45, 60]  # degrees
p = plot(xlabel="2θ (degrees)", 
        ylabel="Normalized Intensity",
        title="Peak Shapes at Different Angles")

# Calculate each peak over the full range
for θ₀ in deg2rad.(angles)
    peak = Voigt_peak(θ, θ₀, 1.0, w_L, w_G, normalize=true)
    plot!(θ_deg, peak, label="Peak at $(round(Int, rad2deg(θ₀)))°")
end

# Display plot
plot!(size=(800,600), legend=:topleft)
