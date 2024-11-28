using Plots

include("voigt_new_with_docstrings.jl")

θ = collect(range(2.0, 160.0, 1000))

y1 = Voigt_peak(θ, 25.0, 10.0, 4.0, 4.0)
y2 = pseudo_Voigt_peak(θ, 75.0, 10.0, 4.0, 4.0)
# y = y1 + y2

plot(θ, y1)
plot!(θ, y2)
# plot!(θ, y)
