"""
simple_XRD.jl
by Hezy Amiel
January 2023
Julia 1.8.5
"""

using Plots, SpecialFunctions, Random, Distributions, DataFrames, CSV


""" 
=========
Functions
=========
"""

function Lorentzian(x, fwhm)
    γ = fwhm / 2
    return pdf.(Cauchy(0.0, γ), x)
    # equivalent to:
    # return @. (γ / pi) / (x^2 + γ^2)
end


function Gaussian(x, fwhm)
    σ = fwhm / (2√(2log(2)))
    return pdf.(Normal(0.0, σ), x)
    # equivalent to:
    #return @. 1 / √(2π) / σ * exp(-x^2 / 2σ^2)
end


mix_fun(f1, f2, n) = n * f1 + (1 - n) * f2


function Pseudo_Voigt(x, fwhm, n)
    return n * Lorentzian(x, fwhm) + (1 - n) * Gaussian(x, fwhm)
end


pseudo_Voigt(x, fwhm, n) =  mix_fun(Lorentzian(x, fwhm), Gaussian(x, fwhm), n)


function Voigt(x, fwhm_L, fwhm_G)
    γ = fwhm_L / 2
    σ = fwhm_G / (2√(2log(2)))
    z = @. -im * (x + im * γ) / (√2 * σ)
    return @. real(erfcx(z)) / (√(2pi) * σ)
end


function peak(θ, θ₀, A, w, n)
    return @. A * pseudo_Voigt(θ-θ₀, w, n)
    # return @. A * Voigt(θ-θ₀, w, n)
end


function peaks_width(two_θ_deg, U, V, W)
    two_θ_rad = two_θ_deg * π / 180
    return @. √(U * tan(two_θ_rad / 2)^2 + V * tan(two_θ_rad / 2) + W)
end


function bragg_angels(wavelength, d_spacings)
    sinθ = wavelength ./ (2 * d_spacings)
    sinθ_cleaned = [item for item in sinθ if abs(item) <= 1]  # removing values outside (-1,1)
    return 2 * (180 / π) * asin.(sinθ_cleaned)  # *2 for 2θ  
end


function d_list(indices, a)
    return a ./ .√(sum(indices .^ 2, dims = 2))
end


function sum_peaks(θ, two_θ_list, U, V, W)
    y = zeros(size(θ))
    for item in two_θ_list
        y = y + peak(θ, item, 1, peaks_width(θ, U, V, W), 0.5)
    end
    return y
end


 function intensity_vs_angle(θ, indices, λ, a, U, V, W)
     indices = (reduce(hcat, indices))'
     two_θ_list = bragg_angels(λ, d_list(indices, a))
     y = sum_peaks(θ, two_θ_list, U, V, W)
     return y
end



function plot_it(θ, y, Title)
     default(show = true)
     p = plot(θ, y)
     title!(Title)
     xlabel!(raw"2θ (deg)")
     ylabel!(raw"Intensity (arb.)")
     return p
 end


# function plot_it(θ, y; title="", xlabel="2θ (deg)", ylabel="Intensity (arb.)", show_plot=true)
#     p = plot(θ, y, xlabel=xlabel, ylabel=ylabel, title=title)
#     if show_plot
#         default(show=true)
#         display(p)
#     end
#     return p
# end


function Miller_indices(cell_type, min, max)
    if !(cell_type in ["SC", "BCC", "FCC"])
        error("Invalid cell_type: $cell_type. Expected 'SC', 'BCC', or 'FCC'.")
    end
    if min > max
        error("Minimum value cannot be greater than maximum value.")
    end
    if !(isa(min, Int) && isa(max, Int))
        error("Minimum and maximum values must be integers.")
    end
    
    if cell_type == "SC"
        # In simple cubic lattice, all Miller indices are allowed
        return [
            [h, k, l] for h = min:max for k = min:max for
            l = min:max if [h, k, l] != [0, 0, 0]
        ]
    elseif cell_type == "BCC"
        # In body centered cubic lattice, only indices with h+k+l=even are allowed
        return [
            [h, k, l] for h = min:max for k = min:max for
            l = min:max if iseven(h + k + l) && [h, k, l] != [0, 0, 0]
        ]
    elseif cell_type == "FCC"
        # In face centered cubic lattice, h,k,l must all be either odd or even
        return [
            [h, k, l] for h = min:max for k = min:max for l = min:max if
            ((iseven(h) && iseven(k) && iseven(l)) || (isodd(h) && isodd(k) && isodd(l))) &&
            [h, k, l] != [0, 0, 0]
        ]
    end
end



function background(θ)
    return @. 2 + θ * (360 - θ) / 15000
end


function make_noisy(θ, y)
    return (background(θ) + y) .* rand(Normal(1, 0.1), size(θ))
end


function read_file(filename::AbstractString)
    input_data = Dict{AbstractString,Any}()
    lattice_params = Dict{AbstractString,Float64}()

    # read file line by line
    for line in eachline(filename)
        # split the line by whitespace and remove empty strings
        tokens = filter(x -> x ≠ "", split(line))

        if length(tokens) > 0 && tokens[1] ≠ "#"
            if tokens[1] in ["θ_min", "θ_max", "N"]
                input_data[tokens[1]] = parse(Int, tokens[2])
            elseif tokens[1] == "λ"
                input_data[tokens[1]] = parse(Float64, tokens[2])
            elseif tokens[1] in ["U", "V", "W"]
                input_data[tokens[1]] = parse(Float64, tokens[2])
            elseif tokens[1] in ["SC", "BCC", "FCC"]
                lattice_params[tokens[1]] = parse(Float64, tokens[3])
            end
        elseif length(tokens) > 1 && tokens[1] == "BCC" && tokens[2] ≠ "#"
            lattice_params[tokens[1]] = parse(Float64, tokens[3])
        elseif length(tokens) > 1 && tokens[1] == "FCC" && tokens[2] ≠ "#"
            lattice_params[tokens[1]] = parse(Float64, tokens[3])
        elseif length(tokens) > 1 && tokens[1] == "SC" && tokens[2] ≠ "#"
            lattice_params[tokens[1]] = parse(Float64, tokens[3])
        end
    end

    return input_data, lattice_params
end


"""
================
General Settings
================
"""

inst, lattice_par = read_file("simple_XRD.txt")

N = inst["N"]
θ = collect(LinRange(inst["θ_min"], inst["θ_max"], inst["N"]))
y = zeros(inst["N"])
λ = inst["λ"]
U, V, W = inst["U"], inst["V"], inst["W"] 

Random.seed!(347) # Setting the seed for random noise


"""
============
Simple Cubic
============
"""

"""
Lattice parameter for SC Polonium (α-Po)
from https://en.wikipedia.org/wiki/Polonium 
"""
a_SC = 0.3352

indices_SC = Miller_indices("SC", -5, 5)

y_SC =
     (background(θ) + intensity_vs_angle(θ, indices_SC, λ, a_SC, U, V, W)) .*
     rand(Normal(1, 0.1), N)

 plot1 = plot_it(θ, y_SC, "XRD - SC")
display(plot1)
savefig(plot1, "SC")


"""
===================
Body centered Cubic
===================
"""

"""
Lattice parameter for BCC Tantalum (α-Ta)
from https://en.wikipedia.org/wiki/Tantalum
"""
a_BCC = 0.33058

indices_BCC = Miller_indices("BCC", -5, 5)

y_BCC =
    (background(θ) + intensity_vs_angle(θ, indices_BCC, λ, a_BCC, U, V, W)) .*
    rand(Normal(1, 0.1), N)

plot2 = plot_it(θ, y_BCC, "XRD - BCC")
display(plot2)
savefig(plot2, "BCC")


"""
===================
Face Centered Cubic
===================
"""

"""
Lattice parameter for FCC Platinum
from https://periodictable.com/Elements/078/data.html
"""
a_FCC = 0.39242

indices_FCC = Miller_indices("FCC", -5, 5)

y_FCC =
    (background(θ) + intensity_vs_angle(θ, indices_FCC, λ, a_FCC, U, V, W)) .*
    rand(Normal(1, 0.1), N)

plot3 = plot_it(θ, y_FCC, "XRD - FCC")
display(plot3)
savefig(plot3, "FCC")

