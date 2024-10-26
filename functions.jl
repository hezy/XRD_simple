"""
simple_XRD_new.jl
by Hezy Amiel
April 2023
Julia 1.8.5
"""

using Plots; gr()
using SpecialFunctions
#using Random
using Distributions
#using DataFrames
#using CSV


""" 
=========
Functions
=========
"""


"Returns a Voigt peak centered around θ₀, with amplitude A, width w, and mixing factor n "
function Voigt_peak(θ::Vector{Float64},
                    θ₀::Float64,
                    A::Float64,
                    w_L::Vector{Float64},
                    w_G::Vector{Float64},
                    n::Float64
                    )::Vector{Float64}
                    
    γ = w_L / 2
    σ = w_G / (2√(2log(2)))
    z = @. -im * ((θ - θ₀) + im * γ) / (√2 * σ)
    return @. A * real(erfcx(z)) / (√(2π) * σ)
end


"Returns a pseudo Voigt peak centered around θ₀, with amplitude A, width w, and mixing factor n "
function pseudo_Voigt_peak(θ::Vector{Float64},
                           θ₀::Float64,
                           A::Float64,
                           w::Vector{Float64},
                           n::Float64
                           )::Vector{Float64}
                           
    γ = w / 2
    σ = w / (2√(2log(2)))
    return @. A * (n * pdf.(Cauchy(θ₀, γ), θ) + (1 - n) * pdf.(Normal(θ₀, σ), θ))
    # equivalent to:
    # return @. A * (n* (γ / π) / ((θ - θ₀)^2 + γ^2) + (1 - n) * 1 / √(2π) / σ * exp(-(θ - θ₀)^2 / 2σ^2)
end


"Returns the width of a peak as a function of 2θ with U, V, W parameters"
function peaks_width(two_θ_deg::Vector{Float64},
                     U::Float64,
                     V::Float64,
                     W::Float64
                     )::Vector{Float64}

        return @. √(U * tan(two_θ_deg * π / 360)^2 + V * tan(two_θ_deg * π / 360) + W)
end


"Calculating the Bragg angles coresponding to each d-spacing"
function bragg_angels(wavelength::Float64,
                      d_spacings::Vector{Float64}
                      )::Vector{Float64}

    sinθ = wavelength ./ (2 * d_spacings)
    sinθ_cleaned = [item for item in sinθ if abs(item) <= 1]  # removing values outside (-1,1)
    return 2 * (180 / π) * asin.(sinθ_cleaned)  # *2 for 2θ  
end



"Returnes the inter-layers distances as a function of Miller_indices"
function d_list(indices::Vector{Vector{Int8}},
                a::Float64
                )::Vector{Float64}

    return a ./ [√(i^2 + j^2 + k^2) for (i, j, k) in indices]
end


"sums peak functions to return intensity vs angle"
function sum_peaks(θ::Vector{Float64},
                   two_θ_list::Vector{Float64},
                   U::Float64,
                   V::Float64,
                   W::Float64
                   )::Vector{Float64}
    
    y = zeros(size(θ))
    for item in two_θ_list
        y = y + pseudo_Voigt_peak(θ, item, 1.0, peaks_width(θ, U, V, W), 0.5)
    end
    return y
end


"Building the XRD patterns"
function intensity_vs_angle(θ::Vector{Float64},
                            indices::Vector{Vector{Int8}},
                            λ::Float64,
                            a::Float64,
                            U::Float64,
                            V::Float64,
                            W::Float64
                            )::Vector{Float64}
    
    two_θ_list = bragg_angels(λ, d_list(indices, a))
    y = sum_peaks(θ, two_θ_list, U, V, W)
    return y
end


"Returns a list of Miller indices for each one of the cubic symmetries"
function Miller_indices(cell_type::String,
                        min::Int8,
                        max::Int8
                        )::Vector{Vector{Int8}}
    
    if !(cell_type in ["SC", "BCC", "FCC"])
        error("Invalid cell_type: $cell_type. Expected 'SC', 'BCC', or 'FCC'.")
    end
    if min > max
        error("Minimum value cannot be greater than maximum value.")
    end
    if !(isa(min, Int8) && isa(max, Int8))
        error("Minimum and maximum values must be integers.")
    end

    if cell_type == "SC"
        # In simple cubic lattice, all Miller indices are allowed
        return [
            [h, k, l] for h = min:max for k = min:max for l = min:max if [h, k, l] != [0, 0, 0]
        ]

    elseif cell_type == "BCC"
        # In body centered cubic lattice, only indices with h+k+l=even are allowed
        return [
            [h, k, l] for h = min:max for k = min:max for l = min:max if iseven(h + k + l) && [h, k, l] != [0, 0, 0]
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


"background function for the XRD pattern"
function background(θ::Vector{Float64}
                    )::Vector{Float64}
    
    return @. 2 + θ * (360 - θ) / 15000
end


"Adding some noise to the data"
function make_noisy(θ::Vector{Float64},
                    y::Vector{Float64}
                    )::Vector{Float64}
    
    return (background(θ) + y) .* rand(Normal(1, 0.1), size(θ))
end


"Reading a text file with instrument data, and lattice parameters"
function read_file(filename::String
                   )::Tuple{Dict,Dict}
    
    instrument_data = Dict{AbstractString,Any}()
    lattice_params = Dict{AbstractString,Float64}()

    # read file line by line
    for line in eachline(filename)
        # split the line by whitespace and remove empty strings
        tokens = filter(x -> x ≠ "", split(line))

        if length(tokens) > 0 && tokens[1] ≠ "#"
            if tokens[1] in ["θ_min", "θ_max"]
                instrument_data[tokens[1]] = parse(Float64, tokens[2])
            elseif tokens[1] == "N"
                instrument_data[tokens[1]] = parse(Int64, tokens[2])    
            elseif tokens[1] == "λ"
                instrument_data[tokens[1]] = parse(Float64, tokens[2])
            elseif tokens[1] in ["U", "V", "W"]
                instrument_data[tokens[1]] = parse(Float64, tokens[2])
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

    return instrument_data, lattice_params
end


"colecting input data, building the XRD pattern with background and noise, plotting it"
function do_it_zero(file_name::String
                    )::Vector{Float64}
    
    instrument_data, lattice_params = read_file(file_name)
    θ = collect(LinRange(instrument_data["θ_min"], instrument_data["θ_max"], instrument_data["N"]))
    return θ
end


"colecting input data, building the XRD pattern with background and noise, plotting it"
function do_it(file_name::String,
               lattice_type::String,
               plot_theme::Symbol
               )::Tuple{Vector,Vector,String,Plots.Plot}
    
    instrument_data, lattice_params = read_file(file_name)

    N = instrument_data["N"]
    θ = collect(LinRange(instrument_data["θ_min"], instrument_data["θ_max"], instrument_data["N"]))
    y = zeros(instrument_data["N"])
    λ = instrument_data["λ"]
    U, V, W = instrument_data["U"], instrument_data["V"], instrument_data["W"]

    a = lattice_params[lattice_type]

    index_min::Int8 = -5
    index_max::Int8 = 5 
    indices = Miller_indices(lattice_type, index_min, index_max)

    y = (background(θ) +
         intensity_vs_angle(θ, indices, λ, a, U, V, W)) .*
        rand(Normal(1, 0.1), N)

    the_title = "XRD - " * lattice_type

    theme(plot_theme)
    
    the_plot = plot(θ, y, title=the_title, xlabel="2θ (deg)", ylabel="Intensity (arb.)")

    return θ, y, the_title, the_plot
end

