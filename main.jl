"""
simple XRD
by Hezy Amiel
2023 - 2024
"""

#using Plots; gr()
#using SpecialFunctions
using Random
#using Distributions
using DataFrames
using CSV


"""
================
main
================
"""

include("functions.jl")

Random.seed!(347) # Setting the seed for random noise

θ₀ = do_it_zero("data.toml")
df = DataFrame(θ=θ₀, SC=θ₀, BCC=θ₀, FCC=θ₀)

for lattice_type in ("SC", "BCC", "FCC")
    df[:, "θ"], df[:, lattice_type], the_title, the_plot = do_it("data.toml", lattice_type, :dark)

    display(the_plot) 
    sleep(1)
    println("$lattice_type. Press Enter to continue...")
    readline()  # Wait for user input to continue    
    savefig(the_plot, "./results/$the_title")
end

CSV.write("./results/XRD_results.csv", df)

