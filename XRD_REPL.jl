"""
simple_XRD_new.jl
by Hezy Amiel
April 2023
Julia 1.8.5
"""

using Plots
gr()
#using SpecialFunctions
#using Random
#using Distributions
using DataFrames
#using CSV


"""
================
main
================
"""

include("XRD_func.jl")

Random.seed!(347) # Setting the seed for random noise

θ₀ = do_it_zero("simple_XRD.txt")
df = DataFrame(θ=θ₀, SC=θ₀, BCC=θ₀, FCC=θ₀)

for lattice_type in ("SC", "BCC", "FCC")
    df[:, "θ"], df[:, lattice_type], the_title, plot_name = do_it("simple_XRD.txt", lattice_type)

    display(plot_name) 
    sleep(1)
    println("$plot_name . Press Enter to continue...")
    readline()  # Wait for user input to continue    
    savefig(plot_name, the_title)

end

CSV.write("XRD_results.csv", df)

