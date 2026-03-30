"""
simple XRD
by Hezy Amiel
2023 - 2024
"""

using Plots
using Random
using DataFrames
using CSV


"""
================
main
================
"""

include("functions.jl")

Random.seed!(347) # Setting the seed for random noise

isdir("results") || mkdir("results")

_, _, lattice = read_xrd_config("data.toml")
lattice_types = sort(collect(keys(lattice)))

θ₀ = do_it_zero("data.toml")
df = DataFrame(Dict("θ" => θ₀, lt => θ₀ for lt in lattice_types))

for lattice_type in lattice_types
    local twoθ, intensities, the_title, the_plot = do_it("data.toml", lattice_type, :dark)
    df[:, "θ"] = twoθ
    df[:, lattice_type] = intensities

    display(the_plot)
    println("$lattice_type. Press Enter to continue...")
    readline()  # Wait for user input to continue
    savefig(the_plot, "./results/$the_title")
end

closeall()  # Close all plots before finishing
CSV.write("./results/XRD_results.csv", df)
println("Done. Results saved to ./results/")
