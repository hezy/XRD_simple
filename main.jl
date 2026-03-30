"""
simple XRD
by Hezy Amiel
2023 - 2024
"""

using Plots
using Random
using DataFrames
using CSV
using ArgParse

include("functions.jl")


"""
    setup_argparse() -> ArgParseSettings

Create argument parser for XRD simulation.
"""
function setup_argparse()
    s = ArgParseSettings(
        description = "Powder X-ray diffraction simulation for cubic crystal structures.",
    )

    @add_arg_table s begin
        "--config"
            help = "Path to TOML configuration file"
            default = "data.toml"
        "--theme"
            help = "Plot theme (e.g., dark, light, ggplot2)"
            default = "dark"
        "--seed"
            help = "Random seed for reproducibility"
            arg_type = Int
            default = 347
        "--no-interactive"
            help = "Skip interactive pauses between plots"
            action = :store_true
        "--no-plots"
            help = "Skip saving plots (CSV output only)"
            action = :store_true
    end

    return s
end


function main()
    args = ArgParse.parse_args(ARGS, setup_argparse())

    config_file = args["config"]
    plot_theme = Symbol(args["theme"])
    seed = args["seed"]
    interactive = !args["no-interactive"]
    save_plots = !args["no-plots"]

    Random.seed!(seed)

    isdir("results") || mkdir("results")

    _, _, lattice = read_xrd_config(config_file)
    lattice_types = sort(collect(keys(lattice)))

    if isempty(lattice_types)
        error("No active lattice types found in $config_file")
    end

    θ₀ = do_it_zero(config_file)
    df = DataFrame(Dict(lt => θ₀ for lt in lattice_types))
    df[!, "θ"] = θ₀

    for lattice_type in lattice_types
        local twoθ, intensities, title, the_plot = do_it(config_file, lattice_type, plot_theme)
        df[:, "θ"] = twoθ
        df[:, lattice_type] = intensities

        if interactive
            display(the_plot)
            println("$lattice_type. Press Enter to continue...")
            readline()
        end

        if save_plots
            savefig(the_plot, "./results/$title")
        end
    end

    closeall()
    CSV.write("./results/XRD_results.csv", df)
    println("Done. Results saved to ./results/")
end


main()
