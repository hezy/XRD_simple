"""
XRD sim
by Hezy Amiel
2023--2026
"""

using Plots
using Random
using DataFrames
using CSV
using ArgParse

include("functions.jl")


# VS Code's Julia extension loads VSCodeServer into Main and routes plots to a
# persistent plot pane. In that context we skip the between-plots pause and
# don't close plot windows on exit.
const IN_VSCODE = isdefined(Main, :VSCodeServer)


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
            help = "Skip interactive pauses between plots (automatic in VS Code)"
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
    interactive = !args["no-interactive"] && !IN_VSCODE
    save_plots = !args["no-plots"]

    Random.seed!(seed)

    isdir("results") || mkdir("results")

    _, _, samples = read_xrd_config(config_file)

    θ₀ = do_it_zero(config_file)
    df = DataFrame("θ" => θ₀)

    for (structure, element, a) in samples
        local twoθ, intensities, title, the_plot = do_it(config_file, structure, element, a, plot_theme)
        df[:, "θ"] = twoθ
        df[!, title] = intensities

        if interactive || IN_VSCODE
            display(the_plot)
        end
        if interactive
            println("$title. Press Enter to continue...")
            readline()
        end

        if save_plots
            savefig(the_plot, "./results/$title")
        end
    end

    IN_VSCODE || closeall()
    CSV.write("./results/XRD_results.csv", df)
    println("Produced $(length(samples)) samples. Results saved to ./results/")
end


main()
