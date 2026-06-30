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


"""
    write_ring_outputs(instrument, structure, element, a, g, intensities, title)

Render the electron-mode g-profile as a 2D ring image and write the matching
reflection answer key. Outputs to `results/rings/`:
- `{title}.png`              — the Debye–Scherrer ring image (student-facing)
- `{title}_reflections.csv`  — hidden key: hkl, N, g, ring radius (mm), multiplicity
"""
function write_ring_outputs(instrument, structure, element, a, g, intensities, title)
    isdir("results/rings") || mkpath("results/rings")

    camera_constant = Float64(get(instrument, "camera_constant", 50.0))

    ring_plot = render_ring_image(g, intensities, camera_constant;
        image_px    = Int(get(instrument, "image_px", 800)),
        beam_stop_mm = Float64(get(instrument, "beam_stop_mm", 2.5)),
        phosphor    = Bool(get(instrument, "ring_phosphor", true)),
        gamma       = Float64(get(instrument, "ring_gamma", 0.5)),
        noise_level = Float64(get(instrument, "ring_noise", 0.0)))
    savefig(ring_plot, "./results/rings/$title.png")

    rt = reflection_table(structure, a, instrument["g_max"])
    key = DataFrame(
        h = [hkl[1] for hkl in rt.indices],
        k = [hkl[2] for hkl in rt.indices],
        l = [hkl[3] for hkl in rt.indices],
        N = rt.N,
        g_per_A = rt.g,
        r_mm = camera_constant .* rt.g,
        multiplicity = rt.multiplicity,
    )
    CSV.write("./results/rings/$(title)_reflections.csv", key)
    return ring_plot
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

    instrument, _, samples = read_xrd_config(config_file)

    # Electron diffraction is plotted vs scattering vector g (1/Å); X-ray vs 2θ.
    is_electron = get(instrument, "radiation", "xray") == "electron"
    xcol = is_electron ? "g (1/Å)" : "2θ (deg)"

    x₀ = do_it_zero(config_file)
    df = DataFrame(xcol => x₀)

    for (structure, element, a) in samples
        local x, intensities, title, the_plot = do_it(config_file, structure, element, a, plot_theme)
        df[:, xcol] = x
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
            if is_electron
                write_ring_outputs(instrument, structure, element, a, x, intensities, title)
            end
        end
    end

    IN_VSCODE || closeall()
    CSV.write("./results/XRD_results.csv", df)
    println("Produced $(length(samples)) samples. Results saved to ./results/")
end


main()
