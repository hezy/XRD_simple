# XRD sim: plotting
# by Hezy Amiel
# 2023--2026
#
# All Plots.jl calls of the simulation. Include after `XRDSim.jl`.


using Plots; gr()


# Phosphor-green colour ramp (black → dark green → bright green → highlight),
# the look of a fluorescent ED viewing screen.
const PHOSPHOR_RAMP = ["#000000", "#022b06", "#1f9b3a", "#5dff7a", "#e6ffe9"]


"""
    plot_title(mode, title) -> String

Plot title of one pattern. The electron title also states the wavelength.
"""
plot_title(::XRay, title::String) = title
plot_title(m::Electron, title::String) =
    "$title  (e⁻, λ=$(round(electron_wavelength(m.voltage_kV * 1000.0), digits=4)) Å)"


"""
    plot_pattern(mode, x, y, title, plot_theme) -> Plots.Plot

Plot one pattern from `simulate`, with the axis label and title of `mode`.
"""
function plot_pattern(mode::Radiation,
                      x::AbstractVector{<:Real},
                      y::AbstractVector{<:Real},
                      title::String,
                      plot_theme::Symbol
                      )::Plots.Plot
    theme(plot_theme)
    return plot(x, y, title=plot_title(mode, title), xlabel=axis_label(mode),
                ylabel="Intensity (arb.)", show=false)
end


"""
    plot_ring_image(coords, img, mode::Electron) -> Plots.Plot

Draw a ring image from `ring_image` as a square heatmap (no axes or frame),
ready to `savefig`. The colormap is phosphor green if `mode.ring_phosphor`,
else grayscale.
"""
function plot_ring_image(coords::AbstractVector{<:Real},
                         img::AbstractMatrix{<:Real},
                         mode::Electron
                         )::Plots.Plot
    cmap = mode.ring_phosphor ? cgrad(PHOSPHOR_RAMP) : cgrad(:grays)
    return heatmap(coords, coords, img;
                   c = cmap, aspect_ratio = :equal, colorbar = false,
                   axis = false, ticks = false, framestyle = :none,
                   legend = false, grid = false, widen = false,
                   background_color = :black, margin = 0 * Plots.mm,
                   size = (mode.image_px, mode.image_px), clims = (0, 1),
                   show = false)
end
