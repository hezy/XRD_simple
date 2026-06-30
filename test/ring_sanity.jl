# Sanity check for the electron ring-output mode (Phase B0).
# Verifies, for each worked sample, that:
#   (1) the 1D g-profile's peaks coincide with the analytic g = √N/a positions
#       from reflection_table, and
#   (2) the radial image mapping r = camera_constant·g places each ring's
#       intensity maximum at the analytic radius (tests render_ring_image's map).
#
# Run:  julia --project=. test/ring_sanity.jl [config.toml]

include("../functions.jl")

cfg = length(ARGS) ≥ 1 ? ARGS[1] : "data.toml"
instrument, peak_width, samples = read_xrd_config(cfg)
@assert get(instrument, "radiation", "xray") == "electron" "config must be electron mode"

camera_constant = Float64(get(instrument, "camera_constant", 50.0))
g_max = instrument["g_max"]
tol_g = 2 * (g_max - get(instrument, "g_min", 0.0)) / instrument["N"]  # ~2 grid steps

# Find local maxima of v above a floor; returns the x-locations (parabolic-refined).
function peak_locations(x, v; rel_height = 0.02)
    thr = minimum(v) + rel_height * (maximum(v) - minimum(v))
    locs = Float64[]
    @inbounds for i in 2:length(v)-1
        if v[i] > thr && v[i] ≥ v[i-1] && v[i] > v[i+1]
            d = v[i-1] - 2v[i] + v[i+1]
            δ = d == 0 ? 0.0 : 0.5 * (v[i-1] - v[i+1]) / d   # parabolic vertex offset
            push!(locs, x[i] + δ * (x[i+1] - x[i]))
        end
    end
    return locs
end

# Nearest |a-b| match of each reference value into a candidate list.
nearest_err(refs, cands) = [minimum(abs.(cands .- r)) for r in refs]

println("Ring sanity check  (config: $cfg, camera_constant = $camera_constant mm·Å)")
println("="^70)

function run_checks(samples, instrument, peak_width, camera_constant, g_max, tol_g)
all_ok = true
for (structure, element, a) in samples
    g, y, _, _ = do_it_electron(instrument, peak_width, structure, element, a, :dark)
    rt = reflection_table(structure, a, g_max)

    # (1) profile peaks vs analytic g
    gpeaks = peak_locations(g, y)
    eg = nearest_err(rt.g, gpeaks)              # how far each analytic g is from a profile peak
    max_eg = maximum(eg)

    # (2) radial image: build the on-axis radial intensity and find its maxima
    npx = 4000
    rline = collect(LinRange(0.0, camera_constant * g_max, npx))
    vline = [radial_profile_value(y, first(g), last(g), r / camera_constant) for r in rline]
    rpeaks = peak_locations(rline, vline)
    r_analytic = camera_constant .* rt.g
    er = nearest_err(r_analytic, rpeaks)        # mm
    max_er = maximum(er)

    ok = max_eg ≤ tol_g && max_er ≤ 2 * camera_constant * tol_g
    all_ok &= ok
    println(rpad("$element-$structure", 10),
            " | reflections: ", rpad(length(rt.g), 3),
            " | max g-peak err: ", rpad(round(max_eg, sigdigits=2), 9), " 1/Å",
            " | max ring-radius err: ", rpad(round(max_er, sigdigits=2), 8), " mm",
            ok ? "  ✓" : "  ✗ FAIL")
end
return all_ok
end

all_ok = run_checks(samples, instrument, peak_width, camera_constant, g_max, tol_g)
println("="^70)
println(all_ok ? "ALL SANITY CHECKS PASSED" : "SANITY CHECK FAILED")
exit(all_ok ? 0 : 1)
