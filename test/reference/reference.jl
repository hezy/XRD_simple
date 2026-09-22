# Shared by generate.jl and test_reference.jl: compute every sample of one
# reference configuration, with a fixed seed. Returns (x, columns), where
# columns is a vector of (title, intensities) pairs in sample order.

using Random

const REFERENCE_DIR = @__DIR__
const REFERENCE_MODES = ("xray", "electron")
const REFERENCE_SEED = 347

function reference_patterns(mode::String)
    config = joinpath(REFERENCE_DIR, "$mode.toml")
    _, _, samples = read_xrd_config(config)
    Random.seed!(REFERENCE_SEED)

    x = Float64[]
    columns = Pair{String,Vector{Float64}}[]
    for (structure, element, a) in samples
        x, y, title, _ = do_it(config, structure, element, a, :default)
        push!(columns, title => y)
    end
    return x, columns
end
