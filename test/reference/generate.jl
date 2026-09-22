# Regenerate the reference patterns. Run only when a change of the numbers is
# intended, from the project root:
#     julia --project=. test/reference/generate.jl

using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "..", "functions.jl"))
include(joinpath(@__DIR__, "reference.jl"))

for mode in REFERENCE_MODES
    x, columns = reference_patterns(mode)
    df = DataFrame("x" => x, columns...)
    CSV.write(joinpath(REFERENCE_DIR, "$mode.csv"), df)
    println("Wrote $mode.csv: $(nrow(df)) points, $(length(columns)) samples")
end
