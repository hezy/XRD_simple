using Test
using CSV
using DataFrames

include(joinpath(@__DIR__, "reference", "reference.jl"))

# The refactor (REFACTOR_PLAN.md, Phases 3–6) must not change the numbers.
@testset "reference output: $mode" for mode in REFERENCE_MODES
    saved = CSV.read(joinpath(REFERENCE_DIR, "$mode.csv"), DataFrame)
    x, columns = reference_patterns(mode)

    @test isapprox(x, saved.x; rtol=1e-12)
    @test [first(c) for c in columns] == names(saved)[2:end]
    for (title, y) in columns
        @test isapprox(y, saved[!, title]; rtol=1e-12)
    end
end
