using Test

include(joinpath(@__DIR__, "..", "functions.jl"))

@testset "XRD_simple" begin
    include("test_crystal.jl")
    include("test_peak_profiles.jl")
    include("test_widths.jl")
    include("test_pattern.jl")
    include("test_background.jl")
    include("test_config.jl")
    include("test_errors.jl")
end
