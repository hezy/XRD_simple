using Test

include(joinpath(@__DIR__, "..", "functions.jl"))

@testset "read_xrd_config" begin
    instrument, peak_width, lattice = read_xrd_config("data.toml")

    @test haskey(instrument, "two_theta_min")
    @test haskey(instrument, "two_theta_max")
    @test haskey(instrument, "N")
    @test haskey(instrument, "lambda")
    @test haskey(peak_width, "U")
    @test haskey(peak_width, "V")
    @test haskey(peak_width, "W")
    @test haskey(peak_width, "K")
    @test haskey(peak_width, "Epsilon")
    @test haskey(peak_width, "D")

    @test instrument["two_theta_min"] < 1.0
    @test instrument["two_theta_max"] < 3.0

    @test haskey(lattice, "SC") || !haskey(lattice, "SC")
    @test haskey(lattice, "BCC") || !haskey(lattice, "BCC")
    @test haskey(lattice, "FCC") || !haskey(lattice, "FCC")
end
