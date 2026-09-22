# XRD sim
# by Hezy Amiel
# 2023--2026
# Julia > 1.8
#
# Entry file of the simulation (physics only). Include it, then `plotting.jl` if
# plots are wanted.


using SpecialFunctions
using Distributions: Normal
using TOML


include("config.jl")
include("crystal.jl")
include("profiles.jl")
include("xray.jl")
include("electron.jl")
include("simulate.jl")
