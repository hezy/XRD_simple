using Test

@testset "background" begin
    xray = XRay(1.5418, deg2rad(10.0), deg2rad(120.0), 1e-4, 5e-5, 1e-5)
    two_θ = collect(LinRange(deg2rad(10.0), deg2rad(120.0), 1000))

    bg = background(xray, two_θ)
    @test all(bg .>= 0)
    @test bg[1] > bg[end]
    @test bg == background(xray, two_θ)

    electron = Electron(200.0, 0.0, 1.2, 0.005, 50.0, 800, 2.5, true, 0.5, 0.0)
    g = collect(LinRange(0.0, 1.2, 1000))

    bg = background(electron, g)
    @test all(bg .>= 0)
    @test bg[1] > bg[end]
end
