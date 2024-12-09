using SolarDynamo
using Test

@testset "SolarDynamo.jl" begin


    θ = [2.0,  # τ
         3,     # T
         8.5,   # Nd
         0.1,   # sigma
         10]    # Bmax

    r1 = sn(θ, Tobs=100,  seed=34)
    @test length(r1) == 100

    r2 = sn(θ, Tobs=100)
    r3 = sn(θ, Tobs=100, seed=34)

    @test !(r1 ≈ r2)
    @test r1 ≈ r3

end
