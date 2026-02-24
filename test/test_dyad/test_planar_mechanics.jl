using ModelingToolkit, OrdinaryDiffEq, Test
using Multibody
using Multibody: multibody
using ModelingToolkit: t_nounits as t
import Symbolics
import Markdown

# Path to Dyad-generated test component definitions
const GENERATED_PM = joinpath(@__DIR__, "..", "..", "generated", "PlanarMechanics")

# Include generated test component definitions.
include(joinpath(GENERATED_PM, "PendulumTest_definition.jl"))
include(joinpath(GENERATED_PM, "BodyShapePendulumTest_definition.jl"))
include(joinpath(GENERATED_PM, "SpringDamperTest_definition.jl"))
include(joinpath(GENERATED_PM, "SpringDamperSystemTest_definition.jl"))
include(joinpath(GENERATED_PM, "SpringAndDamperTest_definition.jl"))

tspan3 = (0.0, 3.0)
tspan5 = (0.0, 5.0)

@testset "Dyad PlanarMechanics" begin

    # Translated from "Pendulum" in test_PlanarMechanics.jl
    @testset "PendulumTest" begin
        @named model = PendulumTest()
        ssys = multibody(model)
        @test length(unknowns(ssys)) == 2
        prob = ODEProblem(ssys, [], tspan3)
        sol = solve(prob, Rodas5P())
        @test SciMLBase.successful_retcode(sol)
        @test sol(1, idxs=ssys.rod.frame_a.phi) ≈ -2.881383661312169 atol=1e-2
        @test sol(2, idxs=ssys.rod.frame_a.phi) ≈ -1 atol=1e-2
    end

    # Translated from "Pendulum with body shape" in test_PlanarMechanics.jl
    @testset "BodyShapePendulumTest" begin
        @named model = BodyShapePendulumTest()
        ssys = multibody(model)
        @test length(unknowns(ssys)) == 2
        prob = ODEProblem(ssys, [], tspan3)
        sol = solve(prob, Rodas5P())
        @test SciMLBase.successful_retcode(sol)
        @test sol(1, idxs=ssys.rod.frame_a.phi) ≈ -pi atol=1e-2
        @test sol(2, idxs=ssys.rod.frame_a.phi) ≈ 0 atol=1e-2
    end

    # Custom Dyad test: pendulum with a spring
    @testset "SpringDamperTest" begin
        @named model = SpringDamperTest()
        ssys = multibody(model)
        prob = ODEProblem(ssys, [], tspan3)
        sol = solve(prob, Rodas5P())
        @test SciMLBase.successful_retcode(sol)
    end

    # Translated from "SpringDamper" in test_PlanarMechanics.jl
    @testset "SpringDamperSystemTest" begin
        @named model = SpringDamperSystemTest()
        ssys = multibody(model)
        prob = ODEProblem(ssys, [], tspan5)
        sol = solve(prob, Tsit5())
        @test SciMLBase.successful_retcode(sol)
    end

    # Translated from "Spring and damper demo" in test_PlanarMechanics.jl
    @testset "SpringAndDamperTest" begin
        @named model = SpringAndDamperTest()
        ssys = multibody(model)
        prob = ODEProblem(ssys, [], tspan5)
        sol = solve(prob, Tsit5())
        @test SciMLBase.successful_retcode(sol)
    end

end
