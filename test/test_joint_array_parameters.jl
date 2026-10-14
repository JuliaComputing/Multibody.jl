using ModelingToolkit, Multibody, Test

@testset "Joint array parameters" begin
    gear = GearConstraint(; name = :gear, ratio = 10)
    @test gear isa System
    universal = Universal(; name = :universal)
    @test universal isa System
end
