using ModelingToolkit, Multibody, Test

@testset "BodyShape initialization" begin
    @named world = World()
    @named joint = Prismatic(n = [1, 0, 0])
    @named cart = BodyShape(m = 50.0, r = [0.4, 0, 0])
    eqs = [
        connect(world.frame_b, joint.frame_a),
        connect(joint.frame_b, cart.frame_a)
    ]
    @named model = System(eqs, Multibody.t; systems = [world, joint, cart])
    ssys = multibody(model)
    prob = ODEProblem(ssys, [joint.s => 0.0, joint.v => 0.0], (0.0, 1.0))
    @test length(prob.u0) == 2
end
