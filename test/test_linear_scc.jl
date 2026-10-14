using ModelingToolkit, Multibody, OrdinaryDiffEq, OrdinaryDiffEqRosenbrock,
      Test

@testset "Numerical linear algebraic loop at axis-aligned orientation" begin
    @named world = World()
    @named body = Body(m = 1, isroot = true, r_cm = zeros(3), quat = true)
    @named spring = Multibody.Spring(c = 1)
    @named model = System(
        [connect(world.frame_b, spring.frame_a),
            connect(spring.frame_b, body.frame_a)],
        t; systems = [world, spring, body])
    compiled = multibody(model;
        reassemble_alg = ModelingToolkit.StructuralTransformations.DefaultReassembleAlgorithm(
            inline_linear_sccs = true, analytical_linear_scc_limit = 1))
    @test length(unknowns(compiled)) == 13
    initial = [collect(body.r_0) .=> [0, -1e-5, 0];
               collect(body.v_0) .=> zeros(3); collect(body.w_a) .=> [0, 1, 0];
               collect(body.Q̂) .=> [1, 0, 0, 0]]
    prob = ODEProblem(compiled, initial, (0.0, 10.0))
    du = similar(prob.u0)
    prob.f(du, prob.u0, prob.p, 0.0)
    @test all(isfinite, du)
    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)
    @test sol(float(pi); idxs = body.r_0[2]) < -2
end
