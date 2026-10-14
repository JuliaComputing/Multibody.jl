using ModelingToolkit, Multibody, OrdinaryDiffEq, Test
@testset "Joint orientation flags" begin
    for (wfixed, zfixed) in ((false, false), (true, false), (false, true), (true, true))
        @named world = World(g = 0)
        @named joint = FreeMotion(state = false, isroot = false,
            w_rel_a_fixed = wfixed, z_rel_a_fixed = zfixed)
        @named body = Body(state = true, isroot = true)
        @named sys = System(
            [connect(world.frame_b, joint.frame_a),
                connect(joint.frame_b, body.frame_a)],
            t; systems = [world, joint, body])
        compiled = multibody(sys)
        sol = solve(ODEProblem(compiled, [], (0.0, 0.1)), Tsit5())
        @test SciMLBase.successful_retcode(sol)
        @test sol(0.1; idxs = joint.w_rel_b) == zeros(3)
    end
    for damping in (0.0, 0.2)
        @named world = World(g = 0)
        @named joint = Spherical(state = false, d = damping)
        @named body = Body(state = true, isroot = true)
        @named sys = System(
            [connect(world.frame_b, joint.frame_a),
                connect(joint.frame_b, body.frame_a)],
            t; systems = [world, joint, body])
        compiled = multibody(sys)
        sol = solve(ODEProblem(compiled, [], (0.0, 0.1)), Tsit5())
        @test SciMLBase.successful_retcode(sol)
        @test sol(0.1; idxs = joint.w_rel) == zeros(3)
    end
end

@testset "Rotated spherical damping" begin
    for (state, damping) in Iterators.product((false, true), (0.0, 0.002))
        inertia = 0.001
        quaternion = [cos(pi / 4), sin(pi / 4), 0, 0]
        @named world = World(g = 0)
        @named joint = Spherical(; state, quat = state, d = damping)
        @named body = Body(isroot = !state, quat = !state,
            I_11 = inertia, I_22 = inertia, I_33 = inertia,
            Q = state ? nothing : quaternion,
            w_a = state ? nothing : [0, 0, 1])
        @named sys = System(
            [connect(world.frame_b, joint.frame_a),
                connect(joint.frame_b, body.frame_a)],
            t; systems = [world, joint, body])
        compiled = multibody(sys)
        initial = state ?
                  [collect(joint.Q̂) .=> quaternion; collect(joint.w_rel) .=> [0, 0, 1]] :
                  []
        sol = solve(ODEProblem(compiled, initial, (0.0, 0.1)),
            Tsit5(); reltol = 1e-8, abstol = 1e-10)
        @test SciMLBase.successful_retcode(sol)
        omega = sol(0.1; idxs = body.w_a)
        @test sol(0.1; idxs = joint.w_rel)≈omega atol=1e-7
        @test omega≈[0, 0, exp(-damping / inertia * 0.1)] atol=1e-7
        R = reshape(sol(0.1; idxs = vec(ori(body.frame_a).R)), 3, 3)
        @test sol(0.1; idxs = joint.frame_b.tau)≈damping * omega atol=1e-7
        @test sol(0.1; idxs = joint.frame_a.tau)≈-R' * (damping * omega) atol=1e-7
    end
end

@testset "Spherical-spherical scalar constraint" begin
    @named world = World(g = 0)
    @named rod = SphericalSpherical(r_0 = [1, 0, 0], m = 1, kinematic_constraint = false)
    @named body = Body(isroot = true, r_0 = [1, 0, 0])
    @named sys = System(
        [connect(world.frame_b, rod.frame_a),
            connect(rod.frame_b, body.frame_a)],
        t; systems = [world, rod, body])
    compiled = multibody(sys)
    @test compiled isa System
end
