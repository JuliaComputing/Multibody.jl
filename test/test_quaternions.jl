# test utils
import Multibody.Rotations.QuatRotation as Quat
import Multibody.Rotations
function get_R(sol, frame, t)
    reshape(sol(t, idxs=vec(ori(frame).R.mat)), 3, 3)
end
function get_r(sol, frame, t)
    sol(t, idxs=collect(frame.r_0))
end




# ==============================================================================
## Harmonic oscillator with Body as root and quaternions as state variables
# ==============================================================================
using LinearAlgebra
@testset "Harmonic oscillator with Body as root and quaternions as state variables" begin

@named body = Body(; m = 1, isroot = true, r_cm = [0.0, 0, 0], phi0 = [0, 0.9, 0], quat=true) # This time the body isroot since there is no joint containing state
@named spring = Multibody.Spring(c = 1)

connections = [connect(world.frame_b, spring.frame_a)
               connect(spring.frame_b, body.frame_a)]

@named model = ODESystem(connections, t, systems = [world, spring, body])
model = complete(model)
# ssys = structural_simplify(model, allow_parameter = false)

irsys = IRSystem(model)
ssys = structural_simplify(irsys)
@test length(unknowns(ssys)) == 13 # One extra due to quaternions
D = Differential(t)

# du = prob.f.f.f_oop(prob.u0, prob.p, 0)
# @test all(isfinite, du)

# prob = ODEProblem(ssys, ModelingToolkit.missing_variable_defaults(ssys), (0, 10))
prob = ODEProblem(ssys, [collect(body.v_0 .=> [0, 0, 0]); collect(body.w_a .=> [0, 1, 0]); ], (0, 10))
sol = solve(prob, Rodas5P(), u0 = prob.u0 .+ 1e-12 .* randn.())

doplot() &&
    plot(sol, idxs = [collect(body.r_0); collect(body.v_0)], layout = 6) |> display

@test sol(2pi, idxs = body.r_0[1])≈0 atol=1e-3
@test sol(2pi, idxs = body.r_0[2])≈0 atol=1e-3
@test sol(2pi, idxs = body.r_0[3])≈0 atol=1e-3
@test sol(pi, idxs = body.r_0[2]) < -2

# TODO: add more tests
ts = 0:0.1:2pi

ti = 5
for (i, ti) in enumerate(ts)
    R = get_R(sol, body.frame_a, ti)
    Q = sol(ti, idxs=body.Q)
    @test Multibody.from_Q(Q, 0).R ≈ R atol=1e-6 # Test that from_Q yields a rotation matrix consistent with R
    R2 = get_R(sol, spring.frame_b, ti)
    @test norm(R'R2 - I) < 1e-10
end

end




# ==============================================================================
## Simple motion with quaternions===============================================
# ==============================================================================
using LinearAlgebra, ModelingToolkit, Multibody, JuliaSimCompiler
using OrdinaryDiffEq, Test

@testset "Simple motion with quaternions and state in Body" begin

t = Multibody.t
@named joint = Multibody.FreeMotion(isroot = true, state=false)
@named body = Body(; m = 1, r_cm = [0.0, 0, 0], isroot=true, quat=true, w_a=[1,0.5,0.2])

# @named joint = Multibody.FreeMotion(isroot = true, state=true, quat=true)
# @named body = Body(; m = 1, r_cm = [0.0, 0, 0], isroot=false, w_a=[1,1,1])

world = Multibody.world


connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, body.frame_a)]


@named model = ODESystem(connections, t,
                         systems = [world, joint, body])
irsys = IRSystem(model)
ssys = structural_simplify(irsys)

D = Differential(t)
prob = ODEProblem(ssys, [collect(body.w_a) .=> [1,0,0];], (0, 2pi))

sol = solve(prob, Rodas4(); u0 = prob.u0 .+ 0 .* randn.())#, initializealg=ShampineCollocationInit(0.01))
# sol = solve(prob, Rodas4(); u0 = prob.u0 .+ 1e-12 .* randn.(), dtmin=1e-8, force_dtmin=true)
@test SciMLBase.successful_retcode(sol)
doplot() && plot(sol, layout=13)
# end

ts = 0:0.1:2pi
Q = Matrix(sol(ts, idxs = [body.Q...]))
Qh = Matrix(sol(ts, idxs = [body.Q̂...]))
n = Matrix(sol(ts, idxs = [body.n...]))
@test mapslices(norm, Qh, dims=1).^2 ≈ n
@test Q ≈ Qh ./ sqrt.(n) atol=1e-2
@test norm(mapslices(norm, Q, dims=1) .- 1) < 1e-2


ti = 5
for (i, ti) in enumerate(ts)
    R = get_R(sol, body.frame_a, ti)
    @test norm(R'R - I) < 1e-10
    @test Multibody.from_Q(Q[:, i], 0).R ≈ R atol=1e-3
end

# Test that rotational velocity of 1 results in one full rotation in 2π seconds. Test rotation around all three major axes
@test get_R(sol, body.frame_a, 0pi) ≈ I
@test get_R(sol, body.frame_a, pi/2) ≈ [1 0 0; 0 0 -1; 0 1 0]' atol=1e-3
@test get_R(sol, body.frame_a, 1pi) ≈ diagm([1, -1, -1]) atol=1e-3
@test get_R(sol, body.frame_a, 2pi) ≈ I atol=1e-3

prob = ODEProblem(ssys, [collect(body.w_a) .=> [0,1,0];], (0, 2pi))
sol = solve(prob, Rodas4())#, 
@test get_R(sol, body.frame_a, 0pi) ≈ I
@test get_R(sol, body.frame_a, 1pi) ≈ diagm([-1, 1, -1]) atol=1e-3
@test get_R(sol, body.frame_a, 2pi) ≈ I atol=1e-3


prob = ODEProblem(ssys, [collect(body.w_a) .=> [0,0,1];], (0, 2pi))
sol = solve(prob, Rodas4())#, 
@test get_R(sol, body.frame_a, 0pi) ≈ I
@test get_R(sol, body.frame_a, 1pi) ≈ diagm([-1, -1, 1]) atol=1e-3
@test get_R(sol, body.frame_a, 2pi) ≈ I atol=1e-3

end


# ============================================================
## Quaternion states in joint instead of body
# ============================================================

@testset "Quaternions and state in free motion" begin
    using LinearAlgebra, ModelingToolkit, Multibody, JuliaSimCompiler
    t = Multibody.t
    world = Multibody.world

    @named joint = Multibody.FreeMotion(isroot = true, state=true, quat=true)
    @named body = Body(; m = 1, r_cm = [0.0, 0, 0])

    connections = [connect(world.frame_b, joint.frame_a)
                connect(joint.frame_b, body.frame_a)]

    @named model = ODESystem(connections, t,
                            systems = [world, joint, body])
    irsys = IRSystem(model)
    ssys = structural_simplify(irsys)

    D = Differential(t)
    # q0 = randn(4); q0 ./= norm(q0)
    q0 = [1,0,0,0]
    prob = ODEProblem(ssys, [
        collect(body.w_a) .=> [1,1,1];
        collect(joint.Q) .=> q0;
        collect(joint.Q̂) .=> q0;
        ], (0, 2pi))

    using OrdinaryDiffEq, Test
    sol = solve(prob, Rodas4(); u0 = prob.u0 .+ 1e-6 .* randn.())
    @test SciMLBase.successful_retcode(sol)
    # doplot() && plot(sol, layout=21)


    ts = 0:0.1:2pi
    Q = Matrix(sol(ts, idxs = [joint.Q...]))
    Qh = Matrix(sol(ts, idxs = [joint.Q̂...]))
    n = Matrix(sol(ts, idxs = [joint.n...]))
    @test mapslices(norm, Qh, dims=1).^2 ≈ n
    @test Q ≈ Qh ./ sqrt.(n) atol=1e-2
    @test norm(mapslices(norm, Q, dims=1) .- 1) < 1e-2

    @test get_R(sol, joint.frame_b, 0pi) ≈ I
    @test_broken get_R(sol, joint.frame_b, 1pi) ≈ diagm([1, -1, -1]) atol=1e-3
    @test get_R(sol, joint.frame_b, 2pi) ≈ I atol=1e-3

    Matrix(sol(ts, idxs = [joint.w_rel_b...]))
end

# ============================================================
## Spherical joint pendulum with quaternions
# ============================================================

@testset "Spherical joint with quaternion state" begin
    using LinearAlgebra, ModelingToolkit, Multibody, JuliaSimCompiler
    t = Multibody.t
    world = Multibody.world


    @named joint = Multibody.Spherical(isroot=false, state=false, quat=false)
    @named rod = FixedTranslation(; r = [1, 0, 0])
    @named body = Body(; m = 1, isroot=true, quat=true, air_resistance=0.0)

    # @named joint = Multibody.Spherical(isroot=true, state=true, quat=true)
    # @named body = Body(; m = 1, r_cm = [1.0, 0, 0], isroot=false)

    connections = [connect(world.frame_b, joint.frame_a)
                connect(joint.frame_b, rod.frame_a)
                connect(rod.frame_b, body.frame_a)]


    @named model = ODESystem(connections, t,
                            systems = [world, joint, body, rod])
    irsys = IRSystem(model)
    @test_skip begin # https://github.com/JuliaComputing/JuliaSimCompiler.jl/issues/298
        ssys = structural_simplify(irsys)


        D = Differential(t)
        q0 = randn(4); q0 ./= norm(q0)
        # q0 = [1,0,0,0]
        prob = ODEProblem(ssys, [
            # collect(body.w_a) .=> [1,0,0];
            # collect(body.Q) .=> q0;
            # collect(body.Q̂) .=> q0;
            ], (0, 30))

        using OrdinaryDiffEq, Test
        sol = solve(prob, Rodas4(); u0 = prob.u0 .+ 0 .* randn.())
        @test SciMLBase.successful_retcode(sol)
        # doplot() && plot(sol, layout=21)


        ts = 0:0.1:2pi
        Q = Matrix(sol(ts, idxs = [body.Q...]))
        Qh = Matrix(sol(ts, idxs = [body.Q̂...]))
        n = Matrix(sol(ts, idxs = [body.n...]))
        @test mapslices(norm, Qh, dims=1).^2 ≈ n
        @test Q ≈ Qh ./ sqrt.(n) atol=1e-2
        @test norm(mapslices(norm, Q, dims=1) .- 1) < 1e-2

        Matrix(sol(ts, idxs = [body.w_a...]))

        @test get_R(joint.frame_b, 0pi) ≈ I
        @test get_R(joint.frame_b, sqrt(9.81/1)) ≈ diagm([1, -1, -1]) atol=1e-3
        @test get_R(joint.frame_b, 2pi) ≈ I atol=1e-3


        Matrix(sol(ts, idxs = [joint.w_rel_b...]))

        # render(model, sol)
    end
end


# ==============================================================================
## FreeBody ====================================================================
# ==============================================================================
# https://doc.modelica.org/om/Modelica.Mechanics.MultiBody.Examples.Elementary.FreeBody.html
using Multibody
using ModelingToolkit
# using Plots
using JuliaSimCompiler
using OrdinaryDiffEq

@testset "FreeBody" begin
    t = Multibody.t
    D = Differential(t)
    world = Multibody.world

    @named begin
        body = BodyShape(m = 1, I_11 = 1, I_22 = 1, I_33 = 1, r = [0.4, 0, 0],
                        r_0 = [0.2, -0.5, 0.1], r_cm = [0.2, 0, 0], isroot = true, quat=true)
        bar2 = FixedTranslation(r = [0.8, 0, 0])
        spring1 = Multibody.Spring(c = 20, s_unstretched = 0)
        spring2 = Multibody.Spring(c = 20, s_unstretched = 0)
    end

    eqs = [connect(bar2.frame_a, world.frame_b)
        connect(spring1.frame_b, body.frame_a)
        connect(bar2.frame_b, spring2.frame_a)
        connect(spring1.frame_a, world.frame_b)
        connect(body.frame_b, spring2.frame_b)]

    @named model = ODESystem(eqs, t,
                            systems = [
                                world,
                                body,
                                bar2,
                                spring1,
                                spring2,
                            ])
    ssys = structural_simplify(IRSystem(model))#, alias_eliminate = true)
    # ssys = structural_simplify(model, allow_parameters = false)
    prob = ODEProblem(ssys,
                    [world.g => 9.80665;
                        collect(body.body.w_a .=> 0);
                    collect(body.body.v_0 .=> 0);
                    #   collect(body.body.phi .=> deg2rad(10));
                        collect(body.body.Q) .=> vec(QuatRotation(RotXYZ(deg2rad.((10,10,10))...)));
                        collect(body.body.Q̂) .=> vec(QuatRotation(RotXYZ(deg2rad.((10,10,10))...)));
                    ], (0, 10))

    # @test_skip begin # The modelica example uses angles_fixed = true, which causes the body component to run special code for variable initialization. This is not yet supported by MTK
    sol = solve(prob, Rodas5P(), abstol=1e-6, reltol=1e-6)
    @test SciMLBase.successful_retcode(sol)

    @test sol(10, idxs=[
        body.r_0;
    ]) ≈ [
        0.196097
        -0.46243
        0.080562
    ] atol=1e-1

    @test_broken get_rot(sol, body.frame_b, 10)[1,2] ≈ 0.104409 atol=0.01

    doplot() && plot(sol, idxs = [body.r_0...; body.body.w_a; body.body.v_0], layout=(3,3), size=(1000, 1000)) |> display

end