using ModelingToolkit
using Multibody
using Test
using JuliaSimCompiler
using OrdinaryDiffEq
t = Multibody.t
D = Differential(t)
doplot() = false
world = Multibody.world
W(args...; kwargs...) = Multibody.world

## Only body and world
@named body = Body(; m = 1, isroot = false, r_cm = [1, 0, 1])
Multibody.isroot(body)

connections = [
    connect(world.frame_b, body.frame_a),
]

@named model = ODESystem(connections, t, systems = [world, body])

modele = ModelingToolkit.expand_connections(model)

ssys = structural_simplify(model)

@test length(unknowns(ssys)) == 0 # This example is completely rigid and should simplify down to zero state variables

@testset "traj" begin
    @info "Testing traj"
    include("test_traj.jl")
end

@testset "robot" begin
    @info "Testing robot"
    include("test_robot.jl")
end

@testset "orientation_getters" begin
    @info "Testing orientation_getters"
    include("test_orientation_getters.jl")
end

@testset "quaternions" begin
    @info "Testing quaternions"
    include("test_quaternions.jl")
end

@testset "worldforces" begin
    @info "Testing worldforces"
    include("test_worldforces.jl")
end


# ==============================================================================
## Add spring to make a harmonic oscillator ====================================
# ==============================================================================
#=
The multibody paper mentions this as an interesting example, figure 8:
    "The non-standard feature to have potential state
    both in joints and in bodies is especially useful for
    inexperienced users, since they do not have to
    introduce a “virtual” joint with 6 degrees of
    freedom. For example, it is easy to just build up a
    system as in Figure 8, where a body is connected via
    a spring to the environment."
=#
t = Multibody.t
D = Differential(t)
@testset "spring - harmonic oscillator" begin

    @named body = Body(; m = 1, isroot = true, r_cm = [0, -1, 0], quat=true, neg_w=true) # This time the body isroot since there is no joint containing state
    @named spring = Multibody.Spring(c = 1)

    connections = [connect(world.frame_b, spring.frame_a)
                connect(spring.frame_b, body.frame_a)]

    @named model = ODESystem(connections, t, systems = [world, spring, body])

    # ssys = structural_simplify(model, allow_parameter = false)

    irsys = IRSystem(model)
    ssys = structural_simplify(irsys)
    D = Differential(t)

    # du = prob.f.f.f_oop(prob.u0, prob.p, 0)
    # @test all(isfinite, du)

    # @test_skip begin # Yingbo: instability
    prob = ODEProblem(ssys, [
        collect(body.r_0) .=> [0, -1e-5, 0]; # To make sure the spring has non-zero extent
        collect(body.w_a) .=> 0.00;
        collect(body.v_0) .=> 0;
    ], (0, 10))
    sol = solve(prob, Rodas5P(), u0 = prob.u0 .+ 0*1e-5 .* randn.())
    @test SciMLBase.successful_retcode(sol)
    @test sol(2pi, idxs = body.r_0[1])≈0 atol=1e-3
    @test sol(2pi, idxs = body.r_0[2])≈0 atol=1e-3
    @test sol(2pi, idxs = body.r_0[3])≈0 atol=1e-3
    @test sol(pi, idxs = body.r_0[2]) < -2

    doplot() &&
        plot(sol, idxs = [collect(body.r_0); collect(body.v_0)], layout = 6)
end

# ==============================================================================
## Simple pendulum =============================================================
# ==============================================================================
using LinearAlgebra, ModelingToolkit
@testset "Simple pendulum" begin
@named joint = Multibody.Revolute(n = [0, 0, 1], isroot = true, axisflange=true)
@named body = Body(; m = 1, isroot = false, r_cm = [0.5, 0, 0])
@named torksensor = CutTorque()
@named forcesensor = CutForce()
@named powersensor = Multibody.Power()
@named damper = Rotational.Damper(d = 1e-300)

connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, body.frame_a, torksensor.frame_a,
                       forcesensor.frame_a, powersensor.frame_a)]

connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, torksensor.frame_a)
               connect(torksensor.frame_b, forcesensor.frame_a)
               connect(damper.flange_a, joint.axis)
               connect(damper.flange_b, joint.support)
               connect(forcesensor.frame_b, powersensor.frame_a)
               connect(powersensor.frame_b, body.frame_a)]

@named model = ODESystem(connections, t,
                         systems = [world, joint, body, torksensor, forcesensor, powersensor, damper])
modele = ModelingToolkit.expand_connections(model)
# ssys = structural_simplify(model, allow_parameter = false)

irsys = IRSystem(modele)
ssys = structural_simplify(irsys)

D = Differential(t)
defs = Dict()
prob = ODEProblem(ssys, defs, (0, 10))

using OrdinaryDiffEq
sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
@test minimum(sol[joint.phi])≈-π rtol=0.01
@test maximum(sol[joint.phi])≈0 atol=0.01
@test all(x -> abs(x) < 1e-3, reduce(hcat, sol[collect(torksensor.torque.u)]))

@test maximum(norm.(eachcol(reduce(hcat, sol[collect(forcesensor.force.u)])))) ≈
      maximum(norm.(eachcol(reduce(hcat, sol[collect(joint.frame_a.f)]))))
@test norm(sol[powersensor.power.u]) < 1e-14
doplot() && plot(sol, idxs = collect(joint.phi))

# Test power sensos
defs = Dict(damper.d => 10)
prob = ODEProblem(ssys, defs, (0, 1))
sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
@test sol(1, idxs=powersensor.power.u) ≈ -1.94758 atol=1e-2
end

# ==============================================================================
## Point gravity ======================
# ==============================================================================
@mtkmodel PointGrav begin
    @components begin
        world = W()
        body1 = Body(
            m=1,
            I_11=0.1,
            I_22=0.1,
            I_33=0.1,
            r_0=[0,0.6,0],
            isroot=true,
            v_0=[1,0,0])
        body2 = Body(
            m=1,
            I_11=0.1,
            I_22=0.1,
            I_33=0.1,
            r_0=[0.6,0.6,0],
            isroot=true,
            v_0=[0.6,0,0])
    end
end
@named model = PointGrav()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
defs = [
    model.world.mu => 1
    model.world.point_gravity => true
    collect(model.body1.w_a) .=> 0
    collect(model.body2.w_a) .=> 0
    
]
prob = ODEProblem(ssys, defs, (0, 5))
sol = solve(prob, Rodas4())

@test sol(5, idxs=model.body2.r_0) ≈ [0.7867717, 0.478463, 0] atol=1e-1
# plot(sol)


# ==============================================================================
## Simple pendulum from Modelica "First Example" tutorial ======================
# ==============================================================================

@testset "Simple pendula" begin
world = Multibody.world
@named body = Body(; m = 1, isroot = false, r_cm = [0.5, 0, 0])
@named damper = Rotational.Damper(d = 0.1)
@named rev = Multibody.Revolute(n = [0, 0, 1], axisflange = true, isroot = true)

connections = [connect(world.frame_b, rev.frame_a)
               connect(damper.flange_b, rev.axis)
               connect(rev.support, damper.flange_a)
               connect(body.frame_a, rev.frame_b)]

@named model = ODESystem(connections, t, systems = [world, rev, body, damper])
modele = ModelingToolkit.expand_connections(model)
# ssys = structural_simplify(model, allow_parameter = false)

irsys = IRSystem(modele)
ssys = structural_simplify(irsys)
D = Differential(t)
prob = ODEProblem(ssys, [damper.phi_rel => 1, D(rev.phi) => 0, D(D(rev.phi)) => 0],
                  (0, 40))

# du = prob.f.f.f_oop(prob.u0, prob.p, 0)
# @test all(isfinite, du)

using OrdinaryDiffEq
sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
@test minimum(sol[rev.phi]) > -π
@test sol[rev.phi][end]≈-π / 2 rtol=0.01 # pendulum settles at 90 degrees stable equilibrium
doplot() && plot(sol, idxs = collect(rev.phi))

# ==============================================================================
## Simple pendulum with rod ====================================================
# ==============================================================================

world = Multibody.world
@named rod = FixedTranslation(r = [0.5, 0, 0])
@named body = Body(; m = 1, isroot = false, r_cm = [0, 0, 0])
@named damper = Rotational.Damper(d = 0.1)
@named rev = Multibody.Revolute(n = [0, 0, 1], axisflange = true, isroot = true)

connections = [connect(world.frame_b, rev.frame_a)
               connect(damper.flange_b, rev.axis)
               connect(rev.support, damper.flange_a)
               connect(rev.frame_b, rod.frame_a)
               connect(rod.frame_b, body.frame_a)]

@named model = ODESystem(connections, t, systems = [world, rev, body, damper, rod])
# modele = ModelingToolkit.expand_connections(model)
# ssys = structural_simplify(model, allow_parameter = false)

ssys = structural_simplify(IRSystem(model))#, alias_eliminate = false)

D = Differential(t)
prob = ODEProblem(ssys, [damper.phi_rel => 1, D(rev.phi) => 0, D(D(rev.phi)) => 0],
                  (0, 100))
sol2 = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol2)
@test minimum(sol2[rev.phi]) > -π
@test sol2[rev.phi][end]≈-π / 2 rtol=0.01 # pendulum settles at 90 degrees stable equilibrium
doplot() && plot(sol2, idxs = collect(rev.phi))

@test sol2(1:10, idxs=rev.phi).u ≈ sol(1:10, idxs=rev.phi).u atol=1e-2


# ==============================================================================
## Simple pendulum with rod from absolute rotation ====================================================
# ==============================================================================

world = Multibody.world
@named rod = FixedRotation(r = [0.5, 0, 0], n = [0, 0, 1], angle = 0)
@named body = Body(; m = 1, isroot = false, r_cm = [0, 0, 0])
@named damper = Rotational.Damper(d = 0.1)
@named rev = Multibody.Revolute(n = [0, 0, 1], axisflange = true, isroot = true)

connections = [connect(world.frame_b, rev.frame_a)
               connect(damper.flange_b, rev.axis)
               connect(rev.support, damper.flange_a)
               connect(rev.frame_b, rod.frame_a)
               connect(rod.frame_b, body.frame_a)]

@named model = ODESystem(connections, t, systems = [world, rev, body, damper, rod])
modele = ModelingToolkit.expand_connections(model)

# ssys = structural_simplify(model)#, allow_parameter = false)
ssys = structural_simplify(IRSystem(modele))

D = Differential(t)

prob = ODEProblem(ssys, [damper.phi_rel => 1], (0, 40))
sol3 = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol3)
@test minimum(sol3[rev.phi]) > -π
@test sol3[rev.phi][end]≈-π / 2 rtol=0.01 # pendulum settles at 90 degrees stable equilibrium
doplot() && plot(sol3, idxs = rev.phi)
@test sol3(1:10, idxs=rev.phi).u ≈ sol(1:10, idxs=rev.phi).u atol=1e-2
end
# ==============================================================================
## Double pendulum =============================================================
# ==============================================================================

using LinearAlgebra
@testset "Double pendulum" begin
@named rod1 = FixedTranslation(r = [1, 0, 0])
@named rod2 = FixedTranslation(r = [1, 0, 0])
@named body1 = Body(; m = 1, isroot = false, r_cm = [0.0, 0, 0])
@named body2 = Body(; m = 1, isroot = false, r_cm = [0.0, 0, 0])
@named damper1 = Rotational.Damper(d = 5)
@named damper2 = Rotational.Damper(d = 1)
@named rev1 = Multibody.Revolute(n = normalize([0.1, 0, 1]), axisflange = true, isroot = true)
@named rev2 = Multibody.Revolute(n = [0, 0, 1], axisflange = true, isroot = true)

connections = [connect(damper1.flange_b, rev1.axis)
               connect(rev1.support, damper1.flange_a)
               connect(damper2.flange_b, rev2.axis)
               connect(rev2.support, damper2.flange_a)
               connect(world.frame_b, rev1.frame_a)
               connect(rev1.frame_b, rod1.frame_a)
               connect(rod1.frame_b, body1.frame_a)
               connect(body1.frame_a, rev2.frame_a)
               connect(rev2.frame_b, rod2.frame_a)
               connect(rod2.frame_b, body2.frame_a)]

@named model = ODESystem(connections, t,
                         systems = [
                             world,
                             rev1,
                             rev2,
                             rod1,
                             rod2,
                             body1,
                             body2,
                             damper1,
                             damper2,
                         ])
# modele = ModelingToolkit.expand_connections(model)
# ssys = structural_simplify(model, allow_parameter = false)

irsys = IRSystem(model)
ssys = structural_simplify(irsys, alias_eliminate = false)
D = Differential(t)
prob = ODEProblem(ssys,
                  [
                      damper1.phi_rel => 1, D(rev1.phi) => 0, D(D(rev1.phi)) => 0,
                      damper2.phi_rel => 0.5, D(rev2.phi) => 0, D(D(rev2.phi)) => 0,
                      D(damper2.w_rel) => 0,
                  ], (0, 35))

sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
@test minimum(sol[rev1.phi]) > -π
@test sol[rev1.phi][end]≈-π / 2 rtol=0.01 # pendulum settles at 90 degrees stable equilibrium
@test_skip sol[rev2.phi][end]≈0 rtol=0.01 # pendulum settles at 90 degrees stable equilibrium which is 0 relative rotation compared to rev1
@test sol[body2.r_0[2]][end]≈-2 rtol=0.01 # sum of rod lengths = 2
doplot() &&
    plot(sol, idxs = [rev1.phi; rev2.phi; damper2.phi_rel; collect(body2.r_0[1:2])])
end
# ==============================================================================
## Linear mass-spring-damper ===================================================
# ==============================================================================

@testset "Linear mass-spring-damper" begin
@named body = Body(; m = 1, isroot = false, r_cm = [0, 0, 0])
@named damper = Translational.Damper(d=0.5)
@named spring = Translational.Spring(c=1)
@named joint = Prismatic(n = [0, 1, 0], axisflange = true)

connections = [connect(world.frame_b, joint.frame_a)
               connect(damper.flange_b, spring.flange_b, joint.axis)
               connect(joint.support, damper.flange_a, spring.flange_a)
               connect(body.frame_a, joint.frame_b)]

@named model = ODESystem(connections, t, systems = [world, joint, body, damper, spring])
ssys = structural_simplify(IRSystem(model))#, allow_parameter = false)

prob = ODEProblem(ssys, [damper.s_rel => 1, D(joint.s) => 0, D(D(joint.s)) => 0],
                  (0, 30))

sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
@test sol[joint.s][end]≈-9.81 rtol=0.01 # gravitational acceleration since spring stiffness is 1
doplot() && plot(sol, idxs = joint.s)
end
# ==============================================================================
## Spring damper system from https://www.maplesoft.com/documentation_center/online_manuals/modelica/Modelica_Mechanics_MultiBody_Examples_Elementary.html#Modelica.Mechanics.MultiBody.Examples.Elementary.SpringDamperSystem
# ==============================================================================

@testset "Spring damper system" begin
systems = @named begin
    world = W()
    body1 = Body(; m = 1, isroot = true, r_cm = [0.0, 0, 0], I_11 = 0.1, I_22 = 0.1,
                 I_33 = 0.1, r_0 = [0.3, -0.2, 0], quat=false) # This is root since there is no joint parallel to the spring leading to this body
    body2 = Body(; m = 1, isroot = false, r_cm = [0.0, -0.2, 0]) # This is not root since there is a joint parallel to the spring leading to this body
    body3 = Body(; m = 1, isroot = true, r_cm = [0.0, 0, 0], I_11 = 0.1, I_22 = 0.1,
                 I_33 = 0.1, r_0 = [1.8, -0.2, 0], quat=false)
    bar1 = FixedTranslation(r = [0.3, 0, 0])
    bar2 = FixedTranslation(r = [0.6, 0, 0])
    bar3 = FixedTranslation(r = [0.9, 0, 0])
    p2 = Prismatic(n = [0, -1, 0], s0 = 0.1, axisflange = true)
    spring2 = Multibody.Spring(c = 30, s_unstretched = 0.1)
    spring1 = Multibody.Spring(c = 30, s_unstretched = 0.1)
    damper1 = Multibody.Damper(d = 2)
    springdamper = SpringDamperParallel(c=30, d=2, s_unstretched = 0.1)
end
eqs = [connect(world.frame_b, bar1.frame_a)
       connect(bar1.frame_b, bar2.frame_a)
       connect(bar2.frame_b, p2.frame_a)
       connect(p2.frame_b, body2.frame_a)
       connect(bar2.frame_b, spring2.frame_a)
       connect(body2.frame_a, spring2.frame_b)
       connect(damper1.frame_a, bar1.frame_b)
       connect(spring1.frame_a, bar1.frame_b)
       connect(damper1.frame_b, body1.frame_a)
       connect(spring1.frame_b, body1.frame_a)
       
       connect(bar2.frame_b, bar3.frame_a)
       connect(bar3.frame_b, springdamper.frame_a)
       connect(springdamper.frame_b, body3.frame_a)
       ]

@named model = ODESystem(eqs, t; systems)
# ssys = structural_simplify(model, allow_parameter = false)
ssys = structural_simplify(IRSystem(model))#, alias_eliminate = false)

prob = ODEProblem(ssys,
                  [#collect(D.(body1.phid)) .=> 0;
                  collect(body1.v_0 .=> 0)
                  collect(body1.w_a .=> 0)
                  collect(body3.w_a .=> 0)
                  collect(body3.v_0 .=> 0)
                   damper1.d => 0], (0, 10)
)
du = similar(prob.u0)
prob.f(du, prob.u0, prob.p, 0)
sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)

endpoint = sol(sol.t[end], idxs = [spring1.s, spring2.s])
@test_broken endpoint[1]≈endpoint[2] rtol=0.01

prob = ODEProblem(ssys,
                  [collect(body1.v_0 .=> 0)
                  collect(body1.w_a .=> 0)
                  collect(body3.w_a .=> 0)
                  collect(body3.v_0 .=> 0)
                   damper1.d => 2], (0, 10))

sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
@test sol(sol.t[end], idxs = spring1.v)≈0 atol=0.01 # damped oscillation

@test norm([1 -1]*Matrix(sol(0:1:10, idxs=[spring1.s, springdamper.s]))) < 1e-10 # Test that the difference is small

doplot() && plot(sol, idxs = [spring1.s, spring2.s, springdamper.s])
doplot() && plot(sol, idxs = [body1.r_0[2], body2.r_0[2]])
doplot() && plot(sol, idxs = [spring1.f, spring2.f, springdamper.f])
end
# ==============================================================================
## Three springs ===============================================================
# ==============================================================================
# https://doc.modelica.org/om/Modelica.Mechanics.MultiBody.Examples.Elementary.ThreeSprings.html
using Multibody
using ModelingToolkit
using JuliaSimCompiler
using OrdinaryDiffEq


@testset "Three springs" begin
t = Multibody.t
D = Differential(t)
world = Multibody.world

@named begin
    body1 = Body(m = 0.8, I_11 = 0.1, I_22 = 0.1, I_33 = 0.1, r_0 = [0.5, -0.3, 0],
                 r_cm = [0, -0.2, 0], isroot = false)
    bar1 = FixedTranslation(r = [0.3, 0, 0])
    bar2 = FixedTranslation(r = [0, 0, 0.3])
    spring1 = Multibody.Spring(c = 20, m = 0, s_unstretched = 0.1,
                               r_rel_0 = [-0.2, -0.2, 0.2])
    spring2 = Multibody.Spring(c = 40, m = 0, s_unstretched = 0.1,
                               fixed_rotation_at_frame_a = true, fixed_rotation_at_frame_b = true)
    spring3 = Multibody.Spring(c = 20, m = 0, s_unstretched = 0.1)
end
eqs = [connect(world.frame_b, bar1.frame_a)
       connect(world.frame_b, bar2.frame_a)
       connect(bar1.frame_b, spring1.frame_a)
       connect(bar2.frame_b, spring3.frame_a)
       connect(spring2.frame_b, body1.frame_a)
       connect(spring3.frame_b, spring1.frame_b)
       connect(spring2.frame_a, spring1.frame_b)]

@named model = ODESystem(eqs, t,
                         systems = [
                             world,
                             body1,
                             bar1,
                             bar2,
                             spring1,
                             spring2,
                             spring3,
                         ])
ssys = structural_simplify(IRSystem(model))
# ssys = structural_simplify(model, allow_parameters = false)
prob = ODEProblem(ssys, [
    collect(body1.v_0 .=> 0);
], (0, 10))

# @test_skip begin # The modelica example uses angles_fixed = true, which causes the body component to run special code for variable initialization. This is not yet supported by MTK
# Without proper initialization, the example fails most of the time. Random perturbation of u0 can make it work sometimes.
sol = solve(prob, Rodas4())#, u0 = prob.u0 .+ 1e-1 .* rand.())  
@test SciMLBase.successful_retcode(sol)

doplot() && plot(sol, idxs = [body1.r_0...]) |> display
# end
# fixed_rotation_at_frame_a and b = true required
end

@testset "FreeBody" begin
t = Multibody.t
D = Differential(t)
world = Multibody.world

@named begin
    body = BodyShape(m = 1, I_11 = 1, I_22 = 1, I_33 = 1, r = [0.4, 0, 0],
                     r_0 = [0.2, -0.5, 0.1], r_cm = [0.2, 0, 0], isroot = true)
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
                  [collect(body.body.w_a .=> 0);
                  collect(body.body.v_0 .=> 0);], (0, 10))

# @test_skip begin # The modelica example uses angles_fixed = true, which causes the body component to run special code for variable initialization. This is not yet supported by MTK
# Without proper initialization, the example fails most of the time. Random perturbation of u0 can make it work sometimes.
sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)

doplot() && plot(sol, idxs = [body.r_0...]) |> display
# end

end

# ==============================================================================
## Planar joint ================================================================
# ==============================================================================
using LinearAlgebra
@mtkmodel PlanarTest begin
    @components begin
        world = W()
        planar = Planar(n=[0,0,1], n_x=[1,0,0])
        force = Force()
        body = Body(m=1)
    end
    @equations begin
        connect(world.frame_b, planar.frame_a, force.frame_a)
        connect(planar.frame_b, body.frame_a, force.frame_b)
        force.force.u[1] ~ sin(t)
        force.force.u[2] ~ t
        force.force.u[3] ~ t^2/2
        # force.force.u .~ [sin(t), t, t^2]
    end
end
@named sys = PlanarTest()
sys = complete(sys)
ssys = structural_simplify(IRSystem(sys))
prob = ODEProblem(ssys, [
    sys.world.g => 9.80665; # Modelica default
    # collect(sys.body.w_a) .=> 0;
    # collect(sys.body.v_0) .=> 0;
], (0, 2))

sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)

# plot(sol, idxs=sys.force.frame_a.f)
@test sol(2, idxs=sys.body.r_0) ≈ [1.0907, -18.28, 0] atol=1e-3

# plot(sol)

# ==============================================================================
## Sperical-joint pendulum ===================================================
# ==============================================================================
using Multibody
using ModelingToolkit
# using Plots
using JuliaSimCompiler
using OrdinaryDiffEq

@testset "Spherical-joint pendulum" begin
t = Multibody.t
D = Differential(t)
world = Multibody.world

@named begin
    joint = Spherical(state = true, isroot = true, phi = 1)
    bar = FixedTranslation(r = [0, -1, 0])
    body = Body(; m = 1, isroot = false)
end

connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, bar.frame_a)
               connect(bar.frame_b, body.frame_a)]

@named model = ODESystem(connections, t, systems = [world, joint, bar, body])
# ssys = structural_simplify(model, allow_parameters = false)
ssys = structural_simplify(IRSystem(model))
@test length(unknowns(ssys)) == 6


prob = ODEProblem(ssys,
                  [
                   # collect((body.phi)) .=> [0.5, 0.5, 0.5];
                   collect(D.(D.(joint.phi))) .=> 0;
                   collect(D.(joint.phi)) .=> 0
                   #    collect(D.(body.phid)) .=> 0
                   # collect(body.frame_a.w₃) .=> 0;
                   ], (0, 10))

sol = solve(prob, Rodas4())
@assert SciMLBase.successful_retcode(sol)

doplot() && plot(sol, idxs = [body.r_0...])

@test norm(sol(0, idxs = [body.r_0...])) ≈ 1
@test norm(sol(10, idxs = [body.r_0...])) ≈ 1

@test norm(sol(0, idxs = [joint.phi...])) ≈ √(3)
end
# ==============================================================================
## universal pendulum
# ==============================================================================
using LinearAlgebra
@testset "Universal pendulum" begin
t = Multibody.t
D = Differential(t)
world = Multibody.world
@named begin
    joint = Universal(length=0.1)
    bar = FixedTranslation(r = [0, -1, 0])
    body = Body(; m = 1, isroot = false, r_cm=[0.1, 0, 0])
end
connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, bar.frame_a)
               connect(bar.frame_b, body.frame_a)]
@named model = ODESystem(connections, t, systems = [world, joint, bar, body])
model = complete(model)
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys,
                  [
                    joint.phi_b => sqrt(2);
                   joint.revolute_a.phi => sqrt(2);
                   joint.revolute_b.phi => 0;
                   ], (0, 10))
sol2 = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol2)

@test norm(sol2(0, idxs = [body.r_0...])) ≈ 1
@test norm(sol2(10, idxs = [body.r_0...])) ≈ 1

doplot() && plot(sol2, idxs = [body.r_0...])
end
# ==============================================================================
## GearConstraint ===================================================
# ==============================================================================
# https://doc.modelica.org/om/Modelica.Mechanics.MultiBody.Examples.Rotational3DEffects.GearConstraint.html

using Multibody
using ModelingToolkit
using JuliaSimCompiler
using OrdinaryDiffEq

@testset "GearConstraint" begin
t = Multibody.t
D = Differential(t)
world = Multibody.world

systems = @named begin
    gearConstraint = GearConstraint(; ratio = 10)
    cyl1 = Body(; m = 1, r_cm = [0.4, 0, 0])
    cyl2 = Body(; m = 1, r_cm = [0.4, 0, 0])
    torque1 = Torque(resolve_frame = :frame_b)
    fixed = Multibody.Fixed() 
    inertia1 = Rotational.Inertia(J = cyl1.I_11)
    idealGear = Rotational.IdealGear(ratio = 10, use_support = true)
    inertia2 = Rotational.Inertia(J = cyl2.I_11)
    torque2 = Rotational.Torque(use_support = true)
    mounting1D = Mounting1D()
end

eqs = [connect(world.frame_b, gearConstraint.bearing)
       connect(cyl1.frame_a, gearConstraint.frame_a)
       connect(gearConstraint.frame_b, cyl2.frame_a)
       connect(torque1.frame_b, cyl1.frame_a)
       connect(torque1.frame_a, world.frame_b)
       # connect(sine.output, torque1.torque)
       torque1.torque.u .~ [2sin(t), 0, 0]
       connect(inertia1.flange_b, idealGear.flange_a)
       connect(idealGear.flange_b, inertia2.flange_a)
       connect(torque2.flange, inertia1.flange_a)
       # connect(sine.output, torque2.tau)
       torque2.tau.u ~ 2sin(t)
       connect(mounting1D.flange_b, idealGear.support)
       connect(mounting1D.flange_b, torque2.support)
       connect(fixed.frame_b, mounting1D.frame_a)]

@named model = ODESystem(eqs, t, systems = [world; systems])
cm = complete(model)
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, [
    cm.idealGear.phi_b => 0
    D(cm.idealGear.phi_b) => 0
], (0, 10))
sol = solve(prob, Rodas4())

@test sol[cm.inertia1.phi] ≈ 10*sol[cm.inertia2.phi] rtol = 1e-2
@test sol[cm.inertia1.flange_a.tau] ≈ 10*sol[cm.inertia2.flange_a.tau] rtol=1e-2
end

# ==============================================================================
## Rolling wheel ===============================================================
# ==============================================================================
using LinearAlgebra
# The wheel does not need the world
@testset "Rolling wheel" begin
@named wheel = RollingWheel(radius = 0.3, m = 2, I_axis = 0.06,
                        I_long = 0.12,
                        x0 = 0.2,
                        y0 = 0.2,
                        der_angles = [0, 5, 1])


wheel = complete(wheel)

defs = [
    collect(world.n .=> [0, 0, -1]);
    wheel.body.r_0[1] => 0.2;
    wheel.body.r_0[2] => 0.2;
    wheel.body.r_0[3] => 0.3;
    # collect(D.(cwheel.rollingWheel.angles)) .=> [0, 5, 1]
]

ssys = structural_simplify(IRSystem(wheel))
prob = ODEProblem(ssys, defs, (0, 0.8))
sol = solve(prob, RadauIIA5())
@test_broken SciMLBase.successful_retcode(sol)
@info "Write tests"

end

# ==============================================================================
## FreeMotion ==================================================================
# ==============================================================================

@testset "FreeMotion" begin
# Model a free-falling body
# Test 1, enforce state = false
world = Multibody.world
@named freeMotion = FreeMotion(state = false, isroot = false)
@named body = Body(m = 1, isroot = true)

eqs = [connect(world.frame_b, freeMotion.frame_a)
       connect(freeMotion.frame_b, body.frame_a)]

@named model = ODESystem(eqs, t,
                         systems = [world;
                                    freeMotion;
                                    body])
# ssys = structural_simplify(model, allow_parameters = false)
ssys = structural_simplify(IRSystem(model))
@test length(unknowns(ssys)) == 12

prob = ODEProblem(ssys, [world.g=>9.81; collect(body.w_a .=> [0, 0, 0]); collect(body.v_0 .=> [0, 0, 0]); ], (0, 10))

sol = solve(prob, Rodas4())
doplot() && plot(sol, idxs = body.r_0[2], title = "Free falling body")
y = sol(0:0.1:10, idxs = body.r_0[2])
@test y≈-9.81 / 2 .* (0:0.1:10) .^ 2 atol=1e-2 # Analytical solution to acceleration problem



## test 2 enforce state = true
@named freeMotion = FreeMotion(state = true, isroot = true)
@named body = Body(m = 1, isroot = false)

eqs = [connect(world.frame_b, freeMotion.frame_a)
       connect(freeMotion.frame_b, body.frame_a)]

@named model = ODESystem(eqs, t,
                         systems = [world;
                                    freeMotion;
                                    body])
# ssys = structural_simplify(model, allow_parameters = false)
ssys = structural_simplify(IRSystem(model))
@test length(unknowns(ssys)) == 12 

prob = ODEProblem(ssys, [world.g=>9.81; collect(body.w_a .=> [0, 1, 0]); collect(body.v_0 .=> [0, 0, 0]); ], (0, 10))

sol = solve(prob, Rodas4())
doplot() && plot(sol, idxs = body.r_0[2], title = "Free falling body")
y = sol(0:0.1:10, idxs = body.r_0[2])
@test y≈-9.81 / 2 .* (0:0.1:10) .^ 2 atol=1e-2 # Analytical solution to acceleration problem


# without joint
world = Multibody.world
@named body = Body(m = 1, isroot = true)

@named model = ODESystem([], t,
                         systems = [world;
                                    body])
# ssys = structural_simplify(model, allow_parameters = false)
ssys = structural_simplify(IRSystem(model))
@test length(unknowns(ssys)) == 12

prob = ODEProblem(ssys, [
    world.g=>9.81; 
    collect(body.w_a .=> [0, 1, 0]); 
    collect(body.v_0 .=> [0, 0, 0]); 
], (0, 10))

sol = solve(prob, Rodas4())
doplot() && plot(sol, idxs = body.r_0[2], title = "Free falling body")
y = sol(0:0.1:10, idxs = body.r_0[2])
@test y≈-9.81 / 2 .* (0:0.1:10) .^ 2 atol=1e-2 # Analytical solution to acceleration problem
end


# ==============================================================================
## Dzhanibekov effect ==========================================================
# ==============================================================================
@testset "Dzhanibekov effect" begin
world = Multibody.world
@named freeMotion = FreeMotion(state = false, isroot = false)
@named body = BodyShape(m = 1, r=0.1*[1,2,3], isroot = true, I_11 = 1, I_22 = 10, I_33 = 100, quat=false)

eqs = [connect(world.frame_b, freeMotion.frame_a)
       connect(freeMotion.frame_b, body.frame_a)]

@named model = ODESystem(eqs, t,
                         systems = [world;
                                    freeMotion;
                                    body])
# ssys = structural_simplify(model, allow_parameters = false)
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys,
                  Symbolics.scalarize.([
                    body.body.v_0 .=> [0,1,0] # Movement in y direction
                    body.body.w_a .=> [0.01, 3.0, 0.02] # Rotation around y axis
                    body.body.phi .=> 0.0001 .* randn.() 
                    # body.Q .=> 0.001 .* randn.() .+ [0,0,0,1]
                    world.g => 0.0 # In space
                  ]), 
                  (0, 3))

sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
doplot() && plot(sol, idxs = collect(body.body.phi), title = "Dzhanibekov effect") |> display
@info "Write tests"

@test sol(0, idxs = collect(body.body.phi)) != zeros(3) # The problem here is that the initial condition is completely ignored

end


## Actuated joint
using Multibody
using ModelingToolkit
using JuliaSimCompiler
using OrdinaryDiffEq
using ModelingToolkitStandardLibrary.Blocks
import ModelingToolkitStandardLibrary.Mechanical.Rotational

@testset "Actuated joint" begin
@parameters t
D = Differential(t)

# Workarounds for @mtkmodel bugs and limitations
RTorque = Rotational.Torque
BSine = Blocks.Sine
W(args...; kwargs...) = Multibody.world

@mtkmodel ActuatedJoint begin
    @components begin
        world = W()
        torque = RTorque()
        joint = Revolute(axisflange=true) # The axis flange provides an interface to the 1D torque input from ModelingToolkitStandardLibrary.Mechanical.Rotational
        torque_signal = BSine(frequency=1/5)
        body = BodyShape(; m = 1, r = [0.4, 0, 0])
    end
    @equations begin
        connect(world.frame_b, joint.frame_a)
        connect(joint.frame_b, body.frame_a)
        connect(torque_signal.output, torque.tau)
        connect(torque.flange, joint.axis)
    end
end

@named model = ActuatedJoint()
cm = complete(model)
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, [D(D(cm.joint.phi)) => 0], (0, 10))
sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
doplot() && plot(sol) |> display

# using GLMakie
# render(model, sol)


end


# ==============================================================================
## Rope pendulum ===============================================================
# ==============================================================================

# Stiff rope
world = Multibody.world
number_of_links = 6
@named rope = Multibody.Rope(l = 1, m = 1, n=number_of_links, c=0, d=0, air_resistance=0, d_joint=1)
@named body = Body(; m = 1, isroot = false, r_cm = [0, 0, 0])

connections = [connect(world.frame_b, rope.frame_a)
               connect(rope.frame_b, body.frame_a)]

@named stiff_rope = ODESystem(connections, t, systems = [world, body, rope])
# ssys = structural_simplify(model, allow_parameter = false)

@time "Simplify stiff rope pendulum" ssys = structural_simplify(IRSystem(stiff_rope))

D = Differential(t)
prob = ODEProblem(ssys, [
    collect(body.r_0) .=> [1,1,1];
    collect(body.w_a) .=> [1,1,1]
], (0, 5))
@time "Stiff rope pendulum" sol = solve(prob, Rodas4(autodiff=false); u0 = prob.u0 .+ 0.1);
@test SciMLBase.successful_retcode(sol)


if false
    import GLMakie
    @time "render stiff_rope" render(stiff_rope, sol) # Takes very long time for n>=8
end



## Flexible rope
world = Multibody.world
number_of_links = 3
@named rope = Multibody.Rope(l = 1, m = 5, n=number_of_links, c=800.0, d=0.001, d_joint=0.1, air_resistance=0.2)
@named body = Body(; m = 300, isroot = false, r_cm = [0, 0, 0], air_resistance=0)

connections = [connect(world.frame_b, rope.frame_a)
               connect(rope.frame_b, body.frame_a)]

@named flexible_rope = ODESystem(connections, t, systems = [world, body, rope])
# ssys = structural_simplify(model, allow_parameter = false)

@time "Simplify flexible rope pendulum" ssys = structural_simplify(IRSystem(flexible_rope))
D = Differential(t)
prob = ODEProblem(ssys, [
    # D.(D.(collect(rope.r))) .=> 0;
    collect(body.r_0) .=> [1,1,1];
    collect(body.w_a) .=> [1,1,1];
    collect(body.v_0) .=> [10,10,10]
], (0, 10))
@time "Flexible rope pendulum" sol = solve(prob, Rodas4(autodiff=false); u0 = prob.u0 .+ 0.5);
@test SciMLBase.successful_retcode(sol)
if false
    import GLMakie
    @time "render flexible_rope" render(flexible_rope, sol) # Takes very long time for n>=8
end




## Resistance in spherical joint
# This test creates a spherical pendulum and one simple one with Revolute joint. The damping should work the same


systems = @named begin
    joint = Spherical(state=true, isroot=true, phi = [π/2, 0, 0], d = 0.3)
    bar = FixedTranslation(r = [0, -1, 0])
    body = Body(; m = 1, isroot = false)


    bartop = FixedTranslation(r = [1, 0, 0])
    joint2 = Multibody.Revolute(n = [1, 0, 0], axisflange = true, isroot = true)
    bar2 = FixedTranslation(r = [0, 1, 0])
    body2 = Body(; m = 1, isroot = false)
    damper = Rotational.Damper(d = 0.3)
end

connections = [connect(world.frame_b, joint.frame_a)
            connect(joint.frame_b, bar.frame_a)
            connect(bar.frame_b, body.frame_a)
            
            connect(world.frame_b, bartop.frame_a)
            connect(bartop.frame_b, joint2.frame_a)
            connect(joint2.frame_b, bar2.frame_a)
            connect(bar2.frame_b, body2.frame_a)

            connect(joint2.support, damper.flange_a)
            connect(damper.flange_b, joint2.axis)
            ]

@named model = ODESystem(connections, t, systems = [world; systems])
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys, [
                    D.(joint.phi) .=> 0;
                    D.(D.(joint.phi)) .=> 0;
                    joint2.phi => π/2
], (0, 10))

sol = solve(prob, Rodas4())
@assert SciMLBase.successful_retcode(sol)

# plot(sol, idxs = [body.r_0;], layout=3)
# plot!(sol, idxs = [body2.r_0; ], sp=[1 2 3])
# render(model, sol)

tt = 0:0.1:10
@test Matrix(sol(tt, idxs = [collect(body.r_0[2:3]);])) ≈ Matrix(sol(tt, idxs = [collect(body2.r_0[2:3]);]))


@test_skip begin # Produces state with rotation matrix
    number_of_links = 3
    chain_length = 2
    x_dist = 1.5 # Distance between the two mounting points
    systems = @named begin
        chain = Rope(l = chain_length, m = 5, n=number_of_links, c=1, d_joint=0.2, dir=[1, 0, 0], color=[0.5, 0.5, 0.5, 1], radius=0.05, cutprismatic=false, cutspherical=true)
        fixed = FixedTranslation(; r=[x_dist, 0, 0], radius=0.02, color=[0.1,0.1,0.1,1]) # Second mounting point
    end

    connections = [connect(world.frame_b, fixed.frame_a, chain.frame_a)
                connect(chain.frame_b, fixed.frame_b)]

    @named mounted_chain = ODESystem(connections, t, systems = [systems; world])

    ssys = structural_simplify(IRSystem(mounted_chain))
    prob = ODEProblem(ssys, [
        collect(chain.link_3.body.w_a) .=> [0,0,0]; 
        collect(chain.link_3.frame_b.r_0) .=> [x_dist,0,0]; 
    ], (0, 4))
    sol = solve(prob, Rodas4(autodiff=false))
    @test SciMLBase.successful_retcode(sol)

    # Multibody.render(mounted_chain, sol, x=3, filename = "mounted_chain.gif") # May take long time for n>=10
end


# ==============================================================================
## BodyCylinder ================================================================
# ==============================================================================
using LinearAlgebra

@testset "BodyCylinder" begin
    @info "Testing BodyCylinder"
    world = Multibody.world
    @mtkmodel CylinderPend begin
        @components begin
            world = W()
            body = BodyCylinder(r=[1,2,3], diameter=0.1)
            joint = Revolute()
        end
        @equations begin
            connect(world.frame_b, joint.frame_a)
            connect(joint.frame_b, body.frame_a)
        end
    end

    @named model = CylinderPend()
    model = complete(model)
    ssys = structural_simplify(IRSystem(model))

    prob = ODEProblem(ssys, [model.joint.phi => 0], (0, 10))
    sol = solve(prob, Rodas5P(), abstol=1e-8, reltol=1e-8)
    @test sol(10, idxs=model.body.body.m) ≈ 226.27 rtol=1e-3 # Values from open modelica
    @test sol(10, idxs=model.body.body.I_11) ≈ 245.28 rtol=1e-3
    @test sol(10, idxs=model.body.body.I_22) ≈ 188.74 rtol=1e-3
    @test sol(10, idxs=model.body.body.I_33) ≈ 94.515 rtol=1e-3
    @test sol(10, idxs=model.body.body.I_21) ≈ -37.69 rtol=1e-3
    @test sol(10, idxs=model.body.body.I_31) ≈ -56.53 rtol=1e-3
    @test sol(10, idxs=model.body.body.I_32) ≈ -113 rtol=1e-3
    @test sol(10, idxs=model.joint.phi) ≈ -2.1036 atol=1e-2
    # using Plots; plot(sol)


    prob = ODEProblem(ssys, [model.joint.phi => 0; model.body.inner_diameter=>0.05], (0, 10))
    sol = solve(prob, Rodas5P(), abstol=1e-8, reltol=1e-8)
    @test sol(10, idxs=model.body.body.m) ≈ 169.7 rtol=1e-3 # Values from open modelica
    @test sol(10, idxs=model.joint.phi) ≈ -2.0992 atol=1e-2
    @test sol(10, idxs=model.body.body.I_31) ≈ -42.39 rtol=1e-3
# using Plots; plot(sol)
end

#

# ==============================================================================
## BodyBox =====================================================================
# ==============================================================================
using LinearAlgebra

# @testset "BodyBox" begin
    # NOTE: r = [0,1,0] yields unstable simulation due to the commented branch in from_nxy: if n_z_aux' * n_z_aux > 1.0e-6
    # NOTE: for r=[0,0,1], r_shape=[0.1, 0, 0] the render of the box appears to have negative gravity
    @info "Testing BodyBox"
    @mtkmodel BoxPend begin
        @components begin
            world = W()
            body = Multibody.BodyBox(r=[0.1, 1, 0.2], r_shape=[0, 0, 0], width=0.1, height=0.3, inner_width=0.05)
            joint = Revolute()
        end
        @equations begin
            connect(world.frame_b, joint.frame_a)
            connect(joint.frame_b, body.frame_a)
        end
    end

    @named model = BoxPend()
    model = complete(model)
    ssys = structural_simplify(IRSystem(model))

    prob = ODEProblem(ssys, [model.joint.phi => 0], (0, 1))
    sol = solve(prob, Rodas5P(), abstol=1e-8, reltol=1e-8)
    # first(render(model, sol, 0, x=2.5, y=1.5, z=2.5, show_axis=true))
    # @test sol(10, idxs=model.body.body.m) ≈ 226.27 rtol=1e-3 # Values from open modelica
    # @test sol(10, idxs=model.body.body.I_11) ≈ 245.28 rtol=1e-3
    # @test sol(10, idxs=model.body.body.I_22) ≈ 188.74 rtol=1e-3
    # @test sol(10, idxs=model.body.body.I_33) ≈ 94.515 rtol=1e-3
    # @test sol(10, idxs=model.body.body.I_21) ≈ -37.69 rtol=1e-3
    # @test sol(10, idxs=model.body.body.I_31) ≈ -56.53 rtol=1e-3
    # @test sol(10, idxs=model.body.body.I_32) ≈ -113 rtol=1e-3
    # @test sol(10, idxs=model.joint.phi) ≈ -2.1036 atol=1e-2
    # # using Plots; plot(sol)


    # prob = ODEProblem(ssys, [model.joint.phi => 0; model.body.inner_width=>0.05], (0, 10))
    # sol = solve(prob, Rodas5P(), abstol=1e-8, reltol=1e-8)
    # @test sol(10, idxs=model.body.body.m) ≈ 169.7 rtol=1e-3 # Values from open modelica
    # @test sol(10, idxs=model.joint.phi) ≈ -2.0992 atol=1e-2
    # @test sol(10, idxs=model.body.body.I_31) ≈ -42.39 rtol=1e-3
# using Plots; plot(sol)
# end

##

using LinearAlgebra, ModelingToolkit, Multibody, JuliaSimCompiler, OrdinaryDiffEq
using Multibody.Rotations: RotXYZ
t = Multibody.t
D = Multibody.D
world = Multibody.world

@named joint = Multibody.Spherical(isroot=false, state=false, quat=false)
@named rod = FixedTranslation(; r = [1, 0, 0])
@named body = Body(; m = 1, isroot=true, quat=true, neg_w=true)

connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, rod.frame_a)
               connect(rod.frame_b, body.frame_a)]

@named model = ODESystem(connections, t,
                         systems = [world, joint, body, rod])
irsys = IRSystem(model)
ssys = structural_simplify(irsys)
prob = ODEProblem(ssys, [
    # vec(ori(rod.frame_a).R) .=> vec(RotXYZ(0,0,0));
    # D.(body.Q̂) .=> 0;

], (0, 1))
sol1 = solve(prob, FBDF(), abstol=1e-8, reltol=1e-8)

## quat in joint
@named joint = Multibody.Spherical(isroot=true, state=true, quat=true, neg_w=true)
@named rod = FixedTranslation(; r = [1, 0, 0])
@named body = Body(; m = 1, isroot=false, quat=false)

connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, rod.frame_a)
               connect(rod.frame_b, body.frame_a)]

@named model = ODESystem(connections, t,
                         systems = [world, joint, body, rod])
irsys = IRSystem(model)
ssys = structural_simplify(irsys)
prob = ODEProblem(ssys, [
    # vec(ori(rod.frame_a).R) .=> vec(RotXYZ(0,0,0));
    # D.(joint.Q̂) .=> 0;

], (0, 1))
sol2 = solve(prob, FBDF(), abstol=1e-8, reltol=1e-8)

## euler
@named joint = Multibody.Spherical(isroot=true, state=true, quat=false, neg_w=true)
@named rod = FixedTranslation(; r = [1, 0, 0])
@named body = Body(; m = 1, isroot=false, quat=false)

connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, rod.frame_a)
               connect(rod.frame_b, body.frame_a)]

@named model = ODESystem(connections, t,
                         systems = [world, joint, body, rod])
irsys = IRSystem(model)
ssys = structural_simplify(irsys)
prob = ODEProblem(ssys, [
    # vec(ori(rod.frame_a).R) .=> vec(RotXYZ(0,0,0));
    # D.(joint.Q̂) .=> 0;

], (0, 1))
sol3 = solve(prob, FBDF(), abstol=1e-8, reltol=1e-8)

@test Matrix(sol1(0:0.1:1, idxs=collect(body.r_0))) ≈ Matrix(sol2(0:0.1:1, idxs=collect(body.r_0))) rtol=1e-5
@test Matrix(sol1(0:0.1:1, idxs=collect(body.r_0))) ≈ Matrix(sol3(0:0.1:1, idxs=collect(body.r_0))) rtol=1e-5


# ==============================================================================
## SphericalSpherical
# ==============================================================================

@mtkmodel TestSphericalSpherical begin
    @components begin
        world = W()
        ss = SphericalSpherical(r_0 = [1, 0, 0], m = 1, kinematic_constraint=false)
        ss2 = BodyShape(r = [0, 0, 1], m = 1, isroot=true)
        s = Spherical()
        trans = FixedTranslation(r = [1,0,1])
    end
    @equations begin
        connect(world.frame_b, ss.frame_a, trans.frame_a)
        connect(ss.frame_b, ss2.frame_a)
        connect(ss2.frame_b, s.frame_a)
        connect(s.frame_b, trans.frame_b)
    end
end

@named model = TestSphericalSpherical()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, [
    model.ss2.body.phi[1] => 0.1;
    model.ss2.body.phid[3] => 0.0;
], (0, 1.37))
sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
# plot(sol)

# ==============================================================================
## UniversalSpherical
# ==============================================================================

@mtkmodel TestSphericalSpherical begin
    @components begin
        world = W()
        ss = UniversalSpherical(rRod_ia = [1, 0, 0], kinematic_constraint=false, sphere_diameter=0.3)
        ss2 = BodyShape(r = [0, 0, 1], m = 1, isroot=true, neg_w=true)
        s = Spherical()
        trans = FixedTranslation(r = [1,0,1])
        body2 = Body(; m = 1, isroot = false, r_cm=[0.1, 0, 0], neg_w=true)
    end
    @equations begin
        connect(world.frame_b, ss.frame_a, trans.frame_a)
        connect(ss.frame_b, ss2.frame_a)
        connect(ss2.frame_b, s.frame_a)
        connect(s.frame_b, trans.frame_b)
        connect(ss.frame_ia, body2.frame_a)
    end
end

@named model = TestSphericalSpherical()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, [
    model.ss2.body.phi[1] => 0.1;
    model.ss2.body.phid[3] => 0.0;
], (0, 1.37))
sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
# plot(sol)

## =============================================================================
@testset "fourbar" begin
    @info "Testing fourbar"
    include("test_fourbar.jl")
end

## =============================================================================