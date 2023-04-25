using ModelingToolkit
using Multibody
using Test
using SymbolicIR
t = Multibody.t

world = Multibody.world

## Only body and world
@named body = Body(; m = 1, isroot = false, r_cm = [1, 0, 1])
Multibody.isroot(body)

connections = [
    connect(world.frame_b, body.frame_a),
]

@named model = ODESystem(connections, t, systems = [world, body])

modele = ModelingToolkit.expand_connections(model)

ssys = structural_simplify(model)

@test length(states(ssys)) == 0 # This example is completely rigid and should simplify down to zero state variables

# ==============================================================================
## Add spring to make a harmonic oscillator ====================================
# ==============================================================================
#=
The multibody paper mentions this as an interesting example, figure 8:
    "The non-standard feature to have potential states
    both in joints and in bodies is especially useful for
    inexperienced users, since they do not have to
    introduce a “virtual” joint with 6 degrees of
    freedom. For example, it is easy to just build up a
    system as in Figure 8, where a body is connected via
    a spring to the environment."
=#

# With 0 mass, the LineForce should have 

# !isroot
# = 2*(3+3+9+3) // 2*(r_0+f+R.T+t)
#  - 2*(3+3) // 2*(f+t)
#  - 2*(9-3) // 2*(R.T – R.residuals)
# = 12 equations 

# isroot for one of the frames
# = 2*(3+3+9+3) // 2*(r_0+f+R.T+t)
#  - 2*(3+3) // 2*(f+t)
#  - 1*(9-3) // 1*(R.T – R.residuals)
# = 18 equations 

@named body = Body(; m = 1, isroot = true, r_cm = [0, 1, 0], phi0 = [0, 1, 0]) # This time the body isroot since there is no joint containing state
@named spring = Multibody.Spring(c = 1, fixedRotationAtFrame_a = false,
                                 fixedRotationAtFrame_b = false)

connections = [connect(world.frame_b, spring.frame_a)
               connect(spring.frame_b, body.frame_a)]

@named model = ODESystem(connections, t, systems = [world, spring, body])

ModelingToolkit.n_extra_equations(model)

modele = ModelingToolkit.expand_connections(model)
# ssys = structural_simplify(model, allow_parameter = false)
# u0,p = ModelingToolkit.get_u0_p(ssys, [], [])

irsys = IRSystem(modele)
ssys = structural_simplify(irsys)
D = Differential(t)
defs = Dict(collect(spring.r_rel_0 .=> [0, 1, 0])...,
            collect(body.r_0 .=> [0, 0, 0])...,
            collect((D.(body.phi)) .=> [0, 0, 0])...,
            collect(D.(D.(body.phi)) .=> [0, 0, 0])...)
prob = ODEProblem(ssys, defs, (0, 10))

# du = prob.f.f.f_oop(prob.u0, prob.p, 0)
# @test all(isfinite, du)

using OrdinaryDiffEq
sol = solve(prob, Rodas5P())
@test SciMLBase.successful_retcode(sol)
@test sol(2pi, idxs = body.r_0[1])≈0 atol=1e-3
@test sol(2pi, idxs = body.r_0[2])≈1 atol=1e-3
@test sol(2pi, idxs = body.r_0[3])≈0 atol=1e-3
@test sol(pi, idxs = body.r_0[2]) < -2

isinteractive() &&
    plot(sol, idxs = [collect(body.r_0); collect(body.v_0); collect(body.phi)], layout = 9)

# ==============================================================================
## Simple pendulum =============================================================
# ==============================================================================
using LinearAlgebra
@named joint = Multibody.Revolute(n = [0, 0, 1], isroot = true)
@named body = Body(; m = 1, isroot = false, r_cm = [0.5, 0, 0])
@named torksensor = CutTorque()
@named forcesensor = CutForce()

connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, body.frame_a, torksensor.frame_a,
                       forcesensor.frame_a)]

connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, torksensor.frame_a)
               connect(torksensor.frame_b, forcesensor.frame_a)
               connect(forcesensor.frame_b, body.frame_a)]

@named model = ODESystem(connections, t,
                         systems = [world, joint, body, torksensor, forcesensor])
modele = ModelingToolkit.expand_connections(model)
ssys = structural_simplify(model, allow_parameter = false)

# irsys = IRSystem(modele)
# ssys = structural_simplify(irsys)

D = Differential(t)
defs = Dict(collect((D.(joint.phi)) .=> [0, 0, 0])...,
            collect(D.(D.(joint.phi)) .=> [0, 0, 0])...)
prob = ODEProblem(ssys, defs, (0, 10))

using OrdinaryDiffEq
sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
@test minimum(sol[joint.phi])≈-π rtol=0.01
@test maximum(sol[joint.phi])≈0 atol=0.01
@test all(x -> abs(x) < 1e-3, reduce(hcat, sol[torksensor.torque.u]))

@test maximum(norm.(eachcol(reduce(hcat, sol[forcesensor.force.u])))) ≈
      maximum(norm.(eachcol(reduce(hcat, sol[joint.frame_a.f]))))

isinteractive() && plot(sol, idxs = collect(joint.phi))

# ==============================================================================
## Simple pendulum from Modelica "First Example" tutorial ======================
# ==============================================================================

world = Multibody.world
@named body = Body(; m = 1, isroot = false, r_cm = [0.5, 0, 0])
@named damper = Rotational.Damper(d = 0.1)
@named rev = Multibody.Revolute(n = [0, 0, 1], useAxisFlange = true, isroot = true)

connections = [connect(world.frame_b, rev.frame_a)
               connect(damper.flange_b, rev.axis)
               connect(rev.support, damper.flange_a)
               connect(body.frame_a, rev.frame_b)]

@named model = ODESystem(connections, t, systems = [world, rev, body, damper])
modele = ModelingToolkit.expand_connections(model)
ssys = structural_simplify(model, allow_parameter = false)

# irsys = IRSystem(modele)
# ssys = structural_simplify(irsys)
D = Differential(t)
prob = ODEProblem(ssys, [damper.phi_rel => 1, D(rev.phi) => 0, D(D(rev.phi)) => 0],
                  (0, 100))

du = prob.f.f.f_oop(prob.u0, prob.p, 0)
@test all(isfinite, du)

using OrdinaryDiffEq
sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
@test minimum(sol[rev.phi]) > -π
@test sol[rev.phi][end]≈-π / 2 rtol=0.01 # pendulum settles at 90 degrees stable equilibrium
isinteractive() && plot(sol, idxs = collect(rev.phi))

# ==============================================================================
## Simple pendulum with rod ====================================================
# ==============================================================================

world = Multibody.world
@named rod = FixedTranslation(r = [0.5, 0, 0])
@named body = Body(; m = 1, isroot = false, r_cm = [0, 0, 0])
@named damper = Rotational.Damper(d = 0.1)
@named rev = Multibody.Revolute(n = [0, 0, 1], useAxisFlange = true, isroot = true)

connections = [connect(world.frame_b, rev.frame_a)
               connect(damper.flange_b, rev.axis)
               connect(rev.support, damper.flange_a)
               connect(rev.frame_b, rod.frame_a)
               connect(rod.frame_b, body.frame_a)]

@named model = ODESystem(connections, t, systems = [world, rev, body, damper, rod])
modele = ModelingToolkit.expand_connections(model)
# ssys = structural_simplify(model, allow_parameter = false)

ssys = structural_simplify(IRSystem(modele), alias_eliminate = false)

D = Differential(t)
prob = ODEProblem(ssys, [damper.phi_rel => 1, D(rev.phi) => 0, D(D(rev.phi)) => 0],
                  (0, 100))
sol2 = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol2)
@test minimum(sol2[rev.phi]) > -π
@test sol2[rev.phi][end]≈-π / 2 rtol=0.01 # pendulum settles at 90 degrees stable equilibrium
isinteractive() && plot(sol2, idxs = collect(rev.phi))

@info "TODO: write a test that checks that this solution is identical to the one without rod above"

# ==============================================================================
## Simple pendulum with rod from absolute rotation ====================================================
# ==============================================================================

world = Multibody.world
@named rod = FixedRotation(r = [0.5, 0, 0], n = [0, 0, 1], angle = 0)
@named body = Body(; m = 1, isroot = false, r_cm = [0, 0, 0])
@named damper = Rotational.Damper(d = 0.1)
@named rev = Multibody.Revolute(n = [0, 0, 1], useAxisFlange = true, isroot = true)

connections = [connect(world.frame_b, rev.frame_a)
               connect(damper.flange_b, rev.axis)
               connect(rev.support, damper.flange_a)
               connect(rev.frame_b, rod.frame_a)
               connect(rod.frame_b, body.frame_a)]

@named model = ODESystem(connections, t, systems = [world, rev, body, damper, rod])
modele = ModelingToolkit.expand_connections(model)

ssys = structural_simplify(model, allow_parameter = false)
# ssys = structural_simplify(IRSystem(modele)) # Yingbo, this fails with SymbolicIR but not with MTK

D = Differential(t)
@test_skip begin # Needs default dummy der
    prob = ODEProblem(ssys, [damper.phi_rel => 1, D(rev.phi) => 0, D(D(rev.phi)) => 0],
                      (0, 100), default_dummy_der = 0)
    sol2 = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol2)
    @test minimum(sol2[rev.phi]) > -π
    @test sol2[rev.phi][end]≈-π / 2 rtol=0.01 # pendulum settles at 90 degrees stable equilibrium
    isinteractive() && plot(sol2, idxs = collect(rev.phi))
end

# ==============================================================================
## Double pendulum =============================================================
# ==============================================================================

@named rod1 = FixedTranslation(r = [1, 0, 0])
@named rod2 = FixedTranslation(r = [1, 0, 0])
@named body1 = Body(; m = 1, isroot = false, r_cm = [0.0, 0, 0])
@named body2 = Body(; m = 1, isroot = false, r_cm = [0.0, 0, 0])
@named damper1 = Rotational.Damper(d = 0.1)
@named damper2 = Rotational.Damper(d = 0.1)
@named rev1 = Multibody.Revolute(n = [0, 0, 1], useAxisFlange = true, isroot = true)
@named rev2 = Multibody.Revolute(n = [0, 0, 1], useAxisFlange = true, isroot = true)

connections = [connect(damper1.flange_b, rev1.axis)
               connect(rev1.support, damper1.flange_a)
               connect(damper2.flange_b, rev2.axis)
               connect(rev2.support, damper2.flange_a)
               connect(world.frame_b, rev1.frame_a)
               connect(rev1.frame_b, rod1.frame_a)
               connect(rod1.frame_b, body1.frame_a)
               connect(body1.frame_a, rev2.frame_a)
               connect(rev2.frame_a, rod2.frame_a)
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
modele = ModelingToolkit.expand_connections(model)
ssys = structural_simplify(model, allow_parameter = false)

irsys = IRSystem(modele)
ssys = structural_simplify(irsys, alias_eliminate = false)
D = Differential(t)
prob = ODEProblem(ssys,
                  [
                      damper1.phi_rel => 1, D(rev1.phi) => 0, D(D(rev1.phi)) => 0,
                      damper2.phi_rel => 0.5, D(rev2.phi) => 0, D(D(rev2.phi)) => 0,
                      D(damper2.w_rel) => 0,
                  ], (0, 600))

sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
@test minimum(sol[rev1.phi]) > -π
@test sol[rev1.phi][end]≈-π / 2 rtol=0.01 # pendulum settles at 90 degrees stable equilibrium
@test_broken sol[rev2.phi][end]≈0 rtol=0.01 # pendulum settles at 90 degrees stable equilibrium which is 0 relative rotation compared to rev1
@test sol[body2.r_0[2]][end]≈-2 rtol=0.01 # sum of rod lengths = 2
isinteractive() &&
    plot(sol, idxs = [rev1.phi; rev2.phi; damper2.phi_rel; collect(body2.r_0[1:2])])

# ==============================================================================
## Linear mass-spring-damper ===================================================
# ==============================================================================

@named damper = Translational.Damper(0.5)
@named spring = Translational.Spring(1)
@named joint = Prismatic(n = [0, 1, 0], isroot = true, useAxisFlange = true)

connections = [connect(world.frame_b, joint.frame_a)
               connect(damper.flange_b, spring.flange_b, joint.axis)
               connect(joint.support, damper.flange_a, spring.flange_a)
               connect(body.frame_a, joint.frame_b)]

@named model = ODESystem(connections, t, systems = [world, joint, body, damper, spring])
ssys = structural_simplify(model, allow_parameter = false)

prob = ODEProblem(ssys, [damper.s_rel => 1, D(joint.s) => 0, D(D(joint.s)) => 0],
                  (0, 30))

sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
@test sol[joint.s][end]≈-9.81 rtol=0.01 # gravitational acceleration since spring stiffness is 1
isinteractive() && plot(sol, idxs = joint.s)

# ==============================================================================
## Spring damper system from https://www.maplesoft.com/documentation_center/online_manuals/modelica/Modelica_Mechanics_MultiBody_Examples_Elementary.html#Modelica.Mechanics.MultiBody.Examples.Elementary.SpringDamperSystem
# ==============================================================================

world = Multibody.world
@named begin
    body1 = Body(; m = 1, isroot = true, r_cm = [0.0, 0, 0], I_11 = 0.1, I_22 = 0.1,
                 I_33 = 0.1, r_0 = [0.3, -0.2, 0]) # This is root since there is no joint parallel to the spring leading to this body
    body2 = Body(; m = 1, isroot = false, r_cm = [0.0, -0.2, 0]) # This is not root since there is a joint parallel to the spring leading to this body
    bar1 = FixedTranslation(r = [0.3, 0, 0])
    bar2 = FixedTranslation(r = [0.6, 0, 0])
    p2 = Prismatic(n = [0, -1, 0], s0 = 0.1, useAxisFlange = true, isroot = true)
    spring2 = Multibody.Spring(c = 30, s_unstretched = 0.1)
    spring1 = Multibody.Spring(c = 30, s_unstretched = 0.1)
    damper1 = Multibody.Damper(d = 2)
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
       connect(spring1.frame_b, body1.frame_a)]

@named model = ODESystem(eqs, t,
                         systems = [
                             world,
                             body1,
                             body2,
                             bar1,
                             bar2,
                             p2,
                             spring1,
                             spring2,
                             damper1,
                         ])
# ssys = structural_simplify(model, allow_parameter = false)
ssys = structural_simplify(IRSystem(model), alias_eliminate = false)

prob = ODEProblem(ssys,
                  [collect(D.(body1.phid)) .=> 0;
                   D(p2.s) => 0;
                   D(D(p2.s)) => 0;
                   damper1.d => 0], (0, 10))

sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)

endpoint = sol(sol.t[end], idxs = [spring1.s, spring2.s])
@test_broken endpoint[1]≈endpoint[2] rtol=0.01

prob = ODEProblem(ssys,
                  [collect(D.(body1.phid)) .=> 0;
                   D(p2.s) => 0;
                   D(D(p2.s)) => 0;
                   damper1.d => 2], (0, 10))

sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)
@test sol(sol.t[end], idxs = spring1.v)≈0 atol=0.01 # damped oscillation

isinteractive() && plot(sol, idxs = [spring1.s, spring2.s])
isinteractive() && plot(sol, idxs = [body1.r_0[2], body2.r_0[2]])
isinteractive() && plot(sol, idxs = [spring1.f, spring2.f])

#=
model ThreeSprings "3-dim. springs in series and parallel connection"
  extends Modelica.Icons.Example;
  parameter Boolean animation=true "= true, if animation shall be enabled";
  inner Modelica.Mechanics.MultiBody.World world(animateWorld=animation);
  Modelica.Mechanics.MultiBody.Parts.Body body1(
    animation=animation,
    r_CM={0,-0.2,0},
    m=0.8,
    I_11=0.1,
    I_22=0.1,
    I_33=0.1,
    sphereDiameter=0.2,
    r_0(start={0.5,-0.3,0}, fixed=true),
    v_0(fixed=true),
    angles_fixed=true,
    w_0_fixed=true);
  Modelica.Mechanics.MultiBody.Parts.FixedTranslation bar1(animation=animation, r={0.3,0,0});
  Modelica.Mechanics.MultiBody.Forces.Spring spring1(
    lineForce(r_rel_0(start={-0.2,-0.2,0.2})),
    s_unstretched=0.1,
    width=0.1,
    coilWidth=0.005,
    numberOfWindings=5,
    c=20,
    animation=animation);
  Modelica.Mechanics.MultiBody.Parts.FixedTranslation bar2(animation=animation, r={0,0,0.3});
  Modelica.Mechanics.MultiBody.Forces.Spring spring2(
    s_unstretched=0.1,
    width=0.1,
    coilWidth=0.005,
    numberOfWindings=5,
    c=40,
    animation=animation);
  Modelica.Mechanics.MultiBody.Forces.Spring spring3(
    s_unstretched=0.1,
    width=0.1,
    coilWidth=0.005,
    numberOfWindings=5,
    c=20,
    animation=animation);
equation 
  connect(world.frame_b, bar1.frame_a);
  connect(world.frame_b, bar2.frame_a);
  connect(bar1.frame_b, spring1.frame_a);
  connect(bar2.frame_b, spring3.frame_a);
  connect(spring2.frame_b, body1.frame_a);
  connect(spring3.frame_b, spring1.frame_b);
  connect(spring2.frame_a, spring1.frame_b);
end ThreeSprings;
=#

using Multibody
using ModelingToolkit
using SymbolicIR
using OrdinaryDiffEq

# https://doc.modelica.org/om/Modelica.Mechanics.MultiBody.Examples.Elementary.ThreeSprings.html

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
                               fixedRotationAtFrame_a = true, fixedRotationAtFrame_b = true)
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
prob = ODEProblem(ssys,
                  [D.(collect(spring1.frame_b.r_0)) .=> 0;
                   D.(collect(spring3.frame_b.r_0)) .=> 0], (0, 10))

@test_skip begin # The modelica example uses angles_fixed = true, which causes the body component to run special code for variable initialization. This is not yet supported by MTK
    # Without proper initialization, the example fails most of the time. Random perturbation of u0 can make it work sometimes.
    sol = solve(prob, Rodas4(), u0 = prob.u0 .+ 1e-2 .* rand.())
    @test SciMLBase.successful_retcode(sol)

    isinteractive() && plot(sol, idxs = [body1.r_0...])
end
# TODO: add tutorial explaining what interesting things this demos illustrates
# fixedRotationAtFrame_a and b = true required

## FreeBody
#=
model FreeBody 
  "Free flying body attached by two springs to environment"
  extends Modelica.Icons.Example;
  parameter Boolean animation=true "= true, if animation shall be enabled";
  inner Modelica.Mechanics.MultiBody.World world;
  Modelica.Mechanics.MultiBody.Parts.FixedTranslation bar2(r={0.8,0,0}, animation=false);
  Modelica.Mechanics.MultiBody.Forces.Spring spring1(
    width=0.1,
    coilWidth=0.005,
    numberOfWindings=5,
    c=20,
    s_unstretched=0);
  Modelica.Mechanics.MultiBody.Parts.BodyShape body(
    m=1,
    I_11=1,
    I_22=1,
    I_33=1,
    r={0.4,0,0},
    r_CM={0.2,0,0},
    width=0.05,
    r_0(start={0.2,-0.5,0.1}, fixed=true),
    v_0(fixed=true),
    angles_fixed=true,
    w_0_fixed=true,
    angles_start={0.174532925199433,0.174532925199433,0.174532925199433});
  Modelica.Mechanics.MultiBody.Forces.Spring spring2(
    c=20,
    s_unstretched=0,
    width=0.1,
    coilWidth=0.005,
    numberOfWindings=5);
equation 
  connect(bar2.frame_a, world.frame_b);
  connect(spring1.frame_b, body.frame_a);
  connect(bar2.frame_b, spring2.frame_a);
  connect(spring1.frame_a, world.frame_b);
  connect(body.frame_b, spring2.frame_b);
end FreeBody;
=#

# https://doc.modelica.org/om/Modelica.Mechanics.MultiBody.Examples.Elementary.FreeBody.html
using Multibody
using ModelingToolkit
# using Plots
using SymbolicIR
using OrdinaryDiffEq

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
ssys = structural_simplify(IRSystem(model), alias_eliminate = true)
# ssys = structural_simplify(model, allow_parameters = false)
prob = ODEProblem(ssys,
                  [collect(D.(body.body.phid)) .=> 1;
                   collect(D.(body.body.phi)) .=> 1;
                   collect(D.(D.(body.body.phi))) .=> 1], (0, 10))

@test_skip begin # The modelica example uses angles_fixed = true, which causes the body component to run special code for variable initialization. This is not yet supported by MTK
    # Without proper initialization, the example fails most of the time. Random perturbation of u0 can make it work sometimes.
    sol = solve(prob, Rodas4(), u0 = prob.u0 .+ 1e-2 .* rand.())
    @assert SciMLBase.successful_retcode(sol)

    isinteractive() && plot(sol, idxs = [body.r_0...])
end

# ==============================================================================
## Sperical-joint pendulum ===================================================
# ==============================================================================
using Multibody
using ModelingToolkit
using Plots
using SymbolicIR
using OrdinaryDiffEq

t = Multibody.t
D = Differential(t)
world = Multibody.world

@named begin
    joint = Spherical(enforceStates = true, isroot = true)
    bar = FixedTranslation(r = [0, -1, 0])
    body = Body(; m = 1, isroot = false)
end

connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, bar.frame_a)
               connect(bar.frame_b, body.frame_a)]

@named model = ODESystem(connections, t, systems = [world, joint, bar, body])
# ssys = structural_simplify(model, allow_parameters = false)
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys,
                  [
                   # collect((body.phi)) .=> [0.5, 0.5, 0.5];
                   collect(D.(D.(joint.phi))) .=> 0;
                   collect(D.(joint.phi)) .=> 0
                   #    collect(D.(body.phid)) .=> 0
                   # collect(body.frame_a.w₃) .=> 0;
                   ], (0, 10))
@test_skip begin # Codegen leaves symbolic variables in function, Yingbo
    sol = solve(prob, Rodas4())
    @assert SciMLBase.successful_retcode(sol)

    plot(sol, idxs = [body.r_0...])
end

# ==============================================================================
## universal pendulum
# ==============================================================================

t = Multibody.t
D = Differential(t)
world = Multibody.world
@named begin
    joint = Universal()
    bar = FixedTranslation(r = [0, -1, 0])
    body = Body(; m = 1, isroot = false)
end
connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, bar.frame_a)
               connect(bar.frame_b, body.frame_a)]
@named model = ODESystem(connections, t, systems = [world, joint, bar, body])
ssys = structural_simplify(IRSystem(model))

@test_skip begin # NOTE: fails due to rotation matrix used as state
    prob = ODEProblem(ssys,
                      [joint.revolute_a.phi => 0;
                       D(joint.revolute_a.phi) => 0;
                       joint.revolute_b.phi => 0;
                       D(joint.revolute_b.phi) => 0], (0, 10))
    sol = solve(prob, Rodas4())
    @test SciMLBase.successful_retcode(sol)
    plot(sol, idxs = [body.r_0...])
end

# ==============================================================================
## GearConstraint ===================================================
# ==============================================================================
# https://doc.modelica.org/om/Modelica.Mechanics.MultiBody.Examples.Rotational3DEffects.GearConstraint.html

#=
model GearConstraint
   extends Modelica.Icons.Example;
  Joints.GearConstraint gearConstraint(                             ratio=10);
  inner World world(                             driveTrainMechanics3D=true,
      g=0);
  Parts.BodyCylinder cyl1(
    diameter=0.1,
    color={0,128,0},
    r={0.4,0,0});
  Parts.BodyCylinder cyl2(                             r={0.4,0,0}, diameter=
        0.2);
  Forces.Torque torque1;
  Blocks.Sources.Sine sine[         3](amplitude={2,0,0}, freqHz={1,1,1});
  Parts.Fixed fixed;
  Rotational.Components.Inertia inertia1(
      J=cyl1.I[1, 1],
    a(fixed=false),
    phi(fixed=true),
    w(fixed=true));
  Rotational.Components.IdealGear idealGear(        ratio=10, useSupport=true);
  Rotational.Components.Inertia inertia2(        J=cyl2.I[1, 1]);
  Rotational.Sources.Torque torque2(useSupport=true);
  Parts.Mounting1D mounting1D;
equation 
  connect(world.frame_b,gearConstraint. bearing);
  connect(cyl1.frame_b,gearConstraint. frame_a);
  connect(gearConstraint.frame_b,cyl2. frame_a);
  connect(torque1.frame_b,cyl1. frame_a);
  connect(torque1.frame_a,world. frame_b);
  connect(sine.y,torque1. torque);
  connect(inertia1.flange_b,idealGear. flange_a);
  connect(idealGear.flange_b,inertia2. flange_a);
  connect(torque2.flange,inertia1.   flange_a);
  connect(sine[1].y,torque2. tau);
  connect(mounting1D.flange_b,idealGear.support);
  connect(mounting1D.flange_b,torque2.support);
  connect(fixed.frame_b,mounting1D. frame_a);
end GearConstraint;
=#
##
using Multibody
using ModelingToolkit
using SymbolicIR
using OrdinaryDiffEq

t = Multibody.t
D = Differential(t)
world = Multibody.world

@named begin
    gearConstraint = GearConstraint(; ratio = 10)
    cyl1 = Body(; m = 1, r_cm = [0.4, 0, 0])
    cyl2 = Body(; m = 1, r_cm = [0.4, 0, 0])
    torque1 = Torque(resolveInFrame = :frame_b)
    # sine[1:3] = Blocks.Sine(frequency = 1)
    fixed = Fixed() # TODO: implement
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

@named model = ODESystem(eqs, t,
                         systems = [world;
                                    gearConstraint;
                                    cyl1;
                                    cyl2;
                                    torque1;
                                    # sine;
                                    fixed;
                                    inertia1;
                                    idealGear;
                                    inertia2;
                                    torque2;
                                    mounting1D])

# ssys = structural_simplify(model, allow_parameters=false)

@test_skip begin
    ssys = structural_simplify(IRSystem(model)) # Index out of bounds, Yingbo

    prob = ODEProblem(ssys,
                      [
                          D(gearConstraint.actuatedRevolute_b.phi) => 0,
                          D(inertia2.flange_a.phi) => 0,
                          D(D(idealGear.phi_b)) => 0,
                          D(gearConstraint.actuatedRevolute_a) => 0,
                      ], (0, 10))
end

# ==============================================================================
## Rolling wheel ===============================================================
# ==============================================================================
world = Multibody.world
@named wheel = RollingWheel(radius = 0.3, m = 2, I_axis = 0.06,
                            I_long = 0.12,
                            x0 = 0.2,
                            y0 = 0.2,
                            der_angles = [0, 5, 1])

cwheel = complete(wheel)
defs = [
    world.n => [0, 0, -1],
    collect(D.(cwheel.rollingWheel.angles)) => [0, 5, 1], # TODO: redundant since der_angles specified above, Yingbo
]

# ssys = structural_simplify(model, allow_parameters=false)

@test_skip begin # ERROR: AssertionError: ex isa Number Yingbo. MTK simplification works
    ssys = structural_simplify(IRSystem(wheel))
    prob = ODEProblem(ssys, defs, (0, 10))
end
