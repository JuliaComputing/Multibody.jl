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
@named spring = Multibody.Spring(1, fixedRotationAtFrame_a = false,
                                 fixedRotationAtFrame_b = false)

connections = [connect(world.frame_b, spring.frame_a)
               connect(spring.frame_b, body.frame_a)]

@named model = ODESystem(connections, t, systems = [world, spring, body])

ModelingToolkit.n_extra_equations(model)

modele = ModelingToolkit.expand_connections(model)
ssys = structural_simplify(model, allow_parameter = false)
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
@test sol(2pi, idxs=body.r_0[1]) ≈ 0 atol=1e-3
@test sol(2pi, idxs=body.r_0[2]) ≈ 1 atol=1e-3
@test sol(2pi, idxs=body.r_0[3]) ≈ 0 atol=1e-3
@test sol(pi, idxs=body.r_0[2]) > 2

isinteractive() && plot(sol, idxs = [collect(body.r_0); collect(body.v_0); collect(body.phi)], layout=9)

# ==============================================================================
## Simple pendulum =============================================================
# ==============================================================================

@named joint = Multibody.Revolute(n = [0, 0, 1], isroot = true)
@named body = Body(; m = 1, isroot = false, r_cm = [0.5, 0, 0])

connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, body.frame_a)]

@named model = ODESystem(connections, t, systems = [world, joint, body])
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
isinteractive() && plot(sol, idxs = collect(joint.phi))
@test maximum(sol[joint.phi])≈π rtol=0.01

# ==============================================================================
## Simple pendulum from Modelica "First Example" tutorial ======================
# ==============================================================================
using ModelingToolkitStandardLibrary.Mechanical.Rotational

world = Multibody.world
@named body = Body(; m = 1, isroot = false, r_cm = [0.5, 0, 0])
@named damper = Damper(d = 0.1)
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
isinteractive() && plot(sol, idxs = collect(rev.phi))
@test maximum(sol[rev.phi]) < π
@test sol[rev.phi][end]≈π / 2 rtol=0.01
