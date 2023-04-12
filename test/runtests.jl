using ModelingToolkit
using Multibody
t = Multibody.t

world = Multibody.world

## Only body and world
@named body = Body(; m=1, isroot=false, r_cm=[1,0,1])
Multibody.isroot(body)

connections = [
    connect(world.frame_b, body.frame_a)
]

@named model = ODESystem(connections, t, systems=[world, body])

modele = ModelingToolkit.expand_connections(model)

ssys = structural_simplify(model)

@test length(states(ssys)) == 0 # This example is completely rigid and should simplify down to zero state variables

@error("There is a problem with how angular velocities are handled when calling Multibody.ori. When Body.isroot, the body must have additional states corresponding to angular velocities, whereas if it's not root, the states are defined elsewhere. ori just neglects w and recreates from R through differentiation")



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

@named body = Body(; m=1, isroot=true, r_cm=[1,0,1], ϕ0 = [1,1,1]) # This time the body isroot since there is no joint containing state

@named lineForceBase = Multibody.LineForceBase(; length = 0)
@named lineForce = Multibody.LineForceWithMass(; length = 0, m=0, lengthFraction=0.5)
@named spring = Multibody.Spring(40, fixedRotationAtFrame_a=false, fixedRotationAtFrame_b=false)

connections = [
    connect(world.frame_b, spring.frame_a)
    connect(spring.frame_b, body.frame_a)
]

@named model = ODESystem(connections, t, systems=[world, spring, body])

ModelingToolkit.n_extra_equations(model)

modele = ModelingToolkit.expand_connections(model)
ssys = structural_simplify(model)
u0,p = ModelingToolkit.get_u0_p(ssys, [], [])

# prob = ODEProblem(ssys, [body.v_0 => randn(3)], (0, 10))
prob = ODEProblem(ssys, [], (0, 10))

du = prob.f.f.f_oop(u0, p, 0)
@test_broken all(isfinite, du)



using OrdinaryDiffEq
sol = solve(prob, Rodas4())

@test_broken SciMLBase.successful_retcode(sol) # Fails to initialize

# ==============================================================================
## Simple pendulum =============================================================
# ==============================================================================

@named joint = Multibody.Revolute()
@named body = Body(; m=1, isroot=false, r_cm=[1,0,1])

connections = [
    connect(world.frame_b, joint.frame_a)
    connect(joint.frame_b, body.frame_a)
]

@named model = ODESystem(connections, t, systems=[world, joint, body])
modele = ModelingToolkit.expand_connections(model)
ssys = structural_simplify(model)

prob = ODEProblem(ssys, [], (0, 10))

u0,p = ModelingToolkit.get_u0_p(ssys, [], [])

du = prob.f.f.f_oop(u0, p, 0)
@test_broken all(isfinite, du)


using OrdinaryDiffEq
sol = solve(prob, Rodas4())


##
ˍ₋arg1 = u0
ˍ₋arg2 = p
