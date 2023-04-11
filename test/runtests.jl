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

## Add spring to make a harmonic oscillator


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

@named lineForceBase = Multibody.LineForceBase(; length = 0)
@named lineForce = Multibody.LineForceWithMass(; length = 0, m=0, lengthFraction=0.5)
@named spring = Multibody.Spring(40, fixedRotationAtFrame_a=false, fixedRotationAtFrame_b=true)

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



using OrdinaryDiffEq
sol = solve(prob, Rodas4())

@test_broken SciMLBase.successful_retcode(sol) # Fails to initialize

## Simple pendulum

@named joint = Multibody.Revolute()

connections = [
    connect(world.frame_b, joint.frame_a)
    connect(joint.frame_b, body.frame_a)
]

@named model = ODESystem(connections, t, systems=[world, joint, body])
modele = ModelingToolkit.expand_connections(model)
ssys = structural_simplify(model)

prob = ODEProblem(ssys, [], (0, 10))

u0,p = ModelingToolkit.get_u0_p(ssys, [], [])


using OrdinaryDiffEq
sol = solve(prob, Rodas4())


##
ˍ₋arg1 = u0
ˍ₋arg2 = p
