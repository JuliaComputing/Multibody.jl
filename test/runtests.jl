using ModelingToolkit
using Multibody
t = Multibody.t

world = Multibody.world

## Only body and world
@named body = Body(; m=1, isroot=false)
Multibody.isroot(body)

connections = [
    connect(world.frame_b, body.frame_a)
]

@named model = ODESystem(connections, t, systems=[world, body])

modele = ModelingToolkit.expand_connections(model)

ssys = structural_simplify(model)

@test length(states(ssys)) == 0

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
@named spring = Multibody.Spring(40)

connections = [
    connect(world.frame_b, spring.frame_a)
    connect(spring.frame_b, body.frame_a)
]

@named model = ODESystem(connections, t, systems=[world, spring, body])
modele = ModelingToolkit.expand_connections(model)
ssys = structural_simplify(model)



## Simple pendulum

@named joint = Multibody.Revolute()

connections = [
    connect(world.frame_b, joint.frame_a)
    connect(joint.frame_b, body.frame_a)
]

@named model = ODESystem(connections, t, systems=[world, joint, body])
modele = ModelingToolkit.expand_connections(model)
ssys = structural_simplify(model)

##
