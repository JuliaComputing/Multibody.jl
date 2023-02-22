using ModelingToolkit
using Multibody
t = Multibody.t

world = Multibody.world

## Only body and world
@named body = Body(; m=1)

connections = [
    connect(world.frame_b, body.frame_a)
]

@named model = ODESystem(connections, t, systems=[world, body])

modele = ModelingToolkit.expand_connections(model)

ssys = structural_simplify(model)


## Add spring to make a harmonic oscillator
@named spring = Multibody.Spring(40)
connections = [
    connect(world.frame_b, spring.frame_a)
    connect(spring.frame_a, body.frame_a)
]

@named model = ODESystem(connections, t, systems=[world, spring, body])
ssys = structural_simplify(model)



## Simple pendulum

@named joint = Revolute()
