using ModelingToolkit
using Multibody
using ModelingToolkitStandardLibrary.Mechanical.TranslationalPosition
t = Multibody.t

world = Multibody.world
@named spring = Spring(k=40)
@named body = Body(; m=1)

connections = [
    connect(world.frame_b, spring.flange_a)
    connect(spring.flange_a, body.frame_a)
]

@named model = ODESystem(connections, t, systems=[world, spring, body])