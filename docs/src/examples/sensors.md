# Using a sensor
Multibody models are composed out of multibody components, typically by `connect`ing [`Frames`](@ref) together. To add, e.g., a control system to a multibody model, we often make use of sensors that, just like in the real world, translate between a physical quantity in 3D and a signal in the "signal domain". The signal domain consists of the various blocks defined in [`ModelingToolkitStandardLibrary.Blocks`](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/blocks/), and a sensor is simply a component that has two connectors, one [`Frame`](@ref) connector for attaching to the multibody model, and one [`Blocks.RealOutput`](@ref) for connecting to the signal domain.


The example below adds a force and a torque sensor to the pivot point of a pendulum. Note how the two sensors are connected in series with each other, just like how we would typically connect them in practice if they are not integrated into the same component.
```@example sensor
using Multibody
using ModelingToolkit
using Plots
using OrdinaryDiffEq
using LinearAlgebra

t = Multibody.t
D = Differential(t)
world = Multibody.world

@named joint = Multibody.Revolute(n = [0, 0, 1], isroot = true)
@named body = Body(; m = 1, isroot = false, r_cm = [0.5, 0, 0])
@named torquesensor = CutTorque()
@named forcesensor = CutForce()

connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, body.frame_a, torquesensor.frame_a,
                       forcesensor.frame_a)]

connections = [connect(world.frame_b, joint.frame_a)
               connect(joint.frame_b, torquesensor.frame_a)
               connect(torquesensor.frame_b, forcesensor.frame_a)
               connect(forcesensor.frame_b, body.frame_a)]

@named model = ODESystem(connections, t,
                         systems = [world, joint, body, torquesensor, forcesensor])
modele = ModelingToolkit.expand_connections(model)
ssys = structural_simplify(model, allow_parameter = false)


D = Differential(t)
defs = Dict(collect((D.(joint.phi)) .=> [0, 0, 0])...,
            collect(D.(D.(joint.phi)) .=> [0, 0, 0])...)
prob = ODEProblem(ssys, defs, (0, 10))

using OrdinaryDiffEq
sol = solve(prob, Rodas4())
@assert SciMLBase.successful_retcode(sol)

plot(sol, idxs = [collect(forcesensor.force.u); collect(joint.frame_a.f)])
```

Note how the force sensor measures a force that appears to equal the cut-force in the joint in magnitude, but the orientation appears to differ. Frame cut forces and toques are resolved in the world frame by default, while the force sensor measures the force in the frame of the sensor. We can choose which frame to resolve the measurements in by using hte keyword argument `@named forcesensor = CutForce(; resolveInFrame = :world)`. If we do this, the traces in the plot above will overlap.

Since the torque sensor measures a torque in a revolute joint, it should measure zero torque in this case, no torque is transmitted through the revolute joint since the rotational axis is perpendicular to the gravitational force:
```@example sensor
all(x -> abs(x) < 1e-3, reduce(hcat, sol[torquesensor.torque.u]))
```
