# Gyroscopic effects

![animation](gyro.gif)

In this example, we demonstrate how a rotating body creates a reaction torque when its axis of rotation is changed.

The system consists of a pendulum suspended in a spherical joint, a joint without any rotational constraints. The tip of the pendulum is a cylinder that is rotating around a revolute joint in its center. When the pendulum swings, the rotation axis of the rotating tip is changed, this causes the entire pendulum to rotate around the axis through the pendulum rod.


```@example spring_mass_system
using Multibody
using ModelingToolkit
using Plots
using JuliaSimCompiler
using OrdinaryDiffEq

t = Multibody.t
D = Differential(t)
world = Multibody.world

systems = @named begin
    spherical = Spherical(state=true, radius=0.02, color=[1,1,0,1])
    body1 = BodyCylinder(r = [0.25, 0, 0], diameter = 0.05, isroot=false, quat=false)
    rot = FixedRotation(; n = [0,1,0], angle=deg2rad(45))
    revolute = Revolute(n = [1,0,0], radius=0.06, color=[1,0,0,1])
    trans = FixedTranslation(r = [-0.1, 0, 0])
    body2 = BodyCylinder(r = [0.2, 0, 0], diameter = 0.1, color=[0,0,0.5,1])
end

connections = [
    connect(world.frame_b, spherical.frame_a)
    connect(spherical.frame_b, body1.frame_a)
    connect(body1.frame_b, rot.frame_a)
    connect(rot.frame_b, revolute.frame_a)
    connect(revolute.frame_b, trans.frame_a)
    connect(trans.frame_b, body2.frame_a)
]

@named model = ODESystem(connections, t, systems = [world; systems])
model = complete(model)
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys, [model.world.g => 9.80665, model.revolute.w => 10], (0, 5))

sol = solve(prob, Rodas5P(), abstol=1e-6, reltol=1e-6)
@assert SciMLBase.successful_retcode(sol)
using Test # hide
@test sol(5, idxs=collect(model.body2.r_0[1:3])) ≈ [-0.0357364, -0.188245, 0.02076935] atol=1e-3 # hide
# plot(sol, idxs=collect(model.body2.r_0)) # hide

import CairoMakie
Multibody.render(model, sol; x=1, z=1, filename = "gyro.gif") # Use "gyro.mp4" for a video file
nothing # hide
```

![animation](gyro.gif)

Try setting `model.revolute.w => 0` and plot this variable using `plot(sol, idxs=model.revolute.w)` and you will notice that the swinging of the pendulum induces a rotation around this joint, even if it has no rotational velocity from the start.