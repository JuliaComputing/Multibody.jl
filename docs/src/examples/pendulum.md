# Pendulum---The "Hello World of multi-body dynamics"
This beginners tutorial will model a pendulum pivoted around the origin in the world frame. The world frame is a constant that lives inside the Multibody module, all multibody models are "grounded" in the same world. To start, we load the required packages
```@example pendulum
using ModelingToolkit
using Multibody
using OrdinaryDiffEq # Contains the ODE solver we will use
using Plots
```
We then access the world frame and time variable from the Multibody module
```@example pendulum
t = Multibody.t
world = Multibody.world
```
Our simple pendulum consists of a [`Body`](@ref) and a [`Revolute`](@ref) joint (the pivot joint). We construct these elements by calling their constructors
```
@named body = Body(; m=1, isroot=false)
@named joint = Multibody.Revolute(; ϕ0=1)
```
To connect the components together, we create a vector of connections using the `connect` function. A joint typically has two frames, `frame_a` and `frame_b`. The first frame of the joint is attached to the world frame, and the body is attached to the second joint frame. The order of the connections is not important for ModelingToolkit, but it's good practice to follow some convention, here, we start at the world and progress outwards in the kinematic tree.
```@example pendulum
connections = [
    connect(world.frame_b, joint.frame_a)
    connect(joint.frame_b, body.frame_a)
]
```

With all components and connections defined, we can create an `ODESystem` like so:
```@example pendulum
@named model = ODESystem(connections, t, systems=[world, joint, body])
```
The `ODESystem` is the fundamental model type in ModelingToolkit used for multibody-type models.

Before we can simulate the system, we must perform model compilation using [`structural_simplify`](@ref)
```@example pendulum
sys = structural_simplify(model)
```
This results in a simplified model with the minimum required variables and equations to be able to simulate the system efficiently. This step rewrites all `connect` statements into the appropriate equations, and removes any redundant variables and equations.

We are now ready to create an `ODEProblem` and simulate it. We use the `Tsit5` solver from OrdinaryDiffEq.jl, and pass an empty array `[]` for the initial condition, this means that we use the default intial condition defined in the model.
```@example pendulum
prob = ODEProblem(sys, [], (0.0, 10.0))
sol = solve(prob, Tsit5())
plot(sol)
```
The solution `sol` can be plotted directly if the Plots package is loaded. The figure indicates that the pendulum swings back and forth without any damping. To add damping as well, we could add a `Damper` from the `ModelingToolkitStandardLibrary.Mechanical.Rotational` module to the revolute joint.