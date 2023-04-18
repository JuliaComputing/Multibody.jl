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

## Modeling the pendulum
Our simple pendulum will initially consist of a [`Body`](@ref) and a [`Revolute`](@ref) joint (the pivot joint). We construct these elements by calling their constructors
```@example pendulum
@named joint = Revolute(n = [0, 0, 1], isroot = true)
@named body = Body(; m = 1, isroot = false, r_cm = [0.5, 0, 0])
nothing # hide
```
The `n` argument to [`Revolute`](@ref) denotes the rotational axis of the joint, this vector must have `norm(n) == 1`. We also indicate that the revolute joint is the root of the kinematic tree, i.e., the potential states of the joint will serve as the state variables for the system.

The [`Body`](@ref) is constructed by providing its mass, `m`, and the vector `r_cm` from its first frame, `body.frame_a`, to the center of mass. Since the world by default has the gravity field pointing along the negative ``y`` axis, we place the center of mass along the ``x``-axis to make the pendulum swing back and forth. The body is not selected as the root of the kinematic tree, since we have a joint in this system, but if we had attached the body directly to, e.g., a spring, we could set the body to be the root and avoid having to introduce an "artificial joint".

To connect the components together, we create a vector of connections using the `connect` function. A joint typically has two frames, `frame_a` and `frame_b`. The first frame of the joint is attached to the world frame, and the body is attached to the second joint frame. The order of the connections is not important for ModelingToolkit, but it's good practice to follow some convention, here, we start at the world and progress outwards in the kinematic tree.
```@example pendulum
connections = [
    connect(world.frame_b, joint.frame_a)
    connect(joint.frame_b, body.frame_a)
]
nothing # hide
```

With all components and connections defined, we can create an `ODESystem` like so:
```@example pendulum
@named model = ODESystem(connections, t, systems=[world, joint, body])
nothing # hide
```
The `ODESystem` is the fundamental model type in ModelingToolkit used for multibody-type models.

Before we can simulate the system, we must perform model compilation using [`structural_simplify`](@ref)
```@example pendulum
ssys = structural_simplify(model, allow_parameter = false)
```
This results in a simplified model with the minimum required variables and equations to be able to simulate the system efficiently. This step rewrites all `connect` statements into the appropriate equations, and removes any redundant variables and equations.

We are now ready to create an `ODEProblem` and simulate it. We use the `Rodas4` solver from OrdinaryDiffEq.jl, and pass a dictionary for the initial conditions. We specify only initial condition for some variables, for those variables where no initial condition is specified, the default initial condition defined the model will be used.
```@example pendulum
D = Differential(t)
defs = Dict(D(joint.phi) => 0, D(D(joint.phi)) => 0)
prob = ODEProblem(ssys, defs, (0, 10))

sol = solve(prob, Rodas4())
plot(sol, idxs = joint.phi)
```
The solution `sol` can be plotted directly if the Plots package is loaded. The figure indicates that the pendulum swings back and forth without any damping. To add damping as well, we could add a `Damper` from the `ModelingToolkitStandardLibrary.Mechanical.Rotational` module to the revolute joint. We do this below

## Adding damping
To add damping to the pendulum such that the pendulum will eventually come to rest, we add a [`Damper`](@ref) to the revolute joint. The damping coefficient is given by `d`, and the damping force is proportional to the angular velocity of the joint. To add the damper to the revolute joint, we must create the joint with the keyword argument `useAxisFlange = true`, this adds two internal flanges to the joint to which you can attach components from the `ModelingToolkitStandardLibrary.Mechanical.Rotational` module. We then connect one of the flanges of the damper to the axis flange of the joint, and the other damper flange to the support flange which is rigidly attached to the world.
```@example pendulum
using ModelingToolkitStandardLibrary.Mechanical.Rotational

@named damper = Damper(d = 0.1)
@named joint = Revolute(n = [0, 0, 1], isroot = true, useAxisFlange = true)

connections = [connect(world.frame_b, joint.frame_a)
               connect(damper.flange_b, joint.axis)
               connect(joint.support, damper.flange_a)
               connect(body.frame_a, joint.frame_b)]

@named model = ODESystem(connections, t, systems = [world, joint, body, damper])
ssys = structural_simplify(model, allow_parameter = false)


prob = ODEProblem(ssys, [damper.phi_rel => 1, D(joint.phi) => 0, D(D(joint.phi)) => 0],
                  (0, 30))


sol = solve(prob, Rodas4())
plot(sol, idxs = collect(joint.phi))
```
This time we see that the pendulum loses energy and eventually comes to rest at the stable equilibrium point ``\pi / 2``.