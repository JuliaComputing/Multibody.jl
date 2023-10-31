# Pendulum--The "Hello World of multi-body dynamics"
This beginners tutorial will model a pendulum pivoted around the origin in the world frame. The world frame is a constant that lives inside the Multibody module, all multibody models are "grounded" in the same world, i.e., the `world` component must be included in all models.

![Pendulum](https://doc.modelica.org/Modelica%203.2.3/Resources/Images/Mechanics/MultiBody/Examples/Elementary/Pendulum.png)

To start, we load the required packages
```@example pendulum
using ModelingToolkit
using Multibody, JuliaSimCompiler
using OrdinaryDiffEq # Contains the ODE solver we will use
using Plots
```
We then access the world frame and time variable from the Multibody module
```@example pendulum
t = Multibody.t
world = Multibody.world
show(stdout, MIME"text/plain"(), world)
nothing # hide
```

Unless otherwise specified, the world defaults to have a gravitational field pointing in the negative ``y`` direction and a gravitational acceleration of ``9.81``.


## Modeling the pendulum
Our simple pendulum will initially consist of a [`Body`](@ref) (point mass) and a [`Revolute`](@ref) joint (the pivot joint). We construct these elements by calling their constructors
```@example pendulum
@named joint = Revolute(n = [0, 0, 1], isroot = true)
@named body = Body(; m = 1, isroot = false, r_cm = [0.5, 0, 0])
nothing # hide
```
The `n` argument to [`Revolute`](@ref) denotes the rotational axis of the joint, this vector must have `norm(n) == 1`. We also indicate that the revolute joint is the root of the kinematic tree, i.e., the potential state of the joint will serve as the state variables for the system.

The [`Body`](@ref) is constructed by providing its mass, `m`, and the vector `r_cm` from its mounting frame, `body.frame_a`, to the center of mass. Since the world by default has the gravity field pointing along the negative ``y`` axis, we place the center of mass along the ``x``-axis to make the pendulum swing back and forth. The body is not selected as the root of the kinematic tree, since we have a joint in this system, but if we had attached the body directly to, e.g., a spring, we could set the body to be the root and avoid having to introduce an "artificial joint", which is otherwise needed in order to have at least one component that has a potential state.

To connect the components together, we create a vector of connections using the `connect` function. A joint typically has two frames, `frame_a` and `frame_b`. In this example, the first frame of the joint is attached to the world frame, and the body is attached to the second frame of the joint, i.e., the joint allows the body to swing back and forth. The order of the connections is not important for ModelingToolkit, but it's good practice to follow some convention, here, we start at the world and progress outwards in the kinematic tree.
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
plot(sol, idxs = joint.phi, title="Pendulum")
```
The solution `sol` can be plotted directly if the Plots package is loaded. The figure indicates that the pendulum swings back and forth without any damping. To add damping as well, we could add a `Damper` from the `ModelingToolkitStandardLibrary.Mechanical.Rotational` module to the revolute joint. We do this below

## 3D Animation
Multibody.jl supports automatic 3D rendering of mechanisms, we use this feature to illustrate the result of the simulation below:

```@example pendulum
import CairoMakie # GLMakie is another alternative, suitable for interactive plots
Multibody.render(model, sol; z = -5, filename = "pendulum.gif") # Use "pendulum.mp4" for a video file
nothing # hide
```

![animation](pendulum.gif)

By default, the world frame is indicated using the convention x: red, y: green, z: blue. The animation shows how the simple [`Body`](@ref) represents a point mass at a particular distance `r_cm` away from its mounting flange `frame_a`. To model a more physically motivated pendulum rod, we could have used a [`BodyShape`](@ref) component, which has two mounting flanges instead of one. The [`BodyShape`](@ref) component is shown in several of the examples available in the example sections of the documentation.

## Adding damping
To add damping to the pendulum such that the pendulum will eventually come to rest, we add a [`Damper`](@ref) to the revolute joint. The damping coefficient is given by `d`, and the damping force is proportional to the angular velocity of the joint. To add the damper to the revolute joint, we must create the joint with the keyword argument `useAxisFlange = true`, this adds two internal flanges to the joint to which you can attach components from the `ModelingToolkitStandardLibrary.Mechanical.Rotational` module (1-dimensional components). We then connect one of the flanges of the damper to the axis flange of the joint, and the other damper flange to the support flange which is rigidly attached to the world.
```@example pendulum
@named damper = Rotational.Damper(d = 0.1)
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
plot(sol, idxs = joint.phi, title="Damped pendulum")
```
This time we see that the pendulum loses energy and eventually comes to rest at the stable equilibrium point ``\pi / 2``.

## A linear pendulum?
When we think of a pendulum, we typically think of a rotary pendulum that is rotating around a pivot point like in the examples above. 
A mass suspended in a spring can be though of as a linear pendulum (often referred to as a harmonic oscillator rather than a pendulum), and we show here how we can construct a model of such a device. This time around, we make use of a [`Prismatic`](@ref) joint rather than a [`Revolute`](@ref) joint. A [prismatic joint](https://en.wikipedia.org/wiki/Prismatic_joint) has one positional degree of freedom, compared to the single rotational degree of freedom for the revolute joint.

![Spring with mass](https://doc.modelica.org/Modelica%203.2.3/Resources/Images/Mechanics/MultiBody/Examples/Elementary/SpringWithMass.png)

```@example pendulum
@named damper = Translational.Damper(d=0.5)
@named spring = Translational.Spring(c=1)
@named joint = Prismatic(n = [0, 1, 0], isroot = true, useAxisFlange = true)

connections = [connect(world.frame_b, joint.frame_a)
               connect(damper.flange_b, spring.flange_b, joint.axis)
               connect(joint.support, damper.flange_a, spring.flange_a)
               connect(body.frame_a, joint.frame_b)]

@named model = ODESystem(connections, t, systems = [world, joint, body, damper, spring])
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys, [damper.s_rel => 1, D(D(joint.s)) => 0], (0, 30))

sol = solve(prob, Rodas4())
Plots.plot(sol, idxs = joint.s, title="Mass-spring-damper system")
```

As is hopefully evident from the little code snippet above, this linear pendulum model has a lot in common with the rotary pendulum. In this example, we connected both the spring and a damper to the same axis flange in the joint. This time, the components came from the `Translational` submodule of ModelingToolkitStandardLibrary rather than the `Rotational` submodule. Also here do we pass `useAxisFlange` when we create the joint to make sure that it is equipped with the flanges `support` and `axis` needed to connect the translational components.

### Why do we need a joint?
In the example above, we introduced a prismatic joint to model the oscillating motion of the mass-spring system. In reality, we can suspend a mass in a spring without any joint, so why do we need one here? The answer is that we do not, in fact, need the joint, but if we connect the spring directly to the world, we need to make the body (mass) the root object of the kinematic tree instead:
```@example pendulum
@named root_body = Body(; m = 1, isroot = true, r_cm = [0, 1, 0], phi0 = [0, 1, 0])
@named multibody_spring = Multibody.Spring(c=1)

connections = [connect(world.frame_b, multibody_spring.frame_a)
                connect(root_body.frame_a, multibody_spring.frame_b)]

@named model = ODESystem(connections, t, systems = [world, multibody_spring, root_body])
ssys = structural_simplify(IRSystem(model))

defs = Dict(collect(multibody_spring.r_rel_0 .=> [0, 1, 0])...,
            collect(root_body.r_0 .=> [0, 0, 0])...,
            collect((D.(root_body.phi)) .=> [0, 0, 0])...,
            collect(D.(D.(root_body.phi)) .=> [0, 0, 0])...,
            collect(D.(root_body.phid) .=> [0, 0, 0])...,)

prob = ODEProblem(ssys, defs, (0, 30))

sol = solve(prob, Rodas4())
plot(sol, idxs = multibody_spring.r_rel_0[2], title="Mass-spring system without joint")
```
Here, we used a [`Multibody.Spring`](@ref) instead of connecting a `Translational.Spring` to a joint. The `Translational.Spring`, alongside other components from `ModelingToolkitStandardLibrary.Mechanical`, is a 1-dimensional object, whereas multibody components are 3-dimensional objects.

Internally, the [`Multibody.Spring`](@ref) contains a `Translational.Spring`, attached between two flanges, so we could actually add a damper to the system as well:
```@example pendulum
push!(connections, connect(multibody_spring.spring2d.flange_a, damper.flange_a))
push!(connections, connect(multibody_spring.spring2d.flange_b, damper.flange_b))

@named model = ODESystem(connections, t, systems = [world, multibody_spring, root_body, damper])
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, defs, (0, 30))

sol = solve(prob, Rodas4())
plot(sol, idxs = multibody_spring.r_rel_0[2], title="Mass-spring-damper without joint")
```

The figure above should look identical to the simulation of the mass-spring-damper further above.

## Going 3D
The systems we have modeled so far have all been _planar_ mechanisms. We now extend this to a 3-dimensional system, the [_Furuta pendulum_](https://en.wikipedia.org/wiki/Furuta_pendulum).

This pendulum, sometimes referred to as a _rotary pendulum_, has two joints, one in the "shoulder", which is typically configured to rotate around the gravitational axis, and one in the "elbow", which is typically configured to rotate around the axis of the upper arm. The upper arm is attached to the shoulder joint, and the lower arm is attached to the elbow joint. The tip of the pendulum is attached to the lower arm.

```@example pendulum
using ModelingToolkit, Multibody, JuliaSimCompiler, OrdinaryDiffEq, Plots
import ModelingToolkitStandardLibrary.Mechanical.Rotational.Damper as RDamper
import Multibody.Rotations
W(args...; kwargs...) = Multibody.world

@mtkmodel FurutaPendulum begin
    @components begin
        world = W()
        shoulder_joint = Revolute(n = [0, 1, 0], isroot = true, useAxisFlange = true)
        elbow_joint    = Revolute(n = [0, 0, 1], isroot = true, useAxisFlange = true, phi0=0.1)
        upper_arm = BodyShape(; m = 0.1, isroot = false, r = [0, 0, 0.6], radius=0.05)
        lower_arm = BodyShape(; m = 0.1, isroot = false, r = [0, 0.6, 0], radius=0.05)
        tip = Body(; m = 0.3, isroot = false)

        damper1 = RDamper(d = 0.07)
        damper2 = RDamper(d = 0.07)
    end
    @equations begin
        connect(world.frame_b, shoulder_joint.frame_a)
        connect(shoulder_joint.frame_b, upper_arm.frame_a)
        connect(upper_arm.frame_b, elbow_joint.frame_a)
        connect(elbow_joint.frame_b, lower_arm.frame_a)
        connect(lower_arm.frame_b, tip.frame_a)

        connect(shoulder_joint.axis, damper1.flange_a)
        connect(shoulder_joint.support, damper1.flange_b)

        connect(elbow_joint.axis, damper2.flange_a)
        connect(elbow_joint.support, damper2.flange_b)

    end
end

@named model = FurutaPendulum()
model = complete(model)
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys, [model.shoulder_joint.phi => 0.0, model.elbow_joint.phi => 0.1], (0, 12))
sol = solve(prob, Rodas4())
plot(sol, layout=4)
```

```@example pendulum
import CairoMakie
Multibody.render(model, sol, z=-5, framerate=60, R = Rotations.RotXYZ(0.2, -0.2, 0), filename = "furuta.gif")
nothing # hide
```
![furuta](furuta.gif)


### Orientations and directions
Let's break down how to think about directions and orientations when building 3D mechanisms. In the example above, we started with the shoulder joint, this joint rotated around the gravitational axis, `n = [0, 1, 0]`. When this joint is positioned in joint coordinate `shoulder_joint.phi = 0`, its `frame_a` and `frame_b` will coincide. When the joint rotates, `frame_b` will rotate around the axis `n` of `frame_a`. The `frame_a` of the joint is attached to the world, so the joint will rotate around the world's `y`-axis:

```@example pendulum
function get_R(frame, t)
    reshape(sol(t, idxs=vec(ori(frame).R.mat')), 3, 3)
end
function get_r(frame, t)
    sol(t, idxs=collect(frame.r_0))
end
get_R(model.shoulder_joint.frame_b)
```
we see that at time $t = 0$, we have no rotation of `frame_b` around the $y$ axis of the world (frames are always resolved in the world frame), but a second into the simulation, we do:
```@example pendulum
get_R(model.shoulder_joint.frame_b, 1)
```

The next body is the upper arm. This body has an extent of `0.6` in the $z$ direction, as measured in its local `frame_a`
```@example pendulum
get_r(model.upper_arm.frame_b, 0)
```
One second into the simulation, the upper arm has rotated around the $y$ axis of the world
```@example pendulum
get_r(model.upper_arm.frame_b, 1)
```

If we look at the variable `model.upper_arm.r`, we do not see this rotation!
```@example pendulum
arm_r = sol(1, idxs=collect(model.upper_arm.r))
```
The reason is that this variable is resolved in the local `frame_a` and not in the world frame. To transform this variable to the world frame, we may multiply with the rotation matrix of `frame_a`
```@example pendulum
get_R(model.upper_arm.frame_a, 1)*arm_r
```


Slightly more formally, let $R_A^B$ denote the rotation matrix that rotates a vector expressed in a frame $A$ into one that is expressed in frame $B$, i.e., $r_B = R_B^A r_A$. We have then just performed the transformation $r_W = R_W^A r_A$, where $W$ denotes the world frame, and $A$ denotes `body.frame_a`.

The next joint, the elbow joint, has the rotational axis `n = [0, 0, 1]`. This indicates that the joint rotates around the $z$-axis of its `frame_a`. Since the upper arm was oriented along the $z$ direction, the joint is rotating around the axis that coincides with the upper arm. 

The lower arm is finally having an extent in the $y$-axis. At the final time when the pendulum motion has been full damped, we see that the second frame of this body ends up with an $y$-coordinate of `-0.6`:
```@example pendulum
get_r(model.lower_arm.frame_b, 12)
```

If we rotate the vector of extent of the lower arm to the world frame, we indeed see that the only coordinate that is nonzero is the $y$ coordinate:
```@example pendulum
get_R(model.lower_arm.frame_a, 12)*sol(12, idxs=collect(model.lower_arm.r))
```

The reason that the latter vector differs from `get_r(model.lower_arm.frame_b, 12)` above is that `get_r(model.lower_arm.frame_b, 12)` has been _translated_ as well. To both translate and rotate `model.lower_arm.r` into the world frame, we must use the full transformation matrix $T_W_A \in SE(3)$:

```@example pendulum
function get_T(frame, t)
    R = get_R(frame, t)
    r = get_r(frame, t)
    [R r; 0 0 0 1]
end

r_A = sol(12, idxs=collect(model.lower_arm.r))
r_A = [r_A; 1] # Homogeneous coordinates

get_T(model.lower_arm.frame_a, 12)*r_A
```
the vector is now coinciding with `get_r(model.lower_arm.frame_b, 12)`.