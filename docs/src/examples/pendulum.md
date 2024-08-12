# Pendulum--The "Hello World of multi-body dynamics"
This beginners tutorial will start by modeling a pendulum pivoted around the origin in the world frame. The world frame is a constant that lives inside the Multibody module, all multibody models are "grounded" in the same world, i.e., the `world` component must be included in all models.

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
show(stdout, MIME"text/plain"(), world) # hide
nothing # hide
```

Unless otherwise specified, the world defaults to having a gravitational field pointing in the negative ``y`` direction and a gravitational acceleration of ``9.81`` (See [Bodies in space](@ref) for more options).


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
model = complete(model)
nothing # hide
```
The `ODESystem` is the fundamental model type in ModelingToolkit used for multibody-type models.

Before we can simulate the system, we must perform model compilation using [`structural_simplify`](@ref)
```@example pendulum
ssys = structural_simplify(IRSystem(model))
```
This results in a simplified model with the minimum required variables and equations to be able to simulate the system efficiently. This step rewrites all `connect` statements into the appropriate equations, and removes any redundant variables and equations. To simulate the pendulum, we require two state variables, one for angle and one for angular velocity, we can see above that these state variables have indeed been chosen.

We are now ready to create an `ODEProblem` and simulate it. We use the `Rodas4` solver from OrdinaryDiffEq.jl, and pass a dictionary for the initial conditions. We specify only initial condition for some variables, for those variables where no initial condition is specified, the default initial condition defined the model will be used.
```@example pendulum
D = Differential(t)
defs = Dict() # We may specify the initial condition here
prob = ODEProblem(ssys, defs, (0, 3.35))

sol = solve(prob, Rodas4())
plot(sol, idxs = joint.phi, title="Pendulum")
```
The solution `sol` can be plotted directly if the Plots package is loaded. The figure indicates that the pendulum swings back and forth without any damping. To add damping as well, we could add a `Damper` from the `ModelingToolkitStandardLibrary.Mechanical.Rotational` module to the revolute joint. We do this below

## 3D Animation
Multibody.jl supports automatic 3D rendering of mechanisms, we use this feature to illustrate the result of the simulation below:

```@example pendulum
import GLMakie # GLMakie is another alternative, suitable for interactive plots
Multibody.render(model, sol; filename = "pendulum.gif") # Use "pendulum.mp4" for a video file
nothing # hide
```

![animation](pendulum.gif)

By default, the world frame is indicated using the convention x: red, y: green, z: blue. The animation shows how the simple [`Body`](@ref) represents a point mass with inertial properties at a particular distance `r_cm` away from its mounting flange `frame_a`. The cylinder that is shown connecting the pivot point to the body is for visualization purposes only, it does not have any inertial properties. To model a more physically motivated pendulum rod, we could have used a [`BodyShape`](@ref) component, which has two mounting flanges instead of one. The [`BodyShape`](@ref) component is shown in several of the examples available in the example sections of the documentation.

## Adding damping
To add damping to the pendulum such that the pendulum will eventually come to rest, we add a [`Damper`](@ref) to the revolute joint. The damping coefficient is given by `d`, and the damping force is proportional to the angular velocity of the joint. To add the damper to the revolute joint, we must create the joint with the keyword argument `axisflange = true`, this adds two internal flanges to the joint to which you can attach components from the `ModelingToolkitStandardLibrary.Mechanical.Rotational` module (1-dimensional components). We then connect one of the flanges of the damper to the axis flange of the joint, and the other damper flange to the support flange which is rigidly attached to the world.
```@example pendulum
@named damper = Rotational.Damper(d = 0.3)
@named joint = Revolute(n = [0, 0, 1], isroot = true, axisflange = true)

connections = [connect(world.frame_b, joint.frame_a)
               connect(damper.flange_b, joint.axis)
               connect(joint.support, damper.flange_a)
               connect(body.frame_a, joint.frame_b)]

@named model = ODESystem(connections, t, systems = [world, joint, body, damper])
model = complete(model)
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys, [damper.phi_rel => 1], (0, 10))

sol = solve(prob, Rodas4())
plot(sol, idxs = joint.phi, title="Damped pendulum")
```
This time we see that the pendulum loses energy and eventually comes to rest at the stable equilibrium point ``\pi / 2``.

```@example pendulum
Multibody.render(model, sol; filename = "pendulum_damped.gif")
nothing # hide
```
![animation damped](pendulum_damped.gif)

## A linear pendulum?
When we think of a pendulum, we typically think of a rotary pendulum that is rotating around a pivot point like in the examples above. 
A mass suspended in a spring can be though of as a linear pendulum (often referred to as a harmonic oscillator rather than a pendulum), and we show here how we can construct a model of such a device. This time around, we make use of a [`Prismatic`](@ref) joint rather than a [`Revolute`](@ref) joint. A [prismatic joint](https://en.wikipedia.org/wiki/Prismatic_joint) has one positional degree of freedom, compared to the single rotational degree of freedom for the revolute joint.


```@example pendulum
@named body_0 = Body(; m = 1, isroot = false, r_cm = [0, 0, 0])
@named damper = Translational.Damper(d=1)
@named spring = Translational.Spring(c=10)
@named joint = Prismatic(n = [0, 1, 0], axisflange = true)

connections = [connect(world.frame_b, joint.frame_a)
               connect(damper.flange_b, spring.flange_b, joint.axis)
               connect(joint.support, damper.flange_a, spring.flange_a)
               connect(body_0.frame_a, joint.frame_b)]

@named model = ODESystem(connections, t, systems = [world, joint, body_0, damper, spring])
model = complete(model)
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys, [], (0, 10))

sol = solve(prob, Rodas4())
Plots.plot(sol, idxs = joint.s, title="Mass-spring-damper system")
```

As is hopefully evident from the little code snippet above, this linear pendulum model has a lot in common with the rotary pendulum. In this example, we connected both the spring and a damper to the same axis flange in the joint. This time, the components came from the `Translational` submodule of ModelingToolkitStandardLibrary rather than the `Rotational` submodule. Also here do we pass `axisflange` when we create the joint to make sure that it is equipped with the flanges `support` and `axis` needed to connect the translational components.

```@example pendulum
Multibody.render(model, sol; filename = "linear_pend.gif", framerate=24)
nothing # hide
```
![linear pendulum](linear_pend.gif)

### Why do we need a joint?
In the example above, we introduced a prismatic joint to model the oscillating motion of the mass-spring system. In reality, we can suspend a mass in a spring without any joint, so why do we need one here? The answer is that we do not, in fact, need the joint, but if we connect the spring directly to the world, we need to make the body (mass) the root object of the kinematic tree instead:
```@example pendulum
@named root_body = Body(; m = 1, isroot = true, r_cm = [0, 1, 0], phi0 = [0, 1, 0])
@named multibody_spring = Multibody.Spring(c=10)

connections = [connect(world.frame_b, multibody_spring.frame_a)
                connect(root_body.frame_a, multibody_spring.frame_b)]

@named model = ODESystem(connections, t, systems = [world, multibody_spring, root_body])
model = complete(model)
ssys = structural_simplify(IRSystem(model))

defs = Dict(collect(root_body.r_0) .=> [0, 1e-3, 0]) # The spring has a singularity at zero length, so we start some distance away

prob = ODEProblem(ssys, defs, (0, 10))

sol = solve(prob, Rodas4())
plot(sol, idxs = multibody_spring.r_rel_0[2], title="Mass-spring system without joint")
```
Here, we used a [`Multibody.Spring`](@ref) instead of connecting a `Translational.Spring` to a joint. The `Translational.Spring`, alongside other components from `ModelingToolkitStandardLibrary.Mechanical`, is a 1-dimensional object, whereas multibody components are 3-dimensional objects.

Internally, the [`Multibody.Spring`](@ref) contains a `Translational.Spring`, attached between two flanges, so we could actually add a damper to the system as well:
```@example pendulum
push!(connections, connect(multibody_spring.spring2d.flange_a, damper.flange_a))
push!(connections, connect(multibody_spring.spring2d.flange_b, damper.flange_b))

@named model = ODESystem(connections, t, systems = [world, multibody_spring, root_body, damper])
model = complete(model)
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, defs, (0, 10))

sol = solve(prob, Rodas4(), u0 = prob.u0 .+ 1e-5 .* randn.())
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
        shoulder_joint = Revolute(n = [0, 1, 0], axisflange = true)
        elbow_joint    = Revolute(n = [0, 0, 1], axisflange = true, phi0=0.1)
        upper_arm = BodyShape(; m = 0.1, r = [0, 0, 0.6], radius=0.04)
        lower_arm = BodyShape(; m = 0.1, r = [0, 0.6, 0], radius=0.04)
        tip = Body(; m = 0.3)

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

prob = ODEProblem(ssys, [model.shoulder_joint.phi => 0.0, model.elbow_joint.phi => 0.1], (0, 10))
sol = solve(prob, Rodas4())
plot(sol, layout=4)
```

In the animation below, we visualize the path that the origin of the pendulum tip traces by providing the tip frame in a vector of frames passed to `traces`
```@example pendulum
import GLMakie
Multibody.render(model, sol, filename = "furuta.gif", traces=[model.tip.frame_a])
nothing # hide
```
![furuta](furuta.gif)


## Orientations and directions
Let's break down how to think about directions and orientations when building 3D mechanisms. In the example above, we started with the shoulder joint, this joint rotated around the gravitational axis, `n = [0, 1, 0]`. When this joint is positioned in joint coordinate `shoulder_joint.phi = 0`, its `frame_a` and `frame_b` will coincide. When the joint rotates, `frame_b` will rotate around the axis `n` of `frame_a`. The `frame_a` of the joint is attached to the world, so the joint will rotate around the world's `y`-axis:

```@example pendulum
get_rot(sol, model.shoulder_joint.frame_b, 0)
```
we see that at time $t = 0$, we have no rotation of `frame_b` around the $y$ axis of the world (frames are always resolved in the world frame), but a second into the simulation, we have:
```@example pendulum
R1 = get_rot(sol, model.shoulder_joint.frame_b, 1)
```
Here, the `frame_b` has rotated around the $y$ axis of the world (if you are not familiar with rotation matrices, we can ask for the rotation axis and angle)
```@example pendulum
using Multibody.Rotations
rotation_axis(R1), rotation_angle(R1)
```

This rotation axis and angle should correspond to the joint coordinate (the orientation described by an axis and an angle is invariant to a multiplication of both by -1)
```@example pendulum
sol(1, idxs=model.shoulder_joint.phi)
```

!!! note "Convention"
    The convention used in [`get_rot`](@ref) is to return the rotation matrix ``R_{World}^{Local}`` that rotates a coordinate from the local frame to the world frame, ``r_{World} = R_{World}^{Local} r_{Local}``.



Here, we made use of the function [`get_rot`](@ref), we will now make use of also [`get_trans`](@ref) and [`get_frame`](@ref).

The next body is the upper arm. This body has an extent of `0.6` in the $z$ direction, as measured in its local `frame_a`
```@example pendulum
get_trans(sol, model.upper_arm.frame_b, 0)
```
One second into the simulation, the upper arm has rotated around the $y$ axis of the world
```@example pendulum
rb1 = get_trans(sol, model.upper_arm.frame_b, 1)
```

If we look at the variable `model.upper_arm.r`, we do not see this rotation!
```@example pendulum
arm_r = sol(1, idxs=collect(model.upper_arm.r))
```
The reason is that this variable is resolved in the local `frame_a` and not in the world frame. To transform this variable to the world frame, we may multiply with the rotation matrix of `frame_a` which is always resolved in the world frame:
```@example pendulum
get_rot(sol, model.upper_arm.frame_a, 1)*arm_r
```
We now get the same result has when we asked for the translation vector of `frame_b` above.
```@example pendulum
using Test # hide
get_rot(sol, model.upper_arm.frame_a, 1)*arm_r ≈ rb1 # hide
nothing # hide
```


Slightly more formally, let $R_A^B$ denote the rotation matrix that rotates a vector expressed in a frame $A$ into one that is expressed in frame $B$, i.e., $r_B = R_B^A r_A$. We have then just performed the transformation $r_W = R_W^A r_A$, where $W$ denotes the world frame, and $A$ denotes `body.frame_a`.

The next joint, the elbow joint, has the rotational axis `n = [0, 0, 1]`. This indicates that the joint rotates around the $z$-axis of its `frame_a`. Since the upper arm was oriented along the $z$ direction, the joint is rotating around the axis that coincides with the upper arm. 

The lower arm is finally having an extent along the $y$-axis. At the final time when the pendulum motion has been fully damped, we see that the second frame of this body ends up with an $y$-coordinate of `-0.6`:
```@example pendulum
get_trans(sol, model.lower_arm.frame_b, 12)
```

If we rotate the vector of extent of the lower arm to the world frame, we indeed see that the only coordinate that is nonzero is the $y$ coordinate:
```@example pendulum
get_rot(sol, model.lower_arm.frame_a, 12)*sol(12, idxs=collect(model.lower_arm.r))
```

The reason that the latter vector differs from `get_trans(sol, model.lower_arm.frame_b, 12)` above is that `get_trans(sol, model.lower_arm.frame_b, 12)` has been _translated_ as well. To both translate and rotate `model.lower_arm.r` into the world frame, we must use the full transformation matrix $T_W^A \in SE(3)$:

```@example pendulum
r_A = sol(12, idxs=collect(model.lower_arm.r))
r_A = [r_A; 1] # Homogeneous coordinates

get_frame(sol, model.lower_arm.frame_a, 12)*r_A
```
the vector is now coinciding with `get_trans(sol, model.lower_arm.frame_b, 12)`.


## Control-design example: Pendulum on cart
We will now demonstrate a complete workflow including
- Modeling
- Linearizaiton
- Control design

We will continue the pendulum theme and design an inverted pendulum on cart. The cart is modeled as [`BodyShape`](@ref) with specified mass, and `shape = "box"` to render it as a box in animations. The cart is moving along the $x$-axis by means of a [`Prismatic`](@ref) joint, and the pendulum is attached to the cart by means of a [`Revolute`](@ref) joint. The pendulum is a [`BodyCylinder`](@ref) with a diameter of `0.015`, the mass and inertia properties are automatically computed using the geometrical dimensions and the density (which defaults to that of steel). A force is applied to the cart by means of a `TranslationalModelica.Force` component.

We start by putting the model together and control it in open loop using a simple periodic input signal:

```@example pendulum
import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica
import ModelingToolkitStandardLibrary.Blocks
using Plots
W(args...; kwargs...) = Multibody.world
gray = [0.5, 0.5, 0.5, 1]
@mtkmodel Cartpole begin
    @components begin
        world = W()
        cart = BodyShape(m = 1, r = [0.2, 0, 0], color=[0.2, 0.2, 0.2, 1], shape="box")
        mounting_point = FixedTranslation(r = [0.1, 0, 0])
        prismatic = Prismatic(n = [1, 0, 0], axisflange = true, color=gray, state_priority=100)
        revolute = Revolute(n = [0, 0, 1], axisflange = false, state_priority=100)
        pendulum = BodyCylinder(r = [0, 0.5, 0], diameter = 0.015, color=gray)
        motor = TranslationalModelica.Force(use_support = true)
        tip = Body(m = 0.05)
    end
    @variables begin
        u(t) = 0
        x(t)
        v(t)
        phi(t)
        w(t)
    end
    @equations begin
        connect(world.frame_b, prismatic.frame_a)
        connect(prismatic.frame_b, cart.frame_a, mounting_point.frame_a)
        connect(mounting_point.frame_b, revolute.frame_a)
        connect(revolute.frame_b, pendulum.frame_a)
        connect(pendulum.frame_b, tip.frame_a)
        connect(motor.flange, prismatic.axis)
        connect(prismatic.support, motor.support)
        u ~ motor.f.u
        x ~ prismatic.s
        v ~ prismatic.v
        phi ~ revolute.phi
        w ~ revolute.w
    end
end
@mtkmodel CartWithInput begin
    @components begin
        cartpole = Cartpole()
        input = Blocks.Cosine(frequency=1, amplitude=1)
    end
    @equations begin
        connect(input.output, :u, cartpole.motor.f)
    end
end
@named model = CartWithInput()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, [model.cartpole.prismatic.s => 0.0, model.cartpole.revolute.phi => 0.1], (0, 10))
sol = solve(prob, Tsit5())
plot(sol, layout=4)
```
As usual, we render the simulation in 3D to get a better feel for the system:
```@example pendulum
import GLMakie
Multibody.render(model, sol, filename = "cartpole.gif", traces=[model.cart.pendulum.frame_b])
nothing # hide
```
![cartpole](cartpole.gif)

### Adding feedback

We will attempt to stabilize the pendulum in the upright position by using feedback control. To design the contorller, we linearize the model in the upward equilibrium position and design an infinite-horizon LQR controller using ControlSystems.jl. We then connect the controller to the motor on the cart. See also [RobustAndOptimalControl.jl: Control design for a pendulum on a cart](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/cartpole/) for a similar example with more detail on the control design.

### Linearization
We start by linearizing the model in the upward equilibrium position using the function `ModelingToolkit.linearize`.
```@example pendulum
import ModelingToolkit: D_nounits as D
using LinearAlgebra
@named cp = Cartpole()
cp = complete(cp)
inputs = [cp.u] # Input to the linearized system
outputs = [cp.x, cp.phi, cp.v, cp.w] # These are the outputs of the linearized system
op = Dict([ # Operating point to linearize in
    cp.u => 0
    cp.revolute.phi => 0 # Pendulum pointing upwards
]
)
matrices, simplified_sys = linearize(IRSystem(cp), inputs, outputs; op)
matrices
```
This gives us the matrices $A,B,C,D$ in a linearized statespace representation of the system. To make these easier to work with, we load the control packages and call `named_ss` instead of `linearize` to get a named statespace object instead:
```@example pendulum
using ControlSystemsMTK
lsys = named_ss(IRSystem(cp), inputs, outputs; op) # identical to linearize, but packages the resulting matrices in a named statespace object for convenience
```

### LQR Control design
With a linear statespace object in hand, we can proceed to design an LQR controller. Since the function `lqr` operates on the state vector, and we have access to the specified output vector, we make use of the system ``C`` matrix to reformulate the problem in terms of the outputs. This relies on the ``C`` matrix being full rank, which is the case here since our outputs include a complete state realization of the system.

To make the simulation interesting, we make a change in the reference position of the cart after a few seconds. 
```@example pendulum
using ControlSystemsBase
C = lsys.C
Q = Diagonal([1, 1, 1, 1])
R = Diagonal([0.1])
Lmat = lqr(lsys, C'Q*C, R)/C # Compute LQR feedback gain. The multiplication by the C matrix is to handle the difference between state and output

@mtkmodel CartWithFeedback begin
    @components begin
        cartpole = Cartpole()
        L = Blocks.MatrixGain(K = Lmat)
        reference = Blocks.Step(start_time = 5, height=0.5)
        control_saturation = Blocks.Limiter(y_max = 10) # To limit the control signal magnitude
    end
    begin
        namespaced_outputs = ModelingToolkit.renamespace.(:cartpole, outputs) # Give outputs correct namespace, they are variables in the cartpole system
    end
    @equations begin
        L.input.u[1] ~ reference.output.u - namespaced_outputs[1] # reference cart position - cartpole.x
        L.input.u[2] ~ 0 - namespaced_outputs[2] # cartpole.phi
        L.input.u[3] ~ 0 - namespaced_outputs[3] # cartpole.v
        L.input.u[4] ~ 0 - namespaced_outputs[4] # cartpole.w
        connect(L.output, control_saturation.input)
        connect(control_saturation.output, cartpole.motor.f)
    end
end
@named model = CartWithFeedback()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, [model.cartpole.prismatic.s => 0.1, model.cartpole.revolute.phi => 0.35], (0, 10))
sol = solve(prob, Tsit5())
cp = model.cartpole
plot(sol, idxs=[cp.prismatic.s, cp.revolute.phi, cp.motor.f.u], layout=3)
plot!(sol, idxs=model.reference.output.u, sp=1, l=(:black, :dash), legend=:bottomright)
```

```@example pendulum
Multibody.render(model, sol, filename = "inverted_cartpole.gif", x=1, z=1)
nothing # hide
```
![inverted cartpole](inverted_cartpole.gif)


### Swing up
Below, we add also an energy-based swing-up controller. For more details this kind of swing-up controller, see [Part 7: Control of rotary pendulum using Julia: Swing up control (YouTube)](https://www.youtube.com/watch?v=RhF2NMCYoiw)
```@example pendulum
"Compute total energy, kinetic + potential, for a body rotating around the z-axis of the world"
function energy(body, w)
    g = world.g
    m = body.m
    d2 = body.r_cm[1]^2 + body.r_cm[2]^2 # Squared distance from 
    I = body.I_33 + m*d2 # Parallel axis theorem
    r_cm_worldframe = Multibody.resolve1(ori(body.frame_a), body.r_cm)[2] # Rotate the distance from frame_a to the center of mass to the world frame
    1/2*I*w^2 + 2m*g*(body.frame_a.r_0[2] + r_cm_worldframe) # Assuming rotation around the z-axis
end

normalize_angle(x::Number) = mod(x+3.1415, 2pi)-3.1415

@mtkmodel CartWithSwingup begin
    @components begin
        cartpole = Cartpole()
        L = Blocks.MatrixGain(K = Lmat)
        control_saturation = Blocks.Limiter(y_max = 12) # To limit the control signal magnitude
    end
    @variables begin
        phi(t)
        w(t)
        E(t), [description = "Total energy of the pendulum"]
        u_swing(t), [description = "Swing-up control signal"]
        switching_condition(t)::Bool, [description = "Switching condition that indicates when stabilizing controller is active"]
    end
    @parameters begin
        Er = 3.825676486352941 # Total energy of the cartpole at the top equilibrium position
    end
    begin
        namespaced_outputs = ModelingToolkit.renamespace.(:cartpole, outputs) # Give outputs correct namespace, they are variables in the cartpole system
    end
    @equations begin
        phi ~ normalize_angle(cartpole.phi)
        w ~ cartpole.w
        E ~ energy(cartpole.pendulum.body, w) + energy(cartpole.tip, w)
        u_swing ~ 100*(E - Er)*sign(w*cos(phi-3.1415))

        L.input.u[1] ~ 0 - namespaced_outputs[1] # - cartpole.x
        L.input.u[2] ~ 0 - phi # cartpole.phi but normalized
        L.input.u[3] ~ 0 - namespaced_outputs[3] # cartpole.v
        L.input.u[4] ~ 0 - namespaced_outputs[4] # cartpole.w
        switching_condition ~ abs(phi) < 0.4
        control_saturation.input.u ~ ifelse(switching_condition, L.output.u, u_swing)
        connect(control_saturation.output, cartpole.motor.f)
    end
end
@named model = CartWithSwingup()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
cp = model.cartpole
prob = ODEProblem(ssys, [cp.prismatic.s => 0.0, cp.revolute.phi => 0.99pi], (0, 5))
sol = solve(prob, Tsit5(), dt = 1e-2, adaptive=false)
plot(sol, idxs=[cp.prismatic.s, cp.revolute.phi, cp.motor.f.u, model.E], layout=4)
hline!([0, 2pi], sp=2, l=(:black, :dash), primary=false)
plot!(sol, idxs=[model.switching_condition], sp=2)
```


```@example pendulum
Multibody.render(model, sol, filename = "swingup.gif", x=2, z=2)
nothing # hide
```
![inverted cartpole](swingup.gif)
