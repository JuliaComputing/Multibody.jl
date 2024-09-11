# Wheels

When modeling wheels, there are several different assumptions that can be made, such as
- Is the wheel able to leave the ground or not?
- Can the wheel slip or does it roll perfectly?

The wheel-related components available are
- [`RollingWheel`](@ref): a wheel that can roll on the ground. It cannot slip and it cannot leave the ground.
- [`RollingWheelJoint`](@ref): a lower-level component used in `RollingWheel` to model the kinematics of the wheel, without inertial or mass properties.
- [`RollingWheelSet`](@ref): a set of two wheels connected by an axis. One of the wheels cannot slip, while the other one slips as required to allow the wheel set to turn (no differential is modeled). No wheel can leave the ground.
- [`RollingWheelSetJoint`](@ref): A lower-level component used in `RollingWheelSet` to model the kinematics of the wheel set, without inertial or mass properties.
- [`RollingConstraintVerticalWheel`](@ref): A low-level constraint that is used to enforce a perfectly rolling wheel that is always vertical, i.e., it can only roll forward and not fall down.
- [`PlanarMechanics.SimpleWheel`](@ref): A 2D wheel component with a simple, linear lateral slip model.
- [`PlanarMechanics.SlipBasedWheelJoint`](@ref): A more advanced 2D wheel component with slip-dependent friction characteristics.

All wheel components are limited to rolling on the ``xz`` plane, i.e., the gravity direction must be the default `[0, -1, 0]`.

## Rolling wheel
```@example WHEEL
using Multibody
using ModelingToolkit
import ModelingToolkitStandardLibrary.Mechanical.Rotational
import ModelingToolkitStandardLibrary.Blocks
using Plots
using OrdinaryDiffEq
using LinearAlgebra
using JuliaSimCompiler
using Test

t = Multibody.t
D = Differential(t)
W(args...; kwargs...) = Multibody.world

@mtkmodel WheelInWorld begin
    @components begin
        world = W()
        wheel = RollingWheel(
            radius = 0.3,
            m = 2,
            I_axis = 0.06,
            I_long = 0.12,
            x0 = 0.2,
            z0 = 0.2,
            der_angles = [0, 5, 1],
        )
    end
end

@named worldwheel = WheelInWorld()
worldwheel = complete(worldwheel)

defs = Dict([
    worldwheel.wheel.body.r_0[1] => 0.2;
    worldwheel.wheel.body.r_0[2] => 0.3;
    worldwheel.wheel.body.r_0[3] => 0.2;
])

ssys = structural_simplify(IRSystem(worldwheel))
prob = ODEProblem(ssys, defs, (0, 4))
sol = solve(prob, Tsit5())
@test SciMLBase.successful_retcode(sol)
```

```@example WHEEL
import GLMakie
Multibody.render(worldwheel, sol; filename = "worldwheel.gif")
nothing # hide
```

![wheel animation](worldwheel.gif)

## Wheel set
A [`RollingWheelSet`](@ref) is comprised out of two wheels mounted on a common axis through their axis of rotation. 
```@example WHEEL
@mtkmodel DrivingWheelSet begin
    @components begin
        sine1 = Blocks.Sine(frequency=1, amplitude=2)
        sine2 = Blocks.Sine(frequency=1, amplitude=2, phase=pi/2)
        torque1 = Rotational.Torque()
        torque2 = Rotational.Torque()
        wheels = RollingWheelSet(radius=0.1, m_wheel=0.5, I_axis=0.01, I_long=0.02, track=0.5, state_priority=100)
        bar = FixedTranslation(r = [0.2, 0, 0])
        body = Body(m=0.01, state_priority=1)
        world = W()
    end
    @equations begin
        connect(sine1.output, torque1.tau)
        connect(sine2.output, torque2.tau)
        connect(torque1.flange, wheels.axis1)
        connect(torque2.flange, wheels.axis2)
        connect(wheels.frame_middle, bar.frame_a)
        connect(bar.frame_b, body.frame_a)
    end
end

@named model = DrivingWheelSet()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
# display(unknowns(ssys))
prob = ODEProblem(ssys, [
    model.wheels.wheelSetJoint.prismatic1.s => 0.1
    model.wheels.wheelSetJoint.prismatic2.s => 0.1
], (0, 3))
sol = solve(prob, Tsit5())
@test SciMLBase.successful_retcode(sol)
```

```@example WHEEL
import GLMakie
Multibody.render(model, sol; filename = "wheelset.gif")
nothing # hide
```

![wheelset animation](wheelset.gif)

The [`RollingWheelSet`](@ref) includes constraints that prevent the wheels from leaving the ground and the connector `frame_middle` from rotating around the wheel axis. This means that if two wheel sets are connected to the same body, the system will be over constrained. To solve this, pass `iscut = true` to one of the wheel sets, like below:

```@example WHEEL
wheel_mass = 15
I_axis=0.01
I_long=0.02
wheel_d = 2
wheel_radius = 0.25
tire_black = [0.1, 0.1, 0.1, 1]

@mtkmodel Car begin
    @structural_parameters begin
        l=4
    end
    @parameters begin
        m=108
        I=10
        g=0
    end
    @components begin
        world = W()

        sine1 = Blocks.Sine(frequency=1, amplitude=150)
        sine2 = Blocks.Sine(frequency=1, amplitude=150, phase=pi/6)
        torque1 = Rotational.Torque()
        torque2 = Rotational.Torque()
        front_wheels = RollingWheelSet(radius=wheel_radius, m_wheel=wheel_mass, I_axis, I_long, track=wheel_d)
        rear_wheels = RollingWheelSet(radius=wheel_radius, m_wheel=wheel_mass, I_axis, I_long, track=wheel_d, iscut=true)

        steering_joint = Revolute(n = [0,1,0], axisflange=true, state_priority=100)
        prefer_straight_ahead = Rotational.SpringDamper(d=10, c=10)

        body = BodyShape(;m, r = [l, 0, 0], I_22 = I, radius=0.3)
    end
    @equations begin
        connect(sine1.output, torque1.tau)
        connect(sine2.output, torque2.tau)
        connect(torque1.flange, front_wheels.axis1)
        connect(torque2.flange, front_wheels.axis2)
        connect(front_wheels.frame_middle, steering_joint.frame_a)

        connect(steering_joint.frame_b, body.frame_a)        
        connect(rear_wheels.frame_middle, body.frame_b)

        connect(prefer_straight_ahead.flange_a, steering_joint.axis)
        connect(prefer_straight_ahead.flange_b, steering_joint.support)
    end
end
@named model = Car()
model = complete(model)
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys, [], (0, 6))
sol = solve(prob, Tsit5())

render(model, sol, framerate=30, filename="car.gif", x=6, z=6, y=5)
nothing # hide
```

![car animation](car.gif)

## Simple planar wheel
This example demonstrates how we can model a simple single-track vehicle with planar (2D or 3DOF) components. 

We will use the component [`PlanarMechanics.SimpleWheel`](@ref), together with a [`PlanarMechanics.Revolute`](@ref) joint to connect the front wheel to the [`PlanarMechanics.BodyShape`](@ref) representing the vehicle body. The revolute joint is used for steering.

```@example WHEEL
import Multibody.PlanarMechanics as Pl

@mtkmodel TestWheel begin
    @components begin
        body = Pl.BodyShape(r = [1.0, 0.0], m=1, I=0.1, gy=0)
        revolute = Pl.Revolute()
        wheel1 = Pl.SimpleWheel(color=tire_black)
        wheel2 = Pl.SimpleWheel(color=tire_black, μ=.5)
        thrust_input1 = Blocks.Constant(k=1)
        thrust_input2 = Blocks.Constant(k=0)
    end
    @equations begin
        connect(body.frame_a, revolute.frame_a)
        connect(revolute.frame_b, wheel1.frame_a)
        connect(thrust_input1.output, wheel1.thrust)
        connect(thrust_input2.output, wheel2.thrust)
        revolute.phi ~ deg2rad(50)*sin(2pi*0.2*t)

        connect(wheel2.frame_a, body.frame_b)
    end
end
@named model = TestWheel()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
defs = Dict(unknowns(ssys) .=> 0)
prob = ODEProblem(ssys, defs, (0.0, 10.0))
sol = solve(prob, Rodas5P())
@test SciMLBase.successful_retcode(sol)
render(model, sol, show_axis=true, x=1, y=-1.8, z=5, lookat=[1,-1.8,0], traces=[model.wheel1.frame_a, model.wheel2.frame_a], filename="drifting.gif")
```

![drifting animation](drifting.gif)

## Slip-based planar wheel
This example demonstrates use of the [`PlanarMechanics.SlipBasedWheelJoint`](@ref) component, which is a more advanced 2D wheel component with slip-dependent friction characteristics. The wheel is being driven by a constant torque, and is connected through a [`PlanarMechanics.Prismatic`](@ref) joint to a [`PlanarMechanics.Revolute`](@ref) joint. This forces the wheel to move in a circular arc around the revolute pivot point, and spiral outwards due to slip. A [`PlanarMechanics.Body](@ref) attached to the end of the prismatic joint is used to add inertial properties.

```@example WHEEL
@mtkmodel TestSlipBasedWheel begin
    @components begin
        slipBasedWheelJoint = Pl.SlipBasedWheelJoint(
            radius = 0.3,
            r = [1,0],              # Driving direction at angle phi = 0
            mu_A = 0.8,             # Friction coefficient at adhesion
            mu_S = 0.4,             # Friction coefficient at sliding
            N = 100,                # Base normal load
            sAdhesion = 0.04,       # Adhesion slippage
            sSlide = 0.12,          # Sliding slippage
            vAdhesion_min = 0.05,   # Minimum adhesion velocity
            vSlide_min = 0.15,      # Minimum sliding velocity
            color = tire_black,
        )
        prismatic = Pl.Prismatic(r = [0,1], s = 1, v = 0)
        revolute = Pl.Revolute(phi = 0, w = 0)
        fixed = Pl.Fixed()
        engineTorque = Rotational.ConstantTorque(tau_constant = 2)
        body = Pl.Body(m = 10, I = 1, gy=0)
        inertia = Rotational.Inertia(J = 1, phi = 0, w = 0)
        constant = Blocks.Constant(k = 0)
    end
    @equations begin
        connect(fixed.frame, revolute.frame_a)
        connect(revolute.frame_b, prismatic.frame_a)
        connect(prismatic.frame_b, body.frame_a)
        connect(prismatic.frame_b, slipBasedWheelJoint.frame_a)
        connect(slipBasedWheelJoint.flange_a, inertia.flange_b)
        connect(constant.output, slipBasedWheelJoint.dynamicLoad)
        connect(engineTorque.flange, inertia.flange_a)
    end
end

@named model = TestSlipBasedWheel()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
display(unknowns(ssys))
defs = ModelingToolkit.defaults(model)
prob = ODEProblem(ssys, [
    model.inertia.w => 1e-10, # This is important, at zero velocity, the friction is ill-defined
    model.revolute.frame_b.phi => 0,
    D(model.revolute.frame_b.phi) => 0,
    D(model.prismatic.r0[2]) => 0,
], (0.0, 15.0))
sol = solve(prob, Rodas5Pr())
render(model, sol, show_axis=false, x=0, y=0, z=4, traces=[model.slipBasedWheelJoint.frame_a], filename="slipwheel.gif")
```

![slipwheel animation](slipwheel.gif)
