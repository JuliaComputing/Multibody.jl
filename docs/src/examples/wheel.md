# Wheels

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

