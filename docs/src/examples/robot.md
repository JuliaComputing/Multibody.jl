# Industrial robot

![animation](robot.gif)

```@example robot
using Multibody
using ModelingToolkit
using Plots
# using JuliaSimCompiler
using OrdinaryDiffEq
using Test

t = Multibody.t
D = Differential(t)
@named robot = Multibody.Robot6DOF()
robot = complete(robot)

length(equations(robot))
```
The robot is a medium sized system with some 2000 equations before simplification.

After simplification, the following states are chosen:
```@example robot
ssys = structural_simplify(multibody(robot))
unknowns(ssys)
```
    
```@example robot
prob = ODEProblem(ssys, Dict([
    robot.mechanics.r1.phi => deg2rad(-60)
    robot.mechanics.r2.phi => deg2rad(20)
    robot.mechanics.r3.phi => deg2rad(90)
    robot.mechanics.r4.phi => deg2rad(0)
    robot.mechanics.r5.phi => deg2rad(-110)
    robot.mechanics.r6.phi => deg2rad(0)

    robot.axis1.motor.Jmotor.phi => deg2rad(-60) *  -105 # Multiply by gear ratio
    robot.axis2.motor.Jmotor.phi => deg2rad(20) *  210
    robot.axis3.motor.Jmotor.phi => deg2rad(90) *  60
]), (0.0, 2.0))
sol = solve(prob, Rodas5P(autodiff=false));
@test SciMLBase.successful_retcode(sol)

plot(sol, idxs = [
    robot.pathPlanning.controlBus.axisControlBus1.angle_ref
    robot.pathPlanning.controlBus.axisControlBus2.angle_ref
    robot.pathPlanning.controlBus.axisControlBus3.angle_ref
    robot.pathPlanning.controlBus.axisControlBus4.angle_ref
    robot.pathPlanning.controlBus.axisControlBus5.angle_ref
    robot.pathPlanning.controlBus.axisControlBus6.angle_ref
], layout=(4,3), size=(800,800), l=(:black, :dash), legend=:outertop, legendfontsize=6)
plot!(sol, idxs = [
    robot.pathPlanning.controlBus.axisControlBus1.angle
    robot.pathPlanning.controlBus.axisControlBus2.angle
    robot.pathPlanning.controlBus.axisControlBus3.angle
    robot.pathPlanning.controlBus.axisControlBus4.angle
    robot.pathPlanning.controlBus.axisControlBus5.angle
    robot.pathPlanning.controlBus.axisControlBus6.angle
], sp=1:6)

plot!(sol, idxs = [
    robot.axis1.controller.feedback1.output.u
    robot.axis2.controller.feedback1.output.u
    robot.axis3.controller.feedback1.output.u
    robot.axis4.controller.feedback1.output.u
    robot.axis5.controller.feedback1.output.u
    robot.axis6.controller.feedback1.output.u
], sp=7:12, lab="Position error", link=:x)
plot!(xlabel=[fill("", 1, 9) fill("Time [s]", 1, 3)])
```
We see that after an initial transient, the robot controller converges to tracking the reference trajectory well. However, since the first three axes of the robot are modeled as slightly flexible, and we are ultimately interested in the tracking performance _after_ the gear box and any flexibilities it may suffer from, we plot also this tracking error
```@example robot
plot(sol, idxs = [
    robot.axis1.controller.feedback1.output.u / ( -105)
    robot.axis2.controller.feedback1.output.u / (210)
    robot.axis3.controller.feedback1.output.u / (60)
], layout=(1,3), lab="Position error, motor side", link=:x)
plot!(sol, idxs = [
            robot.pathPlanning.controlBus.axisControlBus1.angle_ref - robot.mechanics.r1.phi #
            robot.pathPlanning.controlBus.axisControlBus2.angle_ref - robot.mechanics.r2.phi #
            robot.pathPlanning.controlBus.axisControlBus3.angle_ref - robot.mechanics.r3.phi #
], lab="Position error, arm side")
```


## Trajectory planning
In the example, the robot is tracking a reference trajectory generated using the function [`point_to_point`](@ref) and interfaced from the component [`KinematicPTP`](@ref). We can inspect the generated trajectory by plotting the positions, velocities and accelerations (we show one joint only to keep the plot size limited):
```@example robot
plot(sol, idxs = [
            robot.pathPlanning.controlBus.axisControlBus1.angle_ref
            robot.pathPlanning.controlBus.axisControlBus1.speed_ref
            robot.pathPlanning.controlBus.axisControlBus1.acceleration_ref
], layout=(3,1), lab=["\$q\$" "\$\\dot q\$" "\$\\ddot q\$"], xlabel=["" "" "Time [s]"])
```

## 3D animation
Multibody.jl supports automatic 3D rendering of mechanisms, we use this feature to illustrate the result of the simulation below:

```@example robot
import GLMakie
Multibody.render(robot, sol; y=2, lookat=[0,1,0], filename = "robot.gif")
nothing # hide
```

![animation](robot.gif)

## Kinematics
The coordinates of any point on the mechanism may be obtained in the world coordinate frame by either

```@example robot
output = collect(robot.mechanics.b6.frame_b.r_0)
fkine = build_explicit_observed_function(ssys, output)
fkine(prob.u0, prob.p, 0)
```

Query the solution object with the desired output, e.g.,
```@example robot
sol(0, idxs=output)
```
query the problem with the output, in which case the initial condition is used to compute the output
```@example robot
prob[output]
```

or by building an explicit function `(state, parameters, time) -> output`
```@example robot
fkine = build_explicit_observed_function(ssys, output)
fkine(prob.u0, prob.p, 0)
```
!!! note
    The function `fkine` above takes the full state of the robot model, as opposed to only the joint angles.

```@setup
# Temporarily removed due to https://github.com/JuliaComputing/JuliaSimCompiler.jl/issues/366
### Jacobian
# We can compute the Jacobian ``J`` of the forward-kinematics function using the package ForwardDiff (this Jacobian is often referred to as the _analytical Jacobian_, which in the 6DOF case is different from the _geometrical Jacobian_ that is used in the relation ``v = J\dot{q}``). The Jacobian of the end-effector positional coordinates will be a 3×36 matrix, since we have 36-dimensional state of the robot after simplification. Since the end-effector coordinates do not depend on all the state variables, we may ask which variables it depends on by finding non-zero columns of ``J``
# ```@example robot
# using ModelingToolkit.ForwardDiff
# J = ForwardDiff.jacobian(x->fkine(x, prob.p, 0), prob.u0)
# nonzero_inds = findall(any(!iszero, J, dims=1)[:])
# unknowns(ssys)[nonzero_inds]
# ```
# We see that the end-effector position depends on all mechanical angles except for the last one, which is expected since the end-effector origin is on the axis of rotation of joint 6. 
```