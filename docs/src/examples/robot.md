# Industrial robot

![animation](robot.gif)

```@example robot
using Multibody
using ModelingToolkit
using Plots
using JuliaSimCompiler
using OrdinaryDiffEq

t = Multibody.t
D = Differential(t)
world = Multibody.world
@named robot = Robot6DOF(trivial=true)
robot = complete(robot)


ssys = structural_simplify(IRSystem(robot))
prob = ODEProblem(ssys, [], (0.0, 4.0))
sol = solve(prob, Rodas5P(autodiff=false));
@test SciMLBase.successful_retcode(sol)


plot(sol, idxs = [
    robot.pathPlanning.controlBus.axisControlBus1.angle_ref
    robot.pathPlanning.controlBus.axisControlBus2.angle_ref
    robot.pathPlanning.controlBus.axisControlBus3.angle_ref
    robot.pathPlanning.controlBus.axisControlBus4.angle_ref
    robot.pathPlanning.controlBus.axisControlBus5.angle_ref
    robot.pathPlanning.controlBus.axisControlBus6.angle_ref
], layout=6, size=(800,800), l=(:black, :dash), legend=false)
plot!(sol, idxs = [
    robot.pathPlanning.controlBus.axisControlBus1.angle
    robot.pathPlanning.controlBus.axisControlBus2.angle
    robot.pathPlanning.controlBus.axisControlBus3.angle
    robot.pathPlanning.controlBus.axisControlBus4.angle
    robot.pathPlanning.controlBus.axisControlBus5.angle
    robot.pathPlanning.controlBus.axisControlBus6.angle
], layout=6)
display(current())
```



## 3D animation
Multibody.jl supports automatic 3D rendering of mechanisms, we use this feature to illustrate the result of the simulation below:

```@example robot
import CairoMakie
Multibody.render(robot, sol; z = -5, filename = "robot.gif")
nothing # hide
```

![animation](robot.gif)
