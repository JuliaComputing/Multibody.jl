using ModelingToolkit
using Multibody
using Test
using JuliaSimCompiler
using OrdinaryDiffEq
t = Multibody.t
D = Differential(t)

@named robot = Multibody.Robot6DOF()
robot = complete(robot)
println("complete")

ssys = structural_simplify(IRSystem(robot))
println("structural_simplify")


prob = ODEProblem(ssys, [
    robot.mechanics.r1.phi => deg2rad(-60)
    robot.mechanics.r2.phi => deg2rad(20)
    robot.mechanics.r3.phi => deg2rad(90)
    robot.mechanics.r4.phi => deg2rad(0)
    robot.mechanics.r5.phi => deg2rad(-110)
    robot.mechanics.r6.phi => deg2rad(0)

    robot.axis1.motor.Jmotor.phi => deg2rad(-60) * (-105) # Multiply by gear ratio
    robot.axis2.motor.Jmotor.phi => deg2rad(20) * (210)
    robot.axis3.motor.Jmotor.phi => deg2rad(90) * (60)
], (0.0, 2.0))
println("ODEProblem done, calling solve")

sol = solve(prob, Rodas5P(autodiff=false));