using ModelingToolkit
using Multibody
using JuliaSimCompiler
using OrdinaryDiffEq
using OrdinaryDiffEq.SciMLBase: successful_retcode
using BenchmarkTools

using Multibody: Robot6DOF

@info "Building Robot6DOF"
@named robot = Robot6DOF()
robot = complete(robot)

const ic = [
    robot.mechanics.r1.phi => deg2rad(-60),
    robot.mechanics.r2.phi => deg2rad(20),
    robot.mechanics.r3.phi => deg2rad(90),
    robot.mechanics.r4.phi => deg2rad(0),
    robot.mechanics.r5.phi => deg2rad(-110),
    robot.mechanics.r6.phi => deg2rad(0),
    robot.axis1.motor.Jmotor.phi => deg2rad(-60) * (-105),
    robot.axis2.motor.Jmotor.phi => deg2rad(20) * (210),
    robot.axis3.motor.Jmotor.phi => deg2rad(90) * (60),
]

println("\n========== WARM-UP RUN (includes compilation) ==========")
@time "structural_simplify (warm-up)" ssys = structural_simplify(IRSystem(robot))
@time "ODEProblem creation (warm-up)" prob = ODEProblem(ssys, ic, (0.0, 2.0))
@time "solve (warm-up)"               sol  = solve(prob, Rodas5P(autodiff=false))
@assert successful_retcode(sol)

println("\n========== TIMED RUN (cold structural_simplify/ODEProblem) ==========")
@named robot2 = Robot6DOF()
robot2 = complete(robot2)
@time "structural_simplify" ssys2 = structural_simplify(IRSystem(robot2))
ic2 = [
    robot2.mechanics.r1.phi => deg2rad(-60),
    robot2.mechanics.r2.phi => deg2rad(20),
    robot2.mechanics.r3.phi => deg2rad(90),
    robot2.mechanics.r4.phi => deg2rad(0),
    robot2.mechanics.r5.phi => deg2rad(-110),
    robot2.mechanics.r6.phi => deg2rad(0),
    robot2.axis1.motor.Jmotor.phi => deg2rad(-60) * (-105),
    robot2.axis2.motor.Jmotor.phi => deg2rad(20) * (210),
    robot2.axis3.motor.Jmotor.phi => deg2rad(90) * (60),
]
@time "ODEProblem creation" prob2 = ODEProblem(ssys2, ic2, (0.0, 2.0))
@time "first solve" sol2 = solve(prob2, Rodas5P(autodiff=false))
@assert successful_retcode(sol2)

println("\n========== BenchmarkTools @btime (solve only, already-compiled) ==========")
@btime solve($prob2, Rodas5P(autodiff=false)) samples=5 evals=1

println("\nDone.")
