using ModelingToolkit
using Multibody
using Test
using SymbolicIR
t = Multibody.t

cd(@__DIR__)
world = Multibody.world
include("OneAxis.jl")
include("FullRobot.jl")
@named structure = MechanicalStructure()
@named motor = Motor()
@named controller = Controller()
@named axis2 = AxisType2()
@named gear2 = GearType2()
# @named axis1 = AxisType1()
@named gear1 = GearType1()

@named pp = PathPlanning1(;)
@named pp6 = PathPlanning6(;)

@named oneaxis = OneAxis()

ssys = structural_simplify(IRSystem(oneaxis))
ssys = structural_simplify(oneaxis, allow_parameters = false)

@named robot = FullRobot()

ssys = structural_simplify(robot, allow_parameters = false)
ssys = structural_simplify(IRSystem(robot))



dummyder = setdiff(states(ssys), states(oneaxis))
op = merge(ModelingToolkit.defaults(oneaxis), Dict(x => 0.0 for x in dummyder))
prob = ODEProblem(ssys, op, (0.0, 10.0))

using OrdinaryDiffEq
sol = solve(prob, Rodas4(), u0=prob.u0 .+ 0.01.*randn.())


