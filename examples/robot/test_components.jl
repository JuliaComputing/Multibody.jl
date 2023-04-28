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

@named pp = PathPlanning1(; )
@named pp6 = PathPlanning6(; )


@named oneaxis = OneAxis()

ssys = structural_simplify(IRSystem(oneaxis))
ssys = structural_simplify(oneaxis, allow_parameters = false)

@named robot = FullRobot()

# ssys = structural_simplify(IRSystem(robot))
ssys = structural_simplify(robot, allow_parameters = false)

prob = ODEProblem(ssys, [], (0.0, 10.0))

sol = solve(prob, Rodas4())