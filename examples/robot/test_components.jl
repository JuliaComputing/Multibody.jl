using Multibody
cd(@__DIR__)
t = Multibody.t
world = Multibody.world
# include("OneAxis.jl")
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


@named robot = FullRobot()