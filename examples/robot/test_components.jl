@named structure = MechanicalStructure()
@named motor = Motor()
@named controller = Controller()
@named axis2 = AxisType2()
@named gear2 = GearType2()
# @named axis1 = AxisType1()
@named gear1 = GearType1()

PathPlanning1(; name, angleBegDeg = 0, angleEndDeg = 1, speedMax = 3,
              accMax = 2.5, startTime = 0, swingTime = 0.5)
