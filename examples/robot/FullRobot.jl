include("utilities/components.jl")
include("utilities/path_planning.jl")

function FullRobot(; name)
    @parameters begin
        mLoad = 15, [description = "Mass of load"]
        rLoad[1:3] = [0.1, 0.25, 0.1],
                     [description = "Distance from last flange to load mass"]
        g = 9.81, [description = "Gravity acceleration"]
        refStartTime = 0, [description = "Start time of reference motion"]
        refSwingTime = 0.5,
                       [
                           description = "Additional time after reference motion is in rest before simulation is stopped",
                       ]
        # startAngle1 = -60, [description = "Start angle of axis 1"]
        # startAngle2 = 20, [description = "Start angle of axis 2"]
        # startAngle3 = 90, [description = "Start angle of axis 3"]
        # startAngle4 = 0, [description = "Start angle of axis 4"]
        # startAngle5 = -110, [description = "Start angle of axis 5"]
        # startAngle6 = 0, [description = "Start angle of axis 6"]
        # endAngle1 = 60, [description = "End angle of axis 1"]
        # endAngle2 = -70, [description = "End angle of axis 2"]
        # endAngle3 = -35, [description = "End angle of axis 3"]
        # endAngle4 = 45, [description = "End angle of axis 4"]
        # endAngle5 = 110, [description = "End angle of axis 5"]
        # endAngle6 = 45, [description = "End angle of axis 6"]
        refSpeedMax[1:6] = [3, 1.5, 5, 3.1, 3.1, 4.1],
                           [description = "Maximum reference speeds of all joints"]
        refAccMax[1:6] = [15, 15, 15, 60, 60, 60],
                         [description = "Maximum reference accelerations of all joints"]
        kp1 = 5, [description = "Gain of position controller"]
        ks1 = 0.5, [description = "Gain of speed controller"]
        Ts1 = 0.05, [description = "Time constant of integrator of speed controller"]
        kp2 = 5, [description = "Gain of position controller"]
        ks2 = 0.5, [description = "Gain of speed controller"]
        Ts2 = 0.05, [description = "Time constant of integrator of speed controller"]
        kp3 = 5, [description = "Gain of position controller"]
        ks3 = 0.5, [description = "Gain of speed controller"]
        Ts3 = 0.05, [description = "Time constant of integrator of speed controller"]
        kp4 = 5, [description = "Gain of position controller"]
        ks4 = 0.5, [description = "Gain of speed controller"]
        Ts4 = 0.05, [description = "Time constant of integrator of speed controller"]
        kp5 = 5, [description = "Gain of position controller"]
        ks5 = 0.5, [description = "Gain of speed controller"]
        Ts5 = 0.05, [description = "Time constant of integrator of speed controller"]
        kp6 = 5, [description = "Gain of position controller"]
        ks6 = 0.5, [description = "Gain of speed controller"]
        Ts6 = 0.05, [description = "Time constant of integrator of speed controller"]
    end

    startAngle1 = -60 # Can't yet have these as parameters
    startAngle2 = 20 # Can't yet have these as parameters
    startAngle3 = 90 # Can't yet have these as parameters
    startAngle4 = 0 # Can't yet have these as parameters
    startAngle5 = -110 # Can't yet have these as parameters
    startAngle6 = 0 # Can't yet have these as parameters
    endAngle1 = 60 # Can't yet have these as parameters
    endAngle2 = -70 # Can't yet have these as parameters
    endAngle3 = -35 # Can't yet have these as parameters
    endAngle4 = 45 # Can't yet have these as parameters
    endAngle5 = 110 # Can't yet have these as parameters
    endAngle6 = 45 # Can't yet have these as parameters

    systems = @named begin
        mechanics = MechanicalStructure(mLoad = mLoad,
                                        rLoad = rLoad,
                                        g = g)
        pathPlanning = PathPlanning6(naxis = 6,
                                     angleBegDeg = [
                                         startAngle1,
                                         startAngle2,
                                         startAngle3,
                                         startAngle4,
                                         startAngle5,
                                         startAngle6,
                                     ],
                                     angleEndDeg = [
                                         endAngle1,
                                         endAngle2,
                                         endAngle3,
                                         endAngle4,
                                         endAngle5,
                                         endAngle6,
                                     ],
                                     speedMax = refSpeedMax,
                                     accMax = refAccMax,
                                     startTime = refStartTime,
                                     swingTime = refSwingTime)

        axis1 = AxisType1(w = 4590,
                          ratio = -105,
                          c = 43,
                          cd = 0.005,
                          Rv0 = 0.4,
                          Rv1 = (0.13 / 160),
                          kp = kp1,
                          ks = ks1,
                          Ts = Ts1)
        axis2 = AxisType1(w = 5500,
                          ratio = 210,
                          c = 8,
                          cd = 0.01,
                          Rv1 = (0.1 / 130),
                          Rv0 = 0.5,
                          kp = kp2,
                          ks = ks2,
                          Ts = Ts2)

        axis3 = AxisType1(w = 5500,
                          ratio = 60,
                          c = 58,
                          cd = 0.04,
                          Rv0 = 0.7,
                          Rv1 = (0.2 / 130),
                          kp = kp3,
                          ks = ks3,
                          Ts = Ts3)
        axis4 = AxisType2(k = 0.2365,
                          w = 6250,
                          D = 0.55,
                          J = 1.6e-4,
                          ratio = -99,
                          Rv0 = 21.8,
                          Rv1 = 9.8,
                          peak = 26.7 / 21.8,
                          kp = kp4,
                          ks = ks4,
                          Ts = Ts4)
        axis5 = AxisType2(k = 0.2608,
                          w = 6250,
                          D = 0.55,
                          J = 1.8e-4,
                          ratio = 79.2,
                          Rv0 = 30.1,
                          Rv1 = 0.03,
                          peak = 39.6 / 30.1,
                          kp = kp5,
                          ks = ks5,
                          Ts = Ts5)
        axis6 = AxisType2(k = 0.0842,
                          w = 7400,
                          D = 0.27,
                          J = 4.3e-5,
                          ratio = -99,
                          Rv0 = 10.9,
                          Rv1 = 3.92,
                          peak = 16.8 / 10.9,
                          kp = kp6,
                          ks = ks6,
                          Ts = Ts6)

        controlBus = ControlBus()
    end

    eqs = [connect(axis2.flange, mechanics.axis2)
           connect(axis1.flange, mechanics.axis1)
           connect(axis3.flange, mechanics.axis3)
           connect(axis4.flange, mechanics.axis4)
           connect(axis5.flange, mechanics.axis5)
           connect(axis6.flange, mechanics.axis6)
           connect(controlBus, pathPlanning.controlBus)
           connect(controlBus.axisControlBus1, axis1.axisControlBus)
           connect(controlBus.axisControlBus2, axis2.axisControlBus)
           connect(controlBus.axisControlBus3, axis3.axisControlBus)
           connect(controlBus.axisControlBus4, axis4.axisControlBus)
           connect(controlBus.axisControlBus5, axis5.axisControlBus)
           connect(controlBus.axisControlBus6, axis6.axisControlBus)]

    ODESystem(eqs, t; systems, name)
end
