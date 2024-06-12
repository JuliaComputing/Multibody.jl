"""
    RobotAxis(; name, mLoad = 15, kp = 5, ks = 0.5, Ts = 0.05, q0 = 0, q1 = 120, swingtime = 0.5, refSpeedMax = 3, refAccMax = 10, kwargs)

A single robot axis.

- `mLoad`: Mass of load
- `kp`: Proportional gain of position controller
- `ks`: Proportional gain of velocity controller
- `Ts`: Time constant of integrator of velocity controller
- `q0`: Start angle in degrees
- `q1`: End angle in degrees
"""
function RobotAxis(; name, mLoad = 15, kp = 5.0, ks = 0.5, Ts = 0.05, q0 = 0,
    q1 = 120, swingtime = 0.5, refSpeedMax = 3, refAccMax = 10, kwargs...)

    # parameter SI.Mass mLoad(min=0)=15 "Mass of load";
    # parameter Real kp=5 "Gain of position controller of axis";
    # parameter Real ks=0.5 "Gain of speed controller of axis";
    # parameter SI.Time Ts=0.05
    #   "Time constant of integrator of speed controller of axis";
    # parameter Modelica.Units.NonSI.Angle_deg q0 = 0 "Start angle of axis";
    # parameter Modelica.Units.NonSI.Angle_deg q1 = 120 "End angle of axis";

    # parameter SI.Time swingtime=0.5
    #   "Additional time after reference motion is in rest before simulation is stopped";
    # parameter SI.AngularVelocity refSpeedMax=3 "Maximum reference speed";
    # parameter SI.AngularAcceleration refAccMax=10
    #   "Maximum reference acceleration";

    @parameters begin
    mLoad = mLoad, [description = "Mass of load"]
    kp = kp, [description = "Gain of position controller of axis"]
    ks = ks, [description = "Gain of speed controller of axis"]
    Ts = Ts, [description = "Time constant of integrator of speed controller of axis"]
    # q0 = q0, [description = "Start angle of axis"]
    # q1 = q1, [description = "End angle of axis"]
    swingtime = swingtime,
        [
            description = "Additional time after reference motion is in rest before simulation is stopped",
        ]
    # refSpeedMax = refSpeedMax, [description = "Maximum reference speed"]
    # refAccMax = refAccMax, [description = "Maximum reference acceleration"]
    end

    systems = @named begin
    axis = AxisType1(w = 5500,
                ratio = 210,
                c = 8,
                cd = 0.01,
                Rv0 = 0.5,
                Rv1 = (0.1 / 130),
                kp = kp,
                ks = ks,
                Ts = Ts)
    load = Rotational.Inertia(J = 1.3 * mLoad)
    pathPlanning = PathPlanning1(;swingtime = swingtime,
                            q0deg = q0,
                            q1deg = q1,
                            speed_max = refSpeedMax,
                            acc_max = refAccMax,
                            kwargs...
                            )
    controlBus = ControlBus()
    end
    eqs = [
        connect(axis.flange, load.flange_a),
        connect(pathPlanning.controlBus, controlBus),
        connect(controlBus.axisControlBus1, axis.axisControlBus),
    ]
    ODESystem(eqs, t; systems, name)
end


"""
    Robot6DOF(; name, kwargs)

A model of a 6DOF serial industrial robot with revolute joints and cascade P/PI position/velocity controllers.
"""
function Robot6DOF(; name, kwargs...)
    @parameters begin
        # mLoad = 15, [description = "Mass of load"]
        # rLoad[1:3] = [0.1, 0.25, 0.1],
        #              [description = "Distance from last flange to load mass"]
        # g = 9.81, [description = "Gravity acceleration"]
        # refStartTime = 0, [description = "Start time of reference motion"]
        # refSwingTime = 0.5,
        #                [
        #                    description = "Additional time after reference motion is in rest before simulation is stopped",
        #                ]
        # q01 = -60, [description = "Start angle of axis 1"]
        # q02 = 20, [description = "Start angle of axis 2"]
        # q03 = 90, [description = "Start angle of axis 3"]
        # q04 = 0, [description = "Start angle of axis 4"]
        # q05 = -110, [description = "Start angle of axis 5"]
        # q06 = 0, [description = "Start angle of axis 6"]
        # q11 = 60, [description = "End angle of axis 1"]
        # q12 = -70, [description = "End angle of axis 2"]
        # q13 = -35, [description = "End angle of axis 3"]
        # q14 = 45, [description = "End angle of axis 4"]
        # q15 = 110, [description = "End angle of axis 5"]
        # q16 = 45, [description = "End angle of axis 6"]
        # refSpeedMax[1:6] = [3, 1.5, 5, 3.1, 3.1, 4.1],
        #                    [description = "Maximum reference speeds of all joints"]
        # refAccMax[1:6] = [15, 15, 15, 60, 60, 60],
        #                  [description = "Maximum reference accelerations of all joints"]
        # kp1 = 5, [description = "Gain of position controller"]
        # ks1 = 0.5, [description = "Gain of speed controller"]
        # Ts1 = 0.05, [description = "Time constant of integrator of speed controller"]
        # kp2 = 5, [description = "Gain of position controller"]
        # ks2 = 0.5, [description = "Gain of speed controller"]
        # Ts2 = 0.05, [description = "Time constant of integrator of speed controller"]
        # kp3 = 5, [description = "Gain of position controller"]
        # ks3 = 0.5, [description = "Gain of speed controller"]
        # Ts3 = 0.05, [description = "Time constant of integrator of speed controller"]
        # kp4 = 5, [description = "Gain of position controller"]
        # ks4 = 0.5, [description = "Gain of speed controller"]
        # Ts4 = 0.05, [description = "Time constant of integrator of speed controller"]
        # kp5 = 5, [description = "Gain of position controller"]
        # ks5 = 0.5, [description = "Gain of speed controller"]
        # Ts5 = 0.05, [description = "Time constant of integrator of speed controller"]
        # kp6 = 5, [description = "Gain of position controller"]
        # ks6 = 0.5, [description = "Gain of speed controller"]
        # Ts6 = 0.05, [description = "Time constant of integrator of speed controller"]
    end

    mLoad = 15 #, [description = "Mass of load"]
    rLoad = [0.1, 0.25, 0.1] #,                   [description = "Distance from last flange to load mass"]
    g = 9.81 #, [description = "Gravity acceleration"]
    refStartTime = 0 #, [description = "Start time of reference motion"]
    refSwingTime = 0.5 #, [description = "Additional time after reference motion is in rest before simulation is stopped",                   ]
    refSpeedMax = [3, 1.5, 5, 3.1, 3.1, 4.1]    #  [description = "Maximum reference speeds of all joints"]
    refAccMax = [15, 15, 15, 60, 60, 60]  #  [description = "Maximum reference accelerations of all joints"]
    kp1 = 5 #, [description = "Gain of position controller"]
    ks1 = 0.5 #, [description = "Gain of speed controller"]
    Ts1 = 0.05 #, [description = "Time constant of integrator of speed controller"]
    kp2 = 5 #, [description = "Gain of position controller"]
    ks2 = 0.5 #, [description = "Gain of speed controller"]
    Ts2 = 0.05 #, [description = "Time constant of integrator of speed controller"]
    kp3 = 5 #, [description = "Gain of position controller"]
    ks3 = 0.5 #, [description = "Gain of speed controller"]
    Ts3 = 0.05 #, [description = "Time constant of integrator of speed controller"]
    kp4 = 5 #, [description = "Gain of position controller"]
    ks4 = 0.5 #, [description = "Gain of speed controller"]
    Ts4 = 0.05 #, [description = "Time constant of integrator of speed controller"]
    kp5 = 5 #, [description = "Gain of position controller"]
    ks5 = 0.5 #, [description = "Gain of speed controller"]
    Ts5 = 0.05 #, [description = "Time constant of integrator of speed controller"]
    kp6 = 5 #, [description = "Gain of position controller"]
    ks6 = 0.5 #, [description = "Gain of speed controller"]
    Ts6 = 0.05 #, [description = "Time constant of integrator of speed controller"]

    q01 = -60 # Can't yet have these as parameters
    q02 = 20 # Can't yet have these as parameters
    q03 = 90 # Can't yet have these as parameters
    q04 = 0 # Can't yet have these as parameters
    q05 = -110 # Can't yet have these as parameters
    q06 = 0 # Can't yet have these as parameters
    q11 = 60 # Can't yet have these as parameters
    q12 = -70 # Can't yet have these as parameters
    q13 = -35 # Can't yet have these as parameters
    q14 = 45 # Can't yet have these as parameters
    q15 = 110 # Can't yet have these as parameters
    q16 = 45 # Can't yet have these as parameters

    systems = @named begin
        mechanics = MechanicalStructure(mLoad = (mLoad),
                                        rLoad = (rLoad),
                                        g = (g))
        pathPlanning = PathPlanning6(;naxis = 6,
                                     q0deg = [
                                         q01,
                                         q02,
                                         q03,
                                         q04,
                                         q05,
                                         q06,
                                     ],
                                     q1deg = [
                                         q11,
                                         q12,
                                         q13,
                                         q14,
                                         q15,
                                         q16,
                                     ],
                                     speed_max = refSpeedMax,
                                     acc_max = refAccMax,
                                     startTime = refStartTime,
                                     swingtime = refSwingTime,
                                     kwargs...)

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
