include("utilities/components.jl")

function OneAxis(; name, mLoad = 15, kp = 5, ks = 0.5, Ts = 0.05, startAngle = 0,
                 endAngle = 120, swingTime = 0.5, refSpeedMax = 3, refAccMax = 10)

    # parameter SI.Mass mLoad(min=0)=15 "Mass of load";
    # parameter Real kp=5 "Gain of position controller of axis";
    # parameter Real ks=0.5 "Gain of speed controller of axis";
    # parameter SI.Time Ts=0.05
    #   "Time constant of integrator of speed controller of axis";
    # parameter Modelica.Units.NonSI.Angle_deg startAngle = 0 "Start angle of axis";
    # parameter Modelica.Units.NonSI.Angle_deg endAngle = 120 "End angle of axis";

    # parameter SI.Time swingTime=0.5
    #   "Additional time after reference motion is in rest before simulation is stopped";
    # parameter SI.AngularVelocity refSpeedMax=3 "Maximum reference speed";
    # parameter SI.AngularAcceleration refAccMax=10
    #   "Maximum reference acceleration";

    @parameters begin
        mLoad = mLoad, [description = "Mass of load"]
        kp = kp, [description = "Gain of position controller of axis"]
        ks = ks, [description = "Gain of speed controller of axis"]
        Ts = Ts, [description = "Time constant of integrator of speed controller of axis"]
        # startAngle = startAngle, [description = "Start angle of axis"]
        # endAngle = endAngle, [description = "End angle of axis"]
        swingTime = swingTime,
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
        pathPlanning = PathPlanning1(swingTime = swingTime,
                                     angleBegDeg = startAngle,
                                     angleEndDeg = endAngle
                                     #  speedMax = refSpeedMax,
                                     #  accMax = refAccMax
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
