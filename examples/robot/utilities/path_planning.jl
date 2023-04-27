"Generate reference angles for fastest kinematic movement"
function PathPlanning1(; name, angleBegDeg = 0, angleEndDeg = 1, speedMax = 3,
                       accMax = 2.5, startTime = 0, swingTime = 0.5)
    @parameters begin
        angleBegDeg = angleBegDeg, [description = "Start angle"]
        angleEndDeg = angleEndDeg, [description = "End angle"]
        speedMax = speedMax, [description = "Maximum axis speed"]
        accMax = accMax, [description = "Maximum axis acceleration"]
        startTime = startTime, [description = "Start time of movement"]
        swingTime = swingTime,
                    [
                        description = "Additional time after reference motion is in rest before simulation is stopped",
                    ]
        angleBeg = deg2rad(angleBegDeg), [description = "Start angles"]
        angleEnd = deg2rad(angleEndDeg), [description = "End angles"]
    end

    systems = @named begin
        controlBus = ControlBus()
        path = KinematicPTP2(q_end = angleEnd,
                             qd_max = speedMax,
                             qdd_max = accMax,
                             startTime = startTime,
                             q_begin = angleBeg)
        pathToAxis1 = PathToAxisControlBus(nAxis = 1, axisUsed = 1)
        terminateSimulation = TerminateSimulation(condition = time >=
                                                              path.endTime + swingTime)
    end

    eqs = [connect(path.q, pathToAxis1.q)
           connect(path.qd, pathToAxis1.qd)
           connect(path.qdd, pathToAxis1.qdd)
           connect(path.moving, pathToAxis1.moving)
           connect(pathToAxis1.axisControlBus, controlBus.axisControlBus1)]
    ODESystem(eqs, t; name, systems)
end

function PathPlanning6(; name, naxis = 6, angleBegDeg = zeros(naxis),
                       angleEndDeg = ones(naxis), speedMax = fill(3, naxis),
                       accMax = fill(2.5, naxis), startTime = 0, swingTime = 0.5)
    @parameters begin
        naxis = naxis, [description = "Number of driven axis"]
        angleBegDeg[1:naxis] = angleBegDeg, [description = "Start angles"]
        angleEndDeg[1:naxis] = angleEndDeg, [description = "End angles"]
        speedMax[1:naxis] = speedMax, [description = "Maximum axis speed"]
        accMax[1:naxis] = accMax, [description = "Maximum axis acceleration"]
        startTime = startTime, [description = "Start time of movement"]
        swingTime = swingTime,
                    [
                        description = "Additional time after reference motion is in rest before simulation is stopped",
                    ]
        angleBeg[1:6] = deg2rad.(angleBegDeg), [description = "Start angles"]
        angleEnd[1:6] = deg2rad.(angleEndDeg), [description = "End angles"]
    end

    systems = @named begin
        controlBus = ControlBus()
        path = KinematicPTP2(q_end = angleEnd,
                             qd_max = speedMax,
                             qdd_max = accMax,
                             startTime = startTime,
                             q_begin = angleBeg)
        pathToAxis1 = PathToAxisControlBus(nAxis = naxis, axisUsed = 1)
        pathToAxis2 = PathToAxisControlBus(nAxis = naxis, axisUsed = 2)
        pathToAxis3 = PathToAxisControlBus(nAxis = naxis, axisUsed = 3)
        pathToAxis4 = PathToAxisControlBus(nAxis = naxis, axisUsed = 4)
        pathToAxis5 = PathToAxisControlBus(nAxis = naxis, axisUsed = 5)
        pathToAxis6 = PathToAxisControlBus(nAxis = naxis, axisUsed = 6)
        terminateSimulation = TerminateSimulation(condition = time >=
                                                              path.endTime + swingTime)
    end

    eqs = [connect(path.q, pathToAxis1.q)
           connect(path.qd, pathToAxis1.qd)
           connect(path.qdd, pathToAxis1.qdd)
           connect(path.moving, pathToAxis1.moving)
           connect(path.q, pathToAxis2.q)
           connect(path.qd, pathToAxis2.qd)
           connect(path.qdd, pathToAxis2.qdd)
           connect(path.moving, pathToAxis2.moving)
           connect(path.q, pathToAxis3.q)
           connect(path.qd, pathToAxis3.qd)
           connect(path.qdd, pathToAxis3.qdd)
           connect(path.moving, pathToAxis3.moving)
           connect(path.q, pathToAxis4.q)
           connect(path.qd, pathToAxis4.qd)
           connect(path.qdd, pathToAxis4.qdd)
           connect(path.moving, pathToAxis4.moving)
           connect(path.q, pathToAxis5.q)
           connect(path.qd, pathToAxis5.qd)
           connect(path.qdd, pathToAxis5.qdd)
           connect(path.moving, pathToAxis5.moving)
           connect(path.q, pathToAxis6.q)
           connect(path.qd, pathToAxis6.qd)
           connect(path.qdd, pathToAxis6.qdd)
           connect(path.moving, pathToAxis6.moving)
           connect(pathToAxis1.axisControlBus, controlBus.axisControlBus1)
           connect(pathToAxis2.axisControlBus, controlBus.axisControlBus2)
           connect(pathToAxis3.axisControlBus, controlBus.axisControlBus3)
           connect(pathToAxis4.axisControlBus, controlBus.axisControlBus4)
           connect(pathToAxis5.axisControlBus, controlBus.axisControlBus5)
           connect(pathToAxis6.axisControlBus, controlBus.axisControlBus6)]
end

"Map path planning to one axis control bus"
function PathToAxisControlBus(; name, nAxis = 6, axisUsed = 1)
    @parameters begin
        nAxis = nAxis, [description = "Number of driven axis"]
        axisUsed = axisUsed,
                   [description = "Map path planning of axisUsed to axisControlBus"]
    end

    systems = @named begin
        q = Blocks.RealInput(nAxis)
        qd = Blocks.RealInput(nAxis)
        qdd = Blocks.RealInput(nAxis)
        axisControlBus = AxisControlBus()
        q_axisUsed = Blocks.RealPassThrough()
        qd_axisUsed = Blocks.RealPassThrough()
        qdd_axisUsed = Blocks.RealPassThrough()
        moving = Blocks.BooleanInput(nAxis)
        motion_ref_axisUsed = Blocks.BooleanPassThrough()
    end

    eqs = [connect(q_axisUsed.input, q[axisUsed])
           connect(qd_axisUsed.input, qd[axisUsed])
           connect(qdd_axisUsed.input, qdd[axisUsed])
           connect(motion_ref_axisUsed.input, moving[axisUsed])
           connect(motion_ref_axisUsed.output, axisControlBus.motion_ref)
           connect(qdd_axisUsed.output, axisControlBus.acceleration_ref)
           connect(qd_axisUsed.output, axisControlBus.speed_ref)
           connect(q_axisUsed.output, axisControlBus.angle_ref)]
    ODESystem(eqs, t; systems, name)
end
