using DataInterpolations
using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput

include("ptp.jl")

"Generate reference angles for specified kinematic movement"
function PathPlanning1(; name, angleBegDeg = 0, angleEndDeg = 1, time = 0:0.01:10,
                       swingTime = 0.5, kwargs...)
    # @parameters begin
    # angleBegDeg = angleBegDeg, [description = "Start angle"]
    # angleEndDeg = angleEndDeg, [description = "End angle"]
    # speedMax = speedMax, [description = "Maximum axis speed"]
    # accMax = accMax, [description = "Maximum axis acceleration"]
    # startTime = startTime, [description = "Start time of movement"]
    # swingTime = swingTime,
    # [
    #     description = "Additional time after reference motion is in rest before simulation is stopped",
    # ]
    # angleBeg = deg2rad(angleBegDeg), [description = "Start angles"]
    # angleEnd = deg2rad(angleEndDeg), [description = "End angles"]
    # end

    systems = @named begin
        controlBus = ControlBus()
        path = KinematicPTP(; q_end = deg2rad.(angleEndDeg),
                            time,
                            #  qd_max = speedMax,
                            #  qdd_max = accMax,
                            #  startTime = startTime,
                            q_begin = deg2rad.(angleBegDeg), kwargs...)
        pathToAxis1 = PathToAxisControlBus(nAxis = 1, axisUsed = 1)
        # terminateSimulation = TerminateSimulation(condition = time >=
        #                                                       path.endTime + swingTime)
    end

    eqs = [connect(path.q, pathToAxis1.q)
           connect(path.qd, pathToAxis1.qd)
           connect(path.qdd, pathToAxis1.qdd)
           #    connect(path.moving, pathToAxis1.moving)
           connect(pathToAxis1.axisControlBus, controlBus.axisControlBus1)]
    ODESystem(eqs, t; name, systems)
end

function PathPlanning6(; name, naxis = 6, angleBegDeg = zeros(naxis),
                       angleEndDeg = ones(naxis), time = 0:0.01:4,
                       speedMax = fill(3, naxis),
                       accMax = fill(2.5, naxis), startTime = 0, swingTime = 0.5, kwargs...)
    # @parameters begin
    #     naxis = naxis, [description = "Number of driven axis"]
    #     angleBegDeg[1:naxis] = angleBegDeg, [description = "Start angles"]
    #     angleEndDeg[1:naxis] = angleEndDeg, [description = "End angles"]
    #     # speedMax[1:naxis] = speedMax, [description = "Maximum axis speed"]
    #     # accMax[1:naxis] = accMax, [description = "Maximum axis acceleration"]
    #     # startTime = startTime, [description = "Start time of movement"]
    #     swingTime = swingTime,
    #                 [
    #                     description = "Additional time after reference motion is in rest before simulation is stopped",
    #                 ]
    #     angleBeg[1:6] = deg2rad.(angleBegDeg), [description = "Start angles"]
    #     angleEnd[1:6] = deg2rad.(angleEndDeg), [description = "End angles"]
    # end

    systems = @named begin
        controlBus = ControlBus()
        path = KinematicPTP(; q_end = deg2rad.(angleEndDeg),
                            time,
                             qd_max = speedMax,
                             qdd_max = accMax,
                            #  startTime = startTime,
                            q_begin = deg2rad.(angleBegDeg), kwargs...)
        pathToAxis1 = PathToAxisControlBus(nAxis = naxis, axisUsed = 1)
        pathToAxis2 = PathToAxisControlBus(nAxis = naxis, axisUsed = 2)
        pathToAxis3 = PathToAxisControlBus(nAxis = naxis, axisUsed = 3)
        pathToAxis4 = PathToAxisControlBus(nAxis = naxis, axisUsed = 4)
        pathToAxis5 = PathToAxisControlBus(nAxis = naxis, axisUsed = 5)
        pathToAxis6 = PathToAxisControlBus(nAxis = naxis, axisUsed = 6)
        # terminateSimulation = TerminateSimulation(condition = time >=
        #                                                       path.endTime + swingTime)
    end

    eqs = [connect(path.q, pathToAxis1.q)
           connect(path.qd, pathToAxis1.qd)
           connect(path.qdd, pathToAxis1.qdd)
           #    connect(path.moving, pathToAxis1.moving)
           connect(path.q, pathToAxis2.q)
           connect(path.qd, pathToAxis2.qd)
           connect(path.qdd, pathToAxis2.qdd)
           #    connect(path.moving, pathToAxis2.moving)
           connect(path.q, pathToAxis3.q)
           connect(path.qd, pathToAxis3.qd)
           connect(path.qdd, pathToAxis3.qdd)
           #    connect(path.moving, pathToAxis3.moving)
           connect(path.q, pathToAxis4.q)
           connect(path.qd, pathToAxis4.qd)
           connect(path.qdd, pathToAxis4.qdd)
           #    connect(path.moving, pathToAxis4.moving)
           connect(path.q, pathToAxis5.q)
           connect(path.qd, pathToAxis5.qd)
           connect(path.qdd, pathToAxis5.qdd)
           #    connect(path.moving, pathToAxis5.moving)
           connect(path.q, pathToAxis6.q)
           connect(path.qd, pathToAxis6.qd)
           connect(path.qdd, pathToAxis6.qdd)
           #    connect(path.moving, pathToAxis6.moving)
           connect(pathToAxis1.axisControlBus, controlBus.axisControlBus1)
           connect(pathToAxis2.axisControlBus, controlBus.axisControlBus2)
           connect(pathToAxis3.axisControlBus, controlBus.axisControlBus3)
           connect(pathToAxis4.axisControlBus, controlBus.axisControlBus4)
           connect(pathToAxis5.axisControlBus, controlBus.axisControlBus5)
           connect(pathToAxis6.axisControlBus, controlBus.axisControlBus6)]

    ODESystem(eqs, t; name, systems)
end

"Map path planning to one axis control bus"
function PathToAxisControlBus(; name, nAxis = 6, axisUsed = 1)
    # @parameters begin
    #     nAxis = nAxis, [description = "Number of driven axis"]
    #     axisUsed = axisUsed,
    #                [description = "Map path planning of axisUsed to axisControlBus"]
    # end

    systems = @named begin
        q = Blocks.RealInput(nin = nAxis)
        qd = Blocks.RealInput(nin = nAxis)
        qdd = Blocks.RealInput(nin = nAxis)
        axisControlBus = AxisControlBus()
        q_axisUsed = RealPassThrough()
        qd_axisUsed = RealPassThrough()
        qdd_axisUsed = RealPassThrough()
        moving = Blocks.Constant(k = 1) # Blocks.BooleanInput(nAxis) # NOTE
        motion_ref_axisUsed = RealPassThrough() # Blocks.BooleanPassThrough()
    end

    eqs = [(q_axisUsed.input.u ~ q.u[axisUsed])
           (qd_axisUsed.input.u ~ qd.u[axisUsed])
           (qdd_axisUsed.input.u ~ qdd.u[axisUsed])
           connect(motion_ref_axisUsed.input, moving.output)
           (motion_ref_axisUsed.output.u ~ axisControlBus.motion_ref)
           (qdd_axisUsed.output.u ~ axisControlBus.acceleration_ref)
           (qd_axisUsed.output.u ~ axisControlBus.speed_ref)
           (q_axisUsed.output.u ~ axisControlBus.angle_ref)]
    ODESystem(eqs, t; systems, name)
end

"""
    q, qd, qdd = traj5(t; q0, q1, q̇0 = zero(q0), q̇1 = zero(q0), q̈0 = zero(q0), q̈1 = zero(q0))

Generate a 5:th order polynomial trajectory with specified end points, vels and accs.
"""
function traj5(t; q0 = 0.0, q1 = one(q0), q̇0 = zero(q0), q̇1 = zero(q0), q̈0 = zero(q0),
               q̈1 = zero(q0))
    t[1] == 0 || throw(ArgumentError("t must start at 0"))
    tf = t[end]
    a0 = q0
    a1 = q̇0
    a2 = @. q̈0 / 2
    a3 = @. (20q1 - 20q0 - (8q̇1 + 12q̇0) * tf - (3q̈0 - q̈1) * tf^2) / 2tf^3
    a4 = @. (30q0 - 30q1 + (14q̇1 + 16q̇0) * tf + (3q̈0 - 2q̈1) * tf^2) / 2tf^4
    a5 = @. (12q1 - 12q0 - (6q̇1 + 6q̇0) * tf - (q̈0 - q̈1) * tf^2) / 2tf^5
    evalpoly.(t, ((a0, a1, a2, a3, a4, a5),)),
    evalpoly.(t, ((a1, 2a2, 3a3, 4a4, 5a5),)),
    evalpoly.(t, ((2a2, 6a3, 12a4, 20a5),))
end


function traj5(t, tf; q0 = 0.0, q1 = one(q0), q̇0 = zero(q0), q̇1 = zero(q0), q̈0 = zero(q0),
    q̈1 = zero(q0))
    a0 = q0
    a1 = q̇0
    a2 = q̈0 / 2
    a3 = (20q1 - 20q0 - (8q̇1 + 12q̇0) * tf - (3q̈0 - q̈1) * tf^2) / 2tf^3
    a4 = (30q0 - 30q1 + (14q̇1 + 16q̇0) * tf + (3q̈0 - 2q̈1) * tf^2) / 2tf^4
    a5 = (12q1 - 12q0 - (6q̇1 + 6q̇0) * tf - (q̈0 - q̈1) * tf^2) / 2tf^5
    evalpoly(t, (a0, a1, a2, a3, a4, a5)),
    evalpoly(t, (a1, 2a2, 3a3, 4a4, 5a5)),
    evalpoly(t, (2a2, 6a3, 12a4, 20a5))
end

# t = 0:100;
# q, qd, qdd = traj5(t, q0 = 1, q1 = 2, q̇0 = 3, q̇1 = 4, q̈0 = 5, q̈1 = 6);
# plot(t, q, label = "q", layout = (3, 1), sp = 1)
# plot!(t, qd, label = "qd", sp = 2)
# plot!(t, qdd, label = "qdd", sp = 3)
# plot!(t, centraldiff(q), sp = 2)
# plot!(t, centraldiff(centraldiff(q)), sp = 3)

"""
    KinematicPTP(; time, name, q_begin = 0, q_end = 1, qd_begin = 0, qd_end = 0, qdd_begin = 0, qdd_end = 0)

A simple trajectory planner that plans a 5:th order polynomial trajectory between two points, subject to specified boundary conditions on the position, velocity and acceleration.
"""
function KinematicPTP(; time, name, q_begin = 0, q_end = 1, qd_begin = 0, qd_end = 0,
                      qdd_begin = 0, qdd_end = 0, trivial = false, qd_max=1, qdd_max=1)
    nout = max(length(q_begin), length(q_end), length(qd_end), length(qdd_end))

    # @parameters begin
    #     q_begin = q_begin, [description = "Start position"]
    #     q_end = q_end, [description = "End position"]
    #     qd_max = qd_max, [description = "Maximum velocities der(q)"]
    #     qdd_max = qdd_max, [description = "Maximum accelerations der(qd)"]
    #     startTime = startTime, [description = "Time instant at which movement starts"]
    #     p_q_begin[1:nout] = ones(nout) .* q_begin
    #     p_q_end[1:nout] = ones(nout) .* q_end
    #     p_qd_max[1:nout] = ones(nout) .* qd_max
    #     p_qdd_max[1:nout] = ones(nout) .* qdd_max
    #     p_deltaq[1:nout] = p_q_end - p_q_begin
    # end



    systems = @named begin
        q = RealOutput(; nout)
        qd = RealOutput(; nout)
        qdd = RealOutput(; nout)
        # moving = BooleanOutput(; nout)
    end


    startTime = time[1]
    time0 = time .- startTime # traj5 wants time vector to start at 0
    if !trivial
        q_vec, qd_vec, qdd_vec = PTP(time; q0 = q_begin, q1 = q_end, qd_max, qdd_max)
    end

    interp_eqs = map(1:nout) do i
        if trivial
            _q, _qd, _qdd = traj5(t, time[end]; q0 = q_begin[i], q1 = q_end[i],
                                       q̇0 = zero(q_begin[i]),
                                       q̇1 = zero(q_begin[i]),
                                       q̈0 = zero(q_begin[i]),
                                       q̈1 = zero(q_begin[i]))
            [q.u[i] ~ _q 
            qd.u[i] ~ _qd
            qdd.u[i] ~ _qdd]
        else
            # q_vec, qd_vec, qdd_vec = PTP(time; q0 = q_begin[i], q1 = q_end[i], qd_max, qdd_max)
            qfun = CubicSpline(q_vec[:, i], time; extrapolate=true)
            qdfun = LinearInterpolation(qd_vec[:, i], time; extrapolate=true)
            qddfun = ConstantInterpolation(qdd_vec[:, i], time; extrapolate=true)
            [q.u[i] ~ qfun(t) 
            qd.u[i] ~ qdfun(t)
            qdd.u[i] ~ qddfun(t)]
        end
    end
    eqs = reduce(vcat, interp_eqs)

    # push!(eqs, moving.u ~ (time[1] < t < time[end]))

    ODESystem(eqs, t; name, systems)
end

"""
    RealPassThrough(; name)

Pass a Real signal through without modification

# Connectors
- `input`
- `output`
"""
@component function RealPassThrough(; name)
    @named siso = Blocks.SISO()
    @unpack u, y = siso
    eqs = [y ~ u]
    extend(ODESystem(eqs, t; name), siso)
end
