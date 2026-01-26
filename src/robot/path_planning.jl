using DataInterpolations
using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput
using TrajectoryLimiters: JerkLimiter, calculate_trajectory, duration, evaluate_at

include("ptp.jl")

"Generate reference angles for specified kinematic movement"
function PathPlanning1(; name, q0deg = 0, q1deg = 1, time = 0:0.01:10, speed_max=1, acc_max=1, kwargs...)
    systems = @named begin
        controlBus = ControlBus()
        path = KinematicPTP(; q1 = deg2rad.(q1deg),
                            time,
                             qd_max = (speed_max),
                             qdd_max = (acc_max),
                            #  startTime = startTime,
                            q0 = deg2rad.(q0deg), kwargs...)
        pathToAxis1 = PathToAxisControlBus(nAxis = 1, axisUsed = 1)
    end

    eqs = [connect(path.q, pathToAxis1.q)
           connect(path.qd, pathToAxis1.qd)
           connect(path.qdd, pathToAxis1.qdd)
           #    connect(path.moving, pathToAxis1.moving)
           connect(pathToAxis1.axisControlBus, controlBus.axisControlBus1)]
    System(eqs, t; name, systems)
end

function PathPlanning6(; name, naxis = 6, q0deg = zeros(naxis),
                       q1deg = ones(naxis), time = 0:0.01:4,
                       speed_max = fill(3, naxis),
                       acc_max = fill(2.5, naxis), kwargs...)

    systems = @named begin
        controlBus = ControlBus()
        path = KinematicPTP(; q1 = deg2rad.(q1deg),
                            time,
                             qd_max = speed_max,
                             qdd_max = acc_max,
                            q0 = deg2rad.(q0deg), kwargs...)
        pathToAxis1 = PathToAxisControlBus(nAxis = naxis, axisUsed = 1)
        pathToAxis2 = PathToAxisControlBus(nAxis = naxis, axisUsed = 2)
        pathToAxis3 = PathToAxisControlBus(nAxis = naxis, axisUsed = 3)
        pathToAxis4 = PathToAxisControlBus(nAxis = naxis, axisUsed = 4)
        pathToAxis5 = PathToAxisControlBus(nAxis = naxis, axisUsed = 5)
        pathToAxis6 = PathToAxisControlBus(nAxis = naxis, axisUsed = 6)
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

    System(eqs, t; name, systems)
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
        moving = Blocks.Constant(k = 1) # Blocks.BooleanInput(nAxis) 
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
    System(eqs, t; systems, name)
end

"""
    q, qd, qdd = traj5(t; q0, q1, q̇0 = zero(q0), q̇1 = zero(q0), q̈0 = zero(q0), q̈1 = zero(q0))

Generate a 5:th order polynomial trajectory with specified end points, vels and accs.

See also [`point_to_point`](@ref) and [`Kinematic5`](@ref).
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
    KinematicPTP(; time, name, q0 = 0, q1 = 1, qd_max=1, qdd_max=1)

A component emitting a trajectory created by the [`point_to_point`](@ref) trajectory generator.

# Arguments
- `time`: Time vector, e.g., `0:0.01:10`
- `name`: Name of the component
- `q0`: Initial position
- `q1`: Final position
- `qd_max`: Maximum velocity
- `qdd_max`: Maximum acceleration

# Outputs
- `q`: Position
- `qd`: Velocity
- `qdd`: Acceleration

See also [`Kinematic5`](@ref).
"""
function KinematicPTP(; time, name, q0 = 0, q1 = 1, qd_max=1, qdd_max=1)
    nout = max(length(q0), length(q1))

    systems = @named begin
        q = RealOutput(; nout)
        qd = RealOutput(; nout)
        qdd = RealOutput(; nout)
    end

    q_vec, qd_vec, qdd_vec = point_to_point(time; q0 = q0, q1 = q1, qd_max, qdd_max)

    interp_eqs = map(1:nout) do i
        qfun = CubicSpline(q_vec[:, i], time; extrapolation=ExtrapolationType.Constant)
        qdfun = LinearInterpolation(qd_vec[:, i], time; extrapolation=ExtrapolationType.Constant)
        qddfun = ConstantInterpolation(qdd_vec[:, i], time; extrapolation=ExtrapolationType.Constant)
        [q.u[i] ~ qfun(t) 
        qd.u[i] ~ qdfun(t)
        qdd.u[i] ~ qddfun(t)]
    end
    eqs = reduce(vcat, interp_eqs)
    System(eqs, t; name, systems)
end

"""
    KinematicPTPBoundedJerk(; name, q0 = 0, q1 = 1, qd_max=1, qdd_max=1, qddd_max=10)

A component emitting a time-optimal point-to-point trajectory with bounded velocity,
acceleration, and jerk, generated using `JerkLimiter` from TrajectoryLimiters.jl.

When multiple axes are specified, the trajectories are time-synchronized so all axes
reach their targets at the same time.

# Arguments
- `name`: Name of the component
- `q0`: Initial position (scalar or vector)
- `q1`: Final position (scalar or vector)
- `qd_max`: Maximum velocity (scalar or vector)
- `qdd_max`: Maximum acceleration (scalar or vector)
- `qddd_max`: Maximum jerk (scalar or vector)

# Outputs
- `q`: Position
- `qd`: Velocity
- `qdd`: Acceleration
- `qddd`: Jerk

See also [`KinematicPTP`](@ref) and [`Kinematic5`](@ref).
"""
function KinematicPTPBoundedJerk(; name, q0 = 0, q1 = 1, qd_max=1, qdd_max=1, qddd_max=10)
    nout = max(length(q0), length(q1))

    # Broadcast parameters to nout dimensions
    q0_vec = q0 isa Number ? fill(float(q0), nout) : collect(float.(q0))
    q1_vec = q1 isa Number ? fill(float(q1), nout) : collect(float.(q1))
    qd_max_vec = qd_max isa Number ? fill(float(qd_max), nout) : collect(float.(qd_max))
    qdd_max_vec = qdd_max isa Number ? fill(float(qdd_max), nout) : collect(float.(qdd_max))
    qddd_max_vec = qddd_max isa Number ? fill(float(qddd_max), nout) : collect(float.(qddd_max))

    systems = @named begin
        q = RealOutput(; nout)
        qd = RealOutput(; nout)
        qdd = RealOutput(; nout)
        qddd = RealOutput(; nout)
    end

    # Create vector of JerkLimiters for time-synchronized multi-DOF trajectory
    lims = [JerkLimiter(; vmax=qd_max_vec[i], amax=qdd_max_vec[i], jmax=qddd_max_vec[i]) for i in 1:nout]

    # Calculate time-synchronized trajectories for all axes
    profiles = calculate_trajectory(lims; p0=q0_vec, pf=q1_vec)
    [q.u[i] ~ evaluate_at_1(profiles, t)
        qd.u[i] ~ evaluate_at_2(profiles, t)
        qdd.u[i] ~ evaluate_at_3(profiles, t)
        qddd.u[i] ~ evaluate_at_4(profiles, t)]

    System(eqs, t; name, systems)
end

evaluate_at_1(profile, ti::Real) = evaluate_at(profile, ti::Real)[1]
evaluate_at_2(profile, ti::Real) = evaluate_at(profile, ti::Real)[2]
evaluate_at_3(profile, ti::Real) = evaluate_at(profile, ti::Real)[3]
evaluate_at_4(profile, ti::Real) = evaluate_at(profile, ti::Real)[4]

@register_symbolic evaluate_at_1(profile, ti::Real)
@register_symbolic evaluate_at_2(profile, ti::Real)
@register_symbolic evaluate_at_3(profile, ti::Real)
@register_symbolic evaluate_at_4(profile, ti::Real)

"""
    Kinematic5(; time, name, q0 = 0, q1 = 1, qd0 = 0, qd1 = 0, qdd0 = 0, qdd1 = 0)

A component emitting a 5:th order polynomial trajectory created using [`traj5`](@ref). `traj5` is a simple trajectory planner that plans a 5:th order polynomial trajectory between two points, subject to specified boundary conditions on the position, velocity and acceleration.

# Arguments
- `time`: Time vector, e.g., `0:0.01:10`
- `name`: Name of the component

# Outputs
- `q`: Position
- `qd`: Velocity
- `qdd`: Acceleration
"""
function Kinematic5(; time, name, q0 = 0, q1 = 1, qd0 = 0, qd1 = 0,
                      qdd0 = 0, qdd1 = 0)
    nout = max(length(q0), length(q1), length(qd1), length(qdd1))

    systems = @named begin
        q = RealOutput(; nout)
        qd = RealOutput(; nout)
        qdd = RealOutput(; nout)
    end

    interp_eqs = map(1:nout) do i
        _q, _qd, _qdd = traj5(t, time[end]; q0 = q0[i], q1 = q1[i],
                                    q̇0 = zero(q0[i]),
                                    q̇1 = zero(q0[i]),
                                    q̈0 = zero(q0[i]),
                                    q̈1 = zero(q0[i]))
        [q.u[i] ~ _q 
        qd.u[i] ~ _qd
        qdd.u[i] ~ _qdd]

    end
    eqs = reduce(vcat, interp_eqs)
    System(eqs, t; name, systems)
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
    extend(System(eqs, t; name), siso)
end
