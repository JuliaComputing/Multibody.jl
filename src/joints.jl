
"""
    Revolute(; name, phi0 = 0, w0 = 0, n, useAxisFlange = false)

Revolute joint with 1 rotational degree-of-freedom

- `phi0`: Initial angle
- `w0`: Iniitial angular velocity
- `n`: The axis of rotation
- `useAxisFlange`: If true, the joint will have two additional frames from Mechanical.Rotational, `axis` and `support`, between which rotational components such as springs and dampers can be connected.

If `useAxisFlange`, flange connectors for ModelicaStandardLibrary.Mechanics.Rotational are also available:
- `axis`: 1-dim. rotational flange that drives the joint
- `support`: 1-dim. rotational flange of the drive support (assumed to be fixed in the world frame, NOT in the joint)
"""
function Revolute(; name, phi0 = 0, w0 = 0, n = Float64[0, 0, 1], useAxisFlange = false,
                  isroot = false)
    norm(n) ≈ 1 || error("Axis of rotation must be a unit vector")
    @named frame_a = Frame()
    @named frame_b = Frame()
    @parameters n[1:3]=n [description = "axis of rotation"]
    @variables tau(t)=0 [
        connect = Flow,
        state_priority = 2,
        description = "Driving torque in direction of axis of rotation",
    ]
    @variables phi(t)=phi0 [
        state_priority = 20,
        description = "Relative rotation angle from frame_a to frame_b",
    ]
    @variables w(t)=w0 [state_priority = 20, description = "angular velocity (rad/s)"]
    Rrel0 = planar_rotation(n, phi0, w0)
    @named Rrel = NumRotationMatrix(; R = Rrel0.R, w = Rrel0.w)
    n = collect(n)

    if isroot
        eqs = Equation[Rrel ~ planar_rotation(n, phi, w)
                       ori(frame_b) ~ absoluteRotation(ori(frame_a), Rrel)
                       collect(frame_a.f) .~ -resolve1(Rrel, frame_b.f)
                       collect(frame_a.tau) .~ -resolve1(Rrel, frame_b.tau)]
    else
        eqs = Equation[Rrel ~ planar_rotation(-n, phi, w)
                       ori(frame_a) ~ absoluteRotation(ori(frame_b), Rrel)
                       collect(frame_b.f) .~ -resolve1(Rrel, frame_a.f)
                       collect(frame_b.tau) .~ -resolve1(Rrel, frame_a.tau)]
    end
    moreeqs = [collect(frame_a.r_0 .~ frame_b.r_0)
               D(phi) ~ w
               tau ~ -collect(frame_b.tau)'n]
    append!(eqs, moreeqs)
    if useAxisFlange
        # @named internalAxis = Rotational.InternalSupport(tau=tau)
        @named fixed = Rotational.Fixed()

        @named axis = Rotational.Flange()
        @named support = Rotational.Flange()
        # push!(eqs, phi ~ internalAxis.phi)
        push!(eqs, connect(fixed.flange, support))
        push!(eqs, axis.phi ~ phi)
        push!(eqs, axis.tau ~ tau)
        # push!(eqs, connect(internalAxis.flange, axis))
        compose(ODESystem(eqs, t; name), frame_a, frame_b, axis, support, fixed)
    else
        # Modelica Revolute uses a ConstantTorque as well as internalAxis = Rotational.InternalSupport(tau=tau), but it seemed more complicated than required and I couldn't get it to work, likely due to the `input` semantics of modelica not having an equivalent in MTK, so the (tau=tau) input argument caused problems.
        # @named constantTorque = Rotational.ConstantTorque(tau_constant=0, use_support=false) 
        # push!(eqs, connect(constantTorque.flange, internalAxis.flange))
        push!(eqs, tau ~ 0)
        compose(ODESystem(eqs, t; name), frame_a, frame_b)
    end
end

"""
    Prismatic(; name, n = [0, 0, 1], useAxisFlange = false, isroot = true)

Prismatic joint with 1 translational degree-of-freedom

- `n`: The axis of motion (unit vector)
- `useAxisFlange`: If true, the joint will have two additional frames from Mechanical.Translational, `axis` and `support`, between which translational components such as springs and dampers can be connected.
- `isroot`: If true, the joint will be considered the root of the system.

If `useAxisFlange`, flange connectors for ModelicaStandardLibrary.Mechanics.TranslationalModelica are also available:
- `axis`: 1-dim. translational flange that drives the joint
- `support`: 1-dim. translational flange of the drive support (assumed to be fixed in the world frame, NOT in the joint)

The function returns an ODESystem representing the prismatic joint.
"""
function Prismatic(; name, n = Float64[0, 0, 1], useAxisFlange = false,
                   isroot = true, s0 = 0, v0 = 0)
    norm(n) ≈ 1 || error("Axis of motion must be a unit vector")
    @named frame_a = Frame()
    @named frame_b = Frame()
    @parameters n[1:3]=n [description = "axis of motion"]
    n = collect(n)

    @variables s(t)=s0 [
        state_priority = 10,
        description = "Relative distance between frame_a and frame_b",
    ]
    @variables v(t)=v0 [
        state_priority = 10,
        description = "Relative velocity between frame_a and frame_b",
    ]
    @variables a(t)=0 [
        state_priority = 10,
        description = "Relative acceleration between frame_a and frame_b",
    ]
    @variables f(t)=0 [
        connect = Flow,
        description = "Actuation force in direction of joint axis",
    ]

    eqs = [v ~ D(s)
           a ~ D(v)

           # relationships between kinematic quantities of frame_a and of frame_b
           collect(frame_b.r_0) .~ collect(frame_a.r_0) + resolve1(ori(frame_a), n * s)
           ori(frame_b) ~ ori(frame_a)

           # Force and torque balance
           zeros(3) .~ collect(frame_a.f + frame_b.f)
           zeros(3) .~ collect(frame_a.tau + frame_b.tau + cross(n * s, frame_b.f))

           # d'Alemberts principle
           f .~ -n'collect(frame_b.f)]

    if useAxisFlange
        @named fixed = Translational.Fixed(s0=0)
        @named axis = Translational.Flange()
        @named support = Translational.Flange()
        push!(eqs, connect(fixed.flange, support))
        push!(eqs, axis.s ~ s)
        push!(eqs, axis.f ~ f)
        compose(ODESystem(eqs, t; name), frame_a, frame_b, axis, support, fixed)
    else
        push!(eqs, f ~ 0)
        compose(ODESystem(eqs, t; name), frame_a, frame_b)
    end
end

"""
    Spherical(; name, enforceState = false, isroot = true, w_rel_a_fixed = false, z_rel_a_fixed = false, sequence_angleStates, phi = 0, phi_d = 0, phi_dd = 0)

Joint with 3 constraints that define that the origin of `frame_a` and the origin of `frame_b` coincide. By default this joint defines only the 3 constraints without any potential state variables. If parameter `enforceState` is set to true, three states are introduced. The orientation of `frame_b` is computed by rotating `frame_a` along the axes defined in parameter vector `sequence_angleStates` (default = [1,2,3], i.e., the Cardan angle sequence) around the angles used as state. If angles are used as state there is the slight disadvantage that a singular configuration is present leading to a division by zero.

- `isroot`: Indicate that `frame_a` is the root, otherwise `frame_b` is the root. Only relevant if `enforceState = true`.
- `sequence_angleStates`: Rotation sequence
"""
function Spherical(; name, enforceState = false, isroot = true, w_rel_a_fixed = false,
                   z_rel_a_fixed = false, sequence_angleStates = [1, 2, 3], phi = 0,
                   phi_d = 0,
                   phi_dd = 0)
    @named begin
        ptf = PartialTwoFrames()
        R_rel = NumRotationMatrix()
        R_rel_inv = NumRotationMatrix()
    end
    @unpack frame_a, frame_b = ptf
    # @parameters begin # Currently not using parameters due to these appearing in if statements
    #     sequence_angleStates[1:3] = sequence_angleStates
    # end
    @variables begin (w_rel(t)[1:3] = zeros(3)),
                     [
                         description = "relative angular velocity of frame_b with respect to frame_a, resolved in frame_b",
                     ] end

    # torque balance
    eqs = [zeros(3) .~ collect(frame_a.tau)
           zeros(3) .~ collect(frame_b.tau)
           collect(frame_b.r_0) .~ collect(frame_a.r_0)]

    if enforceState
        @variables begin
            (phi(t)[1:3] = phi),
            [state_priority = 10, description = "3 angles to rotate frame_a into frame_b"]
            (phi_d(t)[1:3] = phi_d),
            [state_priority = 10, description = "3 angle derivatives"]
            (phi_dd(t)[1:3] = phi_dd),
            [state_priority = 10, description = "3 angle second derivatives"]
        end
        append!(eqs,
                [R_rel ~ axesRotations(sequence_angleStates, phi, phi_d)
                 collect(w_rel) .~ angularVelocity2(R_rel)
                 collect(phi_d .~ D.(phi))
                 collect(phi_dd .~ D.(phi_d))])
        if isroot
            append!(eqs,
                    [ori(frame_b) ~ absoluteRotation(frame_a, R_rel)
                     zeros(3) .~ collect(frame_a.f) + resolve1(R_rel, frame_b.f)])
        else
            append!(eqs,
                    [R_rel_inv ~ inverseRotation(R_rel)
                     ori(frame_a) ~ absoluteRotation(frame_b, R_rel_inv)
                     zeros(3) .~ collect(frame_b.f) + resolve2(R_rel, frame_a.f)])
        end

    else
        # Spherical joint does not have state
        append!(eqs,
                #frame_b.r_0 ~ transpose(frame_b.R.T)*(frame_b.R.T*(transpose(frame_a.R.T)*(frame_a.R.T*frame_a.r_0)));
                zeros(3) .~ collect(frame_a.f) +
                            resolveRelative(frame_b.f, frame_b, frame_a))
        if w_rel_a_fixed || z_rel_a_fixed
            append!(w_rel .~ angularVelocity2(frame_b) - resolve2(frame_b,
                                      angularVelocity1(frame_a)))
        else
            append!(w_rel .~ zeros(3))
        end
    end

    extend(ODESystem(eqs, t; name), ptf)
end

function Universal(; name, n_a = [1, 0, 0], n_b = [0, 1, 0], phi_a = 0,
                   phi_b = 0,
                   w_a = 0,
                   w_b = 0,
                   a_a = 0,
                   a_b = 0)
    @named begin
        ptf = PartialTwoFrames()
        revolute_a = Revolute(n = n_a, isroot = false)
        revolute_b = Revolute(n = n_b, isroot = false)
    end
    @unpack frame_a, frame_b = ptf
    @parameters begin
        n_a = n_a,
              [
                  description = "Axis of revolute joint 1 resolved in frame_a",
              ]
        n_b = n_b,
              [
                  description = "Axis of revolute joint 2 resolved in frame_b",
              ]
    end
    @variables begin
        (phi_a(t) = phi_a),
        [
            state_priority = 10,
            description = "Relative rotation angle from frame_a to intermediate frame",
        ]
        (phi_b(t) = phi_b),
        [
            state_priority = 10,
            description = "Relative rotation angle from intermediate frame to frame_b",
        ]
        (w_a(t) = w_a),
        [
            state_priority = 10,
            description = "First derivative of angle phi_a (relative angular velocity a)",
        ]
        (w_b(t) = w_b),
        [
            state_priority = 10,
            description = "First derivative of angle phi_b (relative angular velocity b)",
        ]
        (a_a(t) = a_a),
        [
            state_priority = 10,
            description = "Second derivative of angle phi_a (relative angular acceleration a)",
        ]
        (a_b(t) = a_b),
        [
            state_priority = 10,
            description = "Second derivative of angle phi_b (relative angular acceleration b)",
        ]
    end
    eqs = [phi_a ~ revolute_a.phi
           phi_b ~ revolute_b.phi
           w_a ~ D(phi_a)
           w_b ~ D(phi_b)
           a_a ~ D(w_a)
           a_b ~ D(w_b)
           connect(frame_a, revolute_a.frame_a)
           connect(revolute_b.frame_b, frame_b)
           connect(revolute_a.frame_b, revolute_b.frame_a)]
    extend(ODESystem(eqs, t; name, systems = [revolute_a, revolute_b]), ptf)
end

"""
    GearConstraint(; name, ratio, checkTotalPower = false, n_a, n_b, r_a, r_b)

This ideal massless joint provides a gear constraint between frames `frame_a` and `frame_b`. The axes of rotation of `frame_a` and `frame_b` may be arbitrary.

- `ratio`: Gear ratio
- `n_a`: Axis of rotation of `frame_a`
- `n_b`: Axis of rotation of `frame_b`
- `r_a`: Vector from frame `bearing` to `frame_a` resolved in bearing
- `r_b`: Vector from frame `bearing` to `frame_b` resolved in bearing
"""
function GearConstraint(; name, ratio, checkTotalPower = false, n_a = [1, 0, 0],
                        n_b = [1, 0, 0], r_a = [0, 0, 0], r_b = [0, 0, 0])
    @named ptf = PartialTwoFrames()
    systems = @named begin
        bearing = Frame() #"Coordinate system fixed in the bearing"

        actuatedRevolute_a = Revolute(useAxisFlange = true,
                                      n = n_a)
        actuatedRevolute_b = Revolute(useAxisFlange = true,
                                      n = n_b)

        idealGear = Rotational.IdealGear(ratio = ratio)
        fixedTranslation1 = FixedTranslation(r = r_b)
        fixedTranslation2 = FixedTranslation(r = r_a)
    end
    @unpack frame_a, frame_b = ptf

    @parameters begin
        ratio = ratio, [description = "Gear speed ratio"]

        n_a = n_a,
              [
                  description = "Axis of rotation of shaft a (same coordinates in frame_a, frame_b, bearing)",
              ]
        n_b = n_b,
              [
                  description = "Axis of rotation of shaft b (same coordinates in frame_a, frame_b, bearing)",
              ]

        r_a[1:3] = r_a,
                   [
                       description = "Vector from frame bearing to frame_a resolved in bearing",
                   ]
        r_b[1:3] = r_b,
                   [
                       description = "Vector from frame bearing to frame_b resolved in bearing",
                   ]
    end

    @variables begin
        (phi_b(t) = 0),
        [
            state_priority = 10,
            description = "Relative rotation angle of revolute joint at frame_b",
        ]

        (w_b(t) = 0),
        [
            state_priority = 10,
            description = "Relative angular velocity of revolute joint at frame_b",
        ]
        (a_b(t) = 0),
        [
            state_priority = 10,
            description = "Relative angular acceleration of revolute joint at frame_b",
        ]

        (totalPower(t) = 0), [description = "Total power flowing into this element"]
    end

    eqs = [phi_b ~ actuatedRevolute_b.phi
           w_b ~ D(phi_b)
           a_b ~ D(w_b)
           connect(actuatedRevolute_a.axis, idealGear.flange_a)
           connect(idealGear.flange_b, actuatedRevolute_b.axis)
           connect(actuatedRevolute_a.frame_a, fixedTranslation2.frame_b)
           connect(fixedTranslation2.frame_a, bearing)
           connect(fixedTranslation1.frame_a, bearing)
           connect(fixedTranslation1.frame_b, actuatedRevolute_b.frame_a)
           connect(frame_a, actuatedRevolute_a.frame_b)
           connect(actuatedRevolute_b.frame_b, frame_b)]

    # Measure power for test purposes
    if checkTotalPower
        push!(eqs,
              totalPower ~ frame_a.f'resolve2(frame_a, D.(frame_a.r_0)) +
                           frame_b.f'resolve2(frame_b, D.(frame_b.r_0)) +
                           bearing.f'resolve2(bearing, D.(bearing.r_0)) +
                           frame_a.tau'angularVelocity2(frame_a) +
                           frame_b.tau'angularVelocity2(frame_b) +
                           bearing.tau'angularVelocity2(bearing))
    end

    extend(ODESystem(eqs, t; name, systems), ptf)
end

"""
    RollingWheelJoint(; name, radius, angles, x0, y0, z0)

Joint (no mass, no inertia) that describes an ideal rolling wheel (rolling on the plane z=0). See [`RollingWheel`](@ref) for a realistic wheel model with inertia.

A joint for a wheel rolling on the x-y plane of the world frame.
The rolling contact is considered being ideal, i.e. there is no
slip between the wheel and the ground. This is simply
gained by two non-holonomic constraint equations on velocity level
defined for both longitudinal and lateral direction of the wheel.
There is also a holonomic constraint equation on position level
granting a permanent contact of the wheel to the ground, i.e.
the wheel can not take off.

The origin of the frame `frame_a` is placed in the intersection
of the wheel spin axis with the wheel middle plane and rotates
with the wheel itself. The y-axis of `frame_a` is identical with
the wheel spin axis, i.e. the wheel rotates about y-axis of `frame_a`.
A wheel body collecting the mass and inertia should be connected to
this frame.

# Arguments and parameters:

name: Name of the rolling wheel joint component
radius: Radius of the wheel
angles: Angles to rotate world-frame into frame_a around z-, y-, x-axis

# Variables:
- `x`: x-position of the wheel axis
- `y`: y-position of the wheel axis
- `z`: z-position of the wheel axis
- `angles`: Angles to rotate world-frame into `frame_a` around z-, y-, x-axis
- `der_angles`: Derivatives of angles
- `r_road_0`: Position vector from world frame to contact point on road, resolved in world frame
- `f_wheel_0`: Force vector on wheel, resolved in world frame
- `f_n`: Contact force acting on wheel in normal direction
- `f_lat`: Contact force acting on wheel in lateral direction
- `f_long`: Contact force acting on wheel in longitudinal direction
- `err`: Absolute value of `(r_road_0 - frame_a.r_0) - radius` (must be zero; used for checking)
- `e_axis_0`: Unit vector along wheel axis, resolved in world frame
- `delta_0`: Distance vector from wheel center to contact point
- `e_n_0`: Unit vector in normal direction of road at contact point, resolved in world frame
- `e_lat_0`: Unit vector in lateral direction of road at contact point, resolved in world frame
- `e_long_0`: Unit vector in longitudinal direction of road at contact point, resolved in world frame
- `s`: Road surface parameter 1
- `w`: Road surface parameter 2
- `e_s_0`: Road heading at `(s,w)`, resolved in world frame (unit vector)
- `v_0`: Velocity of wheel center, resolved in world frame
- `w_0`: Angular velocity of wheel, resolved in world frame
- `vContact_0`: Velocity of contact point, resolved in world frame

# Connector frames
- `frame_a`: Frame for the wheel joint
"""
function RollingWheelJoint(; name, radius, angles = zeros(3), x0, y0, z0 = 0)
    @named frame_a = Frame()
    @parameters begin radius = radius, [description = "Radius of the wheel"] end
    @variables begin
        (x(t) = x0), [state_priority = 10, description = "x-position of the wheel axis"]
        (y(t) = y0), [state_priority = 10, description = "y-position of the wheel axis"]
        (z(t) = z0), [state_priority = 10, description = "z-position of the wheel axis"]
        (angles(t)[1:3] = angles),
        [description = "Angles to rotate world-frame into frame_a around z-, y-, x-axis"]
        (der_angles(t)[1:3] = zeros(3)), [description = "Derivatives of angles"]
        (r_road_0(t)[1:3] = zeros(3)),
        [
            description = "Position vector from world frame to contact point on road, resolved in world frame",
        ]
        (f_wheel_0(t)[1:3] = zeros(3)),
        [description = "Force vector on wheel, resolved in world frame"]
        (f_n(t) = 0), [description = "Contact force acting on wheel in normal direction"]
        (f_lat(t) = 0), [
            description = "Contact force acting on wheel in lateral direction",
        ]
        (f_long(t) = 0),
        [description = "Contact force acting on wheel in longitudinal direction"]
        (err(t) = 0),
        [
            description = "|r_road_0 - frame_a.r_0| - radius (must be zero; used for checking)",
        ]
        (e_axis_0(t)[1:3] = zeros(3)),
        [description = "Unit vector along wheel axis, resolved in world frame"]
        (delta_0(t)[1:3] = zeros(3)),
        [description = "Distance vector from wheel center to contact point"]
        (e_n_0(t)[1:3] = zeros(3)),
        [
            description = "Unit vector in normal direction of road at contact point, resolved in world frame",
        ]
        (e_lat_0(t)[1:3] = zeros(3)),
        [
            description = "Unit vector in lateral direction of road at contact point, resolved in world frame",
        ]
        (e_long_0(t)[1:3] = zeros(3)),
        [
            description = "Unit vector in longitudinal direction of road at contact point, resolved in world frame",
        ]

        (s(t) = 0), [state_priority = 10, description = "Road surface parameter 1"]
        (w(t) = 0), [state_priority = 10, description = "Road surface parameter 2"]
        (e_s_0(t)[1:3] = zeros(3)),
        [description = "Road heading at (s,w), resolved in world frame (unit vector)"]

        (v_0(t)[1:3] = zeros(3)),
        [description = "Velocity of wheel center, resolved in world frame"]
        (w_0(t)[1:3] = zeros(3)),
        [description = "Angular velocity of wheel, resolved in world frame"]
        (vContact_0(t)[1:3] = zeros(3)),
        [description = "Velocity of contact point, resolved in world frame"]

        (aux(t)[1:3] = zeros(3)), [description = "Auxiliary variable"]
    end

    equations = [
                 # frame_a.R is computed from generalized coordinates
                 collect(frame_a.r_0) .~ [x, y, z]
                 collect(der_angles) .~ D.(angles)
                 ori(frame_a) ~ axesRotations([3, 2, 1], angles, der_angles)

                 # Road description
                 collect(r_road_0) .~ [s, w, 0]
                 collect(e_n_0) .~ [0, 0, 1]
                 collect(e_s_0) .~ [1, 0, 0]

                 # Coordinate system at contact point (e_long_0, e_lat_0, e_n_0)
                 collect(e_axis_0) .~ resolve1(ori(frame_a), [0, 1, 0])
                 collect(aux) .~ collect(cross(e_n_0, e_axis_0))
                 collect(e_long_0) .~ collect(aux ./ norm(aux))
                 collect(e_lat_0) .~ collect(cross(e_long_0, e_n_0))

                 # Determine point on road where the wheel is in contact with the road
                 collect(delta_0) .~ collect(r_road_0 - frame_a.r_0)
                 0 ~ delta_0'e_axis_0
                 0 ~ delta_0'e_long_0

                 # One holonomic positional constraint equation (no penetration in to the ground)
                 0 ~ radius - delta_0'cross(e_long_0, e_axis_0)

                 # only for testing
                 err ~ norm(delta_0) - radius

                 # Slip velocities
                 collect(v_0) .~ D.(frame_a.r_0)
                 collect(w_0) .~ angular_velocity1(ori(frame_a))
                 collect(vContact_0) .~ collect(v_0) + cross(w_0, delta_0)

                 # Two non-holonomic constraint equations on velocity level (ideal rolling, no slippage)
                 0 ~ vContact_0'e_long_0
                 0 ~ vContact_0'e_lat_0

                 # Contact force
                 f_wheel_0 ~ f_n * e_n_0 + f_lat * e_lat_0 + f_long * e_long_0

                 # Force and torque balance at the wheel center
                 zeros(3) .~ collect(frame_a.f) + resolve2(ori(frame_a), f_wheel_0)
                 zeros(3) .~ collect(frame_a.tau) +
                             resolve2(ori(frame_a), cross(delta_0, f_wheel_0))]
    compose(ODESystem(equations, t; name), frame_a)
end

"""
    RollingWheel(; name, radius, m, I_axis, I_long, width=0.035, x0, y0, kwargs...)

Ideal rolling wheel on flat surface z=0 (5 positional, 3 velocity degrees of freedom)

A wheel rolling on the x-y plane of the world frame including wheel mass.
The rolling contact is considered being ideal, i.e. there is no
slip between the wheel and the ground.
The wheel can not take off but it can incline toward the ground.
The frame frame_a is placed in the wheel center point and rotates
with the wheel itself.

# Arguments and parameters:
- `name`: Name of the rolling wheel component
- `radius`: Radius of the wheel
- `m`: Mass of the wheel
- `I_axis`: Moment of inertia of the wheel along its axis
- `I_long`: Moment of inertia of the wheel perpendicular to its axis
- `width`: Width of the wheel (default: 0.035)
- `x0`: Initial x-position of the wheel axis
- `y0`: Initial y-position of the wheel axis
- `kwargs...`: Additional keyword arguments passed to the `RollingWheelJoint` function

# Variables:
- `x`: x-position of the wheel axis
- `y`: y-position of the wheel axis
- `angles`: Angles to rotate world-frame into `frame_a` around z-, y-, x-axis
- `der_angles`: Derivatives of angles

# Named components:
- `frame_a`: Frame for the wheel component
- `rollingWheel`: Rolling wheel joint representing the wheel's contact with the road surface
"""
function RollingWheel(; name, radius, m, I_axis, I_long, width = 0.035, x0, y0,
                      angles = zeros(3), der_angles = zeros(3), kwargs...)
    @named begin
        frame_a = Frame()
        rollingWheel = RollingWheelJoint(; radius, angles, x0, y0, kwargs...)
        body = Body(r_cm = [0, 0, 0],
                    m = m,
                    I_11 = I_long,
                    I_22 = I_axis,
                    I_33 = I_long,
                    I_21 = 0,
                    I_31 = 0,
                    I_32 = 0)
    end
    pars = @parameters begin
        radius = radius, [description = "Radius of the wheel"]
        m = m, [description = "Mass of the wheel"]
        I_axis = I_axis, [description = "Moment of inertia of the wheel along its axis"]
        I_long = I_long,
                 [description = "Moment of inertia of the wheel perpendicular to its axis"]
        width = width, [description = "Width of the wheel"]
    end
    sts = @variables begin
        (x(t) = x0), [state_priority = 10, description = "x-position of the wheel axis"]
        (y(t) = y0), [state_priority = 10, description = "y-position of the wheel axis"]
        (angles(t)[1:3] = angles),
        [description = "Angles to rotate world-frame into frame_a around z-, y-, x-axis"]
        (der_angles(t)[1:3] = der_angles), [description = "Derivatives of angles"]
    end
    sts = reduce(vcat, collect.(sts))

    equations = Equation[rollingWheel.x ~ x
                         rollingWheel.y ~ y
                         collect(rollingWheel.angles) .~ collect(angles)
                         collect(rollingWheel.der_angles) .~ collect(der_angles)
                         connect(body.frame_a, frame_a)
                         connect(rollingWheel.frame_a, frame_a)]
    compose(ODESystem(equations, t, sts, pars; name), frame_a, rollingWheel, body)
end

"""
    FreeMotion(; name, enforceState = true, sequence, isroot = true, w_rel_a_fixed = false, z_rel_a_fixed = false, phi = 0, phi_d = 0, phi_dd = 0, w_rel_b = 0, r_rel_a = 0, v_rel_a = 0, a_rel_a = 0)

Joint which does not constrain the motion between `frame_a` and `frame_b`. Such a joint is only meaningful if the relative distance and orientation between `frame_a` and `frame_b`, and their derivatives, shall be used as state.

Note, that bodies such as [`Body`](@ref), [`BodyShape`](@ref), have potential state variables describing the distance and orientation, and their derivatives, between the world frame and a body fixed frame. Therefore, if these potential state variables are suited, a `FreeMotion` joint is not needed.

The state of the FreeMotion object consits of:

The relative position vector `r_rel_a` from the origin of `frame_a` to the origin of `frame_b`, resolved in `frame_a` and the relative velocity `v_rel_a` of the origin of `frame_b` with respect to the origin of `frame_a`, resolved in `frame_a (= der(r_rel_a))`.

# Arguments

- `enforceState`: Enforce this joint having state, this is often desired and is the default choice.
- `sequence`: Rotation sequence
- `w_rel_a_fixed`: = true, if `w_rel_a_start` are used as initial values, else as guess values
- `z_rel_a_fixed`: = true, if `z_rel_a_start` are used as initial values, else as guess values

# Initial codition arguments:
- `phi`
- `phi_d`
- `phi_dd`
- `w_rel_b`
- `r_rel_a`
- `v_rel_a`
- `a_rel_a`
"""
function FreeMotion(; name, enforceState = true, sequence = [1, 2, 3], isroot = true,
                    w_rel_a_fixed = false, z_rel_a_fixed = false, phi = 0,
                    phi_d = 0,
                    phi_dd = 0,
                    w_rel_b = 0,
                    r_rel_a = 0,
                    v_rel_a = 0,
                    a_rel_a = 0)
    @named begin
        frame_a = Frame()
        frame_b = Frame()
    end
    @variables begin
        (phi(t)[1:3] = phi),
        [state_priority = 10, description = "3 angles to rotate frame_a into frame_b"]
        (phi_d(t)[1:3] = phi_d), [state_priority = 10, description = "Derivatives of phi"]
        (phi_dd(t)[1:3] = phi_dd),
        [state_priority = 10, description = "Second derivatives of phi"]
        (w_rel_b(t)[1:3] = w_rel_b),
        [
            state_priority = 10,
            description = "relative angular velocity of frame_b with respect to frame_a, resolved in frame_b",
        ]
        (r_rel_a(t)[1:3] = r_rel_a),
        [
            description = "Position vector from origin of frame_a to origin of frame_b, resolved in frame_a",
        ]
        (v_rel_a(t)[1:3] = v_rel_a),
        [
            description = "= der(r_rel_a), i.e., velocity of origin of frame_b with respect to origin of frame_a, resolved in frame_a",
        ]
        (a_rel_a(t)[1:3] = a_rel_a), [description = "= der(v_rel_a)"]
    end

    @named R_rel = NumRotationMatrix()
    @named R_rel_inv = NumRotationMatrix()

    eqs = [
           # Cut-forces and cut-torques are zero
           frame_a.f .~ 0
           frame_a.tau .~ 0
           frame_b.f .~ 0
           frame_b.tau .~ 0
           D.(r_rel_a) .~ v_rel_a
           D.(v_rel_a) .~ a_rel_a

           # Kinematic relationships
           frame_b.r_0 .~ frame_a.r_0 .+ resolve1(frame_a, r_rel_a)]

    if enforceState
        if isroot
            append!(eqs,
                    ori(frame_b) ~ absoluteRotation(frame_a, R_rel))
        else
            append!(eqs,
                    [R_rel_inv ~ inverseRotation(R_rel)
                     ori(frame_a) ~ absoluteRotation(frame_b, R_rel_inv)])
        end

        append!(eqs,
                [phi_d .~ D.(phi)
                 phi_dd .~ D.(phi_d)
                 R_rel ~ axesRotations(sequence, phi, phi_d)
                 w_rel_b .~ angularVelocity2(R_rel)])

    else
        # Free motion joint does not have state
        if w_rel_a_fixed || z_rel_a_fixed
            append!(eqs,
                    w_rel_b .~ angularVelocity2(frame_b) - resolve2(frame_b.
                                        R, angularVelocity1(frame_a)))
        end
    end
    compose(ODESystem(eqs, t; name), frame_a, frame_b)
end
