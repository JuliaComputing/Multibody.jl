function add_params(sys, params; name)
    params isa AbstractVector || (params = [params...])
    extend(System(Equation[], t, [], params; name), sys)
end

"""
    Revolute(; name, phi0 = 0, w0 = 0, n, axisflange = false)

Revolute joint with 1 rotational degree-of-freedom

- `phi0`: Initial angle
- `w0`: Iniitial angular velocity
- `n`: The axis of rotation
- `axisflange`: If true, the joint will have two additional frames from Mechanical.Rotational, `axis` and `support`, between which rotational components such as springs and dampers can be connected.

If `axisflange`, flange connectors for ModelicaStandardLibrary.Mechanics.Rotational are also available:
- `axis`: 1-dim. rotational flange that drives the joint
- `support`: 1-dim. rotational flange of the drive support (assumed to be fixed in the world frame, NOT in the joint)

# Rendering options
- `radius = 0.05`: Radius of the joint in animations
- `length = radius`: Length of the joint in animations
- `color`: Color of the joint in animations, a vector of length 4 with values between [0, 1] providing RGBA values
"""
@component function Revolute(; name, phi0 = 0, w0 = 0, n = Float64[0, 0, 1], axisflange = false,
                  isroot = true, iscut = false, radius = 0.05, length = radius, color = [0.5019608f0,0.0f0,0.5019608f0,1.0f0], state_priority = 3.0, render = true)
    if !(eltype(n) <: Num) && !isa(n, Symbolics.Arr{Num, 1})
        norm(n) ≈ 1 || error("Axis of rotation must be a unit vector")
    end
    @named frame_a = Frame()
    @named frame_b = Frame()
    @parameters n[1:3]=n [description = "axis of rotation"]
    pars = @parameters begin
        radius = radius, [description = "radius of the joint in animations"]
        length = length, [description = "length of the joint in animations"]
        color[1:4] = color, [description = "color of the joint in animations (RGBA)"]
        render = render, [description = "render the joint in animations"]
    end
    @variables tau(t) [
        connect = Flow,
        # state_priority = 2,
        description = "Driving torque in direction of axis of rotation",
    ]
    @variables phi(t)=phi0 [
        state_priority = state_priority,
        description = "Relative rotation angle from frame_a to frame_b",
    ]
    @variables w(t)=w0 [state_priority = state_priority, description = "angular velocity (rad/s)"]
    Rrel0 = planar_rotation(n, phi0, w0)
    @named Rrel = NumRotationMatrix(; R = Rrel0.R, w = Rrel0.w)
    n = collect(n)


    if isroot
        eqs = Equation[Rrel ~ planar_rotation(n, phi, w)
                    connect_orientation(ori(frame_b), absolute_rotation(ori(frame_a), Rrel); iscut)
                    collect(frame_a.f) .~ -resolve1(Rrel, frame_b.f)
                    collect(frame_a.tau) .~ -resolve1(Rrel, frame_b.tau)]
    else
        eqs = Equation[Rrel ~ planar_rotation(-n, phi, w)
                    connect_orientation(ori(frame_a), absolute_rotation(ori(frame_b), Rrel); iscut)
                    collect(frame_b.f) .~ -resolve1(Rrel, frame_a.f)
                    collect(frame_b.tau) .~ -resolve1(Rrel, frame_a.tau)]
    end
    moreeqs = [collect(frame_a.r_0 .~ frame_b.r_0)
               D(phi) ~ w
               tau ~ -collect(frame_b.tau)'n]
    append!(eqs, moreeqs)
    sys = if axisflange
        # @named internalAxis = Rotational.InternalSupport(tau=tau)
        @named fixed = Rotational.Fixed()

        @named axis = Rotational.Flange()
        @named support = Rotational.Flange()
        # push!(eqs, phi ~ internalAxis.phi)
        push!(eqs, connect(fixed.flange, support))
        push!(eqs, axis.phi ~ phi)
        push!(eqs, axis.tau ~ tau)
        # push!(eqs, connect(internalAxis.flange, axis))
        System(eqs, t; name=:nothing, systems=[frame_a, frame_b, axis, support, fixed])
    else
        # Modelica Revolute uses a ConstantTorque as well as internalAxis = Rotational.InternalSupport(tau=tau), but it seemed more complicated than required and I couldn't get it to work, likely due to the `input` semantics of modelica not having an equivalent in MTK, so the (tau=tau) input argument caused problems.
        # @named constantTorque = Rotational.ConstantTorque(tau_constant=0, use_support=false) 
        # push!(eqs, connect(constantTorque.flange, internalAxis.flange))
        push!(eqs, tau ~ 0)
        System(eqs, t; name=:nothing, systems=[frame_a, frame_b])
    end
    add_params(sys, pars; name)
end

"""
    Prismatic(; name, n = [0, 0, 1], axisflange = false)

Prismatic joint with 1 translational degree-of-freedom

- `n`: The axis of motion (unit vector)
- `axisflange`: If true, the joint will have two additional frames from Mechanical.Translational, `axis` and `support`, between which translational components such as springs and dampers can be connected.

If `axisflange`, flange connectors for ModelicaStandardLibrary.Mechanics.TranslationalModelica are also available:
- `axis`: 1-dim. translational flange that drives the joint
- `support`: 1-dim. translational flange of the drive support (assumed to be fixed in the world frame, NOT in the joint)

The function returns an System representing the prismatic joint.
"""
@component function Prismatic(; name, n = Float64[0, 0, 1], axisflange = false,
                   s0 = 0, v0 = 0, radius = 0.05, color = [0,0.8,1,1], state_priority=10, iscut=false, render=true)
    if !(eltype(n) <: Num) && !isa(n, Symbolics.Arr{Num, 1})
        norm(n) ≈ 1 || error("Prismatic axis of motion must be a unit vector, got norm(n) = $(norm(n))")
    end
    @named frame_a = Frame()
    @named frame_b = Frame()
    @parameters n[1:3]=_normalize(n) [description = "axis of motion"]
    n = collect(n)

    pars = @parameters begin
        radius = radius, [description = "radius of the joint in animations"]
        color[1:4] = color, [description = "color of the joint in animations (RGBA)"]
        render = render, [description = "render the joint in animations"]
    end

    @variables s(t)=s0 [
        state_priority = state_priority,
        description = "Relative distance between frame_a and frame_b",
    ]
    @variables v(t)=v0 [
        state_priority = state_priority,
        description = "Relative velocity between frame_a and frame_b",
    ]
    @variables a(t)=0 [
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
           connect_orientation(ori(frame_b), ori(frame_a); iscut)

           # Force and torque balance
           zeros(3) .~ collect(frame_a.f + frame_b.f)
           zeros(3) .~ collect(frame_a.tau + frame_b.tau + cross(n * s, frame_b.f))

           # d'Alemberts principle
           f ~ -(n'collect(frame_b.f))[]]

    sys = if axisflange
        @named fixed = Translational.Fixed(s0=0)
        @named axis = Translational.Flange()
        @named support = Translational.Flange()
        push!(eqs, connect(fixed.flange, support))
        push!(eqs, axis.s ~ s)
        push!(eqs, axis.f ~ f)
        compose(System(eqs, t; name=:nothing), frame_a, frame_b, axis, support, fixed)
    else
        push!(eqs, f ~ 0)
        compose(System(eqs, t; name=:nothing), frame_a, frame_b)
    end
    add_params(sys, pars; name)
end

"""
    Spherical(; name, state = false, isroot = true, w_rel_a_fixed = false, z_rel_a_fixed = false, sequence, phi = 0, phid = 0, phidd = 0, d = 0)

Joint with 3 constraints that define that the origin of `frame_a` and the origin of `frame_b` coincide. By default this joint defines only the 3 constraints without any potential state variables. If parameter `state` is set to true, three states are introduced. The orientation of `frame_b` is computed by rotating `frame_a` along the axes defined in parameter vector `sequence` (default = [1,2,3], i.e., the Cardan angle sequence) around the angles used as state. If angles are used as state there is the slight disadvantage that a singular configuration is present leading to a division by zero.

- `isroot`: Indicate that `frame_a` is the root, otherwise `frame_b` is the root. Only relevant if `state = true`.
- `sequence`: Rotation sequence
- `d`: Viscous damping constant. If `d > 0`. the joint dissipates energy due to viscous damping according to ``τ ~ -d*ω``.

# Rendering options
- `radius = 0.1`: Radius of the joint in animations
- `color = [1,1,0,1]`: Color of the joint in animations, a vector of length 4 with values between [0, 1] providing RGBA values
- `render = true`: Render the joint in animations
"""
@component function Spherical(; name, state = false, isroot = true, iscut=false, w_rel_a_fixed = false,
                   z_rel_a_fixed = false, sequence = [1, 2, 3], phi = 0,
                   phid = 0,
                   d = 0,
                   neg_w = true,
                   phidd = 0,
                   color = [1, 1, 0, 1],
                   radius = 0.1,
                   quat = false,
                   render = true,
                   )

    dnum = d
    @named begin
        ptf = PartialTwoFrames()
        Rrel = NumRotationMatrix()
        Rrel_inv = NumRotationMatrix()
    end
    pars = @parameters begin
        radius = radius, [description = "radius of the joint in animations"]
        color[1:4] = color, [description = "color of the joint in animations (RGBA)"]
        d = d
        render = render, [description = "render the joint in animations"]
    end
    @unpack frame_a, frame_b = ptf
    # @parameters begin # Currently not using parameters due to these appearing in if statements
    #     sequence[1:3] = sequence
    # end
    @variables begin (w_rel(t)[1:3] = state ? zeros(3) : nothing),
                     [
                         description = "relative angular velocity of frame_b with respect to frame_a, resolved in frame_b",
                     ] end

    # torque balance
    if dnum <= 0
        eqs = [zeros(3) .~ collect(frame_a.tau)
            zeros(3) .~ collect(frame_b.tau)
            collect(frame_b.r_0) .~ collect(frame_a.r_0)]
    else
        fric = d*w_rel
        eqs = [-fric .~ collect(frame_a.tau)
        fric .~ resolve1(Rrel, collect(frame_b.tau))
        collect(frame_b.r_0) .~ collect(frame_a.r_0)]
    end

    if state
        if quat
            append!(eqs, nonunit_quaternion_equations(Rrel, w_rel; neg_w))
            # append!(eqs, collect(w_rel) .~ angularVelocity2(Rrel))
        else
            @variables begin
                (phi(t)[1:3] = phi),
                [state_priority = 10, description = "3 angles to rotate frame_a into frame_b"]
                (phid(t)[1:3] = phid),
                [state_priority = 10, description = "3 angle derivatives"]
                (phidd(t)[1:3] = phidd),
                [state_priority = 10, description = "3 angle second derivatives"]
            end
            append!(eqs,
                    [Rrel ~ axes_rotations(sequence, phi, phid)
                    collect(w_rel) .~ angular_velocity2(Rrel)
                    collect(phid .~ D.(phi))
                    collect(phidd .~ D.(phid))])
        end
        if isroot
            append!(eqs,
                    [connect_orientation(ori(frame_b), absolute_rotation(frame_a, Rrel); iscut)
                     zeros(3) .~ collect(frame_a.f) + resolve1(Rrel, frame_b.f)])
        else
            # NOTE: this branch should never happen
            append!(eqs,
                    [connect_orientation(Rrel_inv, inverse_rotation(Rrel); iscut)
                     ori(frame_a) ~ absolute_rotation(frame_b, Rrel_inv)
                     zeros(3) .~ collect(frame_b.f) + resolve2(Rrel, frame_a.f)])
        end

    else
        # Spherical joint does not have state
        append!(eqs,
                #frame_b.r_0 ~ transpose(frame_b.R.T)*(frame_b.R.T*(transpose(frame_a.R.T)*(frame_a.R.T*frame_a.r_0)));
                zeros(3) .~ collect(frame_a.f) +
                            resolve_relative(frame_b.f, frame_b, frame_a))
        if w_rel_a_fixed || z_rel_a_fixed
            append!(w_rel .~ angular_velocity2(frame_b) - resolve2(frame_b,
                                      angular_velocity1(frame_a)))
        else
            append!(w_rel .~ zeros(3))
        end
    end

    sys = extend(System(eqs, t; name=:nothing), ptf)
    add_params(sys, pars; name)
end


"""
    Universal(; name, n_a, n_b, phi_a = 0, phi_b = 0, w_a = 0, w_b = 0, a_a = nothing, a_b = nothing, state_priority=10)

Joint where `frame_a` rotates around axis `n_a` which is fixed in `frame_a` and `frame_b` rotates around axis `n_b` which is fixed in `frame_b`. The two frames coincide when `revolute_a.phi=0` and `revolute_b.phi=0`. This joint has the following potential states;

- The relative angle `phi_a = revolute_a.phi` [rad] around axis `n_a`
- the relative angle `phi_b = revolute_b.phi` [rad] around axis `n_b`
- the relative angular velocity `w_a = D(phi_a)`
- the relative angular velocity `w_b = D(phi_b)`
"""
@component function Universal(; name, n_a = [1, 0, 0], n_b = [0, 1, 0], phi_a = 0,
                   phi_b = 0,
                   state_priority = 10,
                   w_a = 0,
                   w_b = 0,
                   a_a = nothing,
                   a_b = nothing,
                   radius = 0.05f0,
                   length = radius, 
                   color = [1,0,0,1]
)
    @named begin
        ptf = PartialTwoFrames()
        revolute_a = Revolute(n = n_a, isroot = false, radius, length, color)
        revolute_b = Revolute(n = n_b, isroot = false, radius, length, color)
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
            state_priority = state_priority,
            description = "Relative rotation angle from frame_a to intermediate frame",
        ]
        (phi_b(t) = phi_b),
        [
            state_priority = state_priority,
            description = "Relative rotation angle from intermediate frame to frame_b",
        ]
        (w_a(t) = w_a),
        [
            state_priority = state_priority,
            description = "First derivative of angle phi_a (relative angular velocity a)",
        ]
        (w_b(t) = w_b),
        [
            state_priority = state_priority,
            description = "First derivative of angle phi_b (relative angular velocity b)",
        ]
        (a_a(t) = a_a),
        [
            state_priority = state_priority,
            description = "Second derivative of angle phi_a (relative angular acceleration a)",
        ]
        (a_b(t) = a_b),
        [
            state_priority = state_priority,
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
    extend(System(eqs, t; name, systems = [revolute_a, revolute_b]), ptf)
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
@component function GearConstraint(; name, ratio, checkTotalPower = false, n_a = [1, 0, 0],
                        n_b = [1, 0, 0], r_a = [0, 0, 0], r_b = [0, 0, 0])
    @named ptf = PartialTwoFrames()
    systems = @named begin
        bearing = Frame() #"Coordinate system fixed in the bearing"

        actuatedRevolute_a = Revolute(axisflange = true,
                                      n = n_a)
        actuatedRevolute_b = Revolute(axisflange = true,
                                      n = n_b)

        idealGear = Rotational.IdealGear(ratio = ratio)
        translation1 = FixedTranslation(r = r_b)
        translation2 = FixedTranslation(r = r_a)
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
        (a_b(t)),
        [
            state_priority = 10,
            description = "Relative angular acceleration of revolute joint at frame_b",
        ]

        (totalPower(t)), [description = "Total power flowing into this element"]
    end

    eqs = [phi_b ~ actuatedRevolute_b.phi
           w_b ~ D(phi_b)
           a_b ~ D(w_b)
           connect(actuatedRevolute_a.axis, idealGear.flange_a)
           connect(idealGear.flange_b, actuatedRevolute_b.axis)
           connect(actuatedRevolute_a.frame_a, translation2.frame_b)
           connect(translation2.frame_a, bearing)
           connect(translation1.frame_a, bearing)
           connect(translation1.frame_b, actuatedRevolute_b.frame_a)
           connect(frame_a, actuatedRevolute_a.frame_b)
           connect(actuatedRevolute_b.frame_b, frame_b)]

    # Measure power for test purposes
    if checkTotalPower
        push!(eqs,
              totalPower ~ frame_a.f'resolve2(frame_a, D.(frame_a.r_0)) +
                           frame_b.f'resolve2(frame_b, D.(frame_b.r_0)) +
                           bearing.f'resolve2(bearing, D.(bearing.r_0)) +
                           frame_a.tau'angular_velocity2(frame_a) +
                           frame_b.tau'angular_velocity2(frame_b) +
                           bearing.tau'angular_velocity2(bearing))
    end

    extend(System(eqs, t; name, systems), ptf)
end


"""
    FreeMotion(; name, state = true, sequence, isroot = true, w_rel_a_fixed = false, z_rel_a_fixed = false, phi = 0, phid = 0, phidd = 0, w_rel_b = 0, r_rel_a = 0, v_rel_a = 0, a_rel_a = 0)

Joint which _does not_ constrain the motion between `frame_a` and `frame_b`. Such a joint is only meaningful if the relative distance and orientation between `frame_a` and `frame_b`, and their derivatives, shall be used as state.

Note, that bodies such as [`Body`](@ref), [`BodyShape`](@ref), have potential state variables describing the distance and orientation, and their derivatives, between the world frame and a body fixed frame. Therefore, if these potential state variables are suited, a `FreeMotion` joint is not needed.

The state of the FreeMotion object consits of:

The relative position vector `r_rel_a` from the origin of `frame_a` to the origin of `frame_b`, resolved in `frame_a` and the relative velocity `v_rel_a` of the origin of `frame_b` with respect to the origin of `frame_a`, resolved in `frame_a (= D(r_rel_a))`.

# Arguments

- `state`: Enforce this joint having state, this is often desired and is the default choice.
- `sequence`: Rotation sequence, defaults to `[1, 2, 3]`
- `w_rel_a_fixed`: = true, if `w_rel_a_start` are used as initial values, else as guess values
- `z_rel_a_fixed`: = true, if `z_rel_a_start` are used as initial values, else as guess values

# Initial condition arguments:
- `phi`
- `phid`
- `phidd`
- `w_rel_b`
- `r_rel_a`
- `v_rel_a`
- `a_rel_a`
"""
@component function FreeMotion(; name, state = true, sequence = [1, 2, 3], isroot = true,
                    quat = false,
                    w_rel_a_fixed = false, z_rel_a_fixed = false, phi = 0,
                    iscut = false,
                    state_priority = 4,
                    phid = 0,
                    phidd = nothing,
                    neg_w = true,
                    w_rel_b = 0,
                    r_rel_a = 0,
                    v_rel_a = 0,
                    a_rel_a = nothing)
    @named begin
        frame_a = Frame()
        frame_b = Frame()
    end
    @variables begin
        (phi(t)[1:3] = phi),
        [state_priority = state_priority, description = "3 angles to rotate frame_a into frame_b"]
        (phid(t)[1:3] = phid), [state_priority = state_priority, description = "Derivatives of phi"]
        (phidd(t)[1:3] = phidd),
        [state_priority = state_priority, description = "Second derivatives of phi"]
        (w_rel_b(t)[1:3] = w_rel_b),
        [
            state_priority = quat ? state_priority : 1.0,
            description = "relative angular velocity of frame_b with respect to frame_a, resolved in frame_b",
        ]
        (r_rel_a(t)[1:3] = r_rel_a),
        [
            state_priority = state_priority,
            description = "Position vector from origin of frame_a to origin of frame_b, resolved in frame_a",
        ]
        (v_rel_a(t)[1:3] = v_rel_a),
        [
            state_priority = state_priority,
            description = "= D(r_rel_a), i.e., velocity of origin of frame_b with respect to origin of frame_a, resolved in frame_a",
        ]
        (a_rel_a(t)[1:3] = a_rel_a), [description = "= D(v_rel_a)"]
    end

    @named Rrel_f = Frame()
    @named Rrel_inv_f = Frame()
    Rrel = ori(Rrel_f)
    Rrel_inv = ori(Rrel_inv_f)

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

    if state
        if isroot
            append!(eqs,
                    connect_orientation(ori(frame_b), absolute_rotation(frame_a, Rrel); iscut))
        else
            append!(eqs,
                    [Rrel_inv ~ inverse_rotation(Rrel)
                     connect_orientation(ori(frame_a), absolute_rotation(frame_b, Rrel_inv); iscut)])
        end

        if quat
            append!(eqs, nonunit_quaternion_equations(Rrel, w_rel_b; neg_w))

        else
            append!(eqs,
                    [phid .~ D.(phi)
                    phidd .~ D.(phid)
                    Rrel ~ axes_rotations(sequence, phi, phid)
                    w_rel_b .~ (neg_w ? -1 : 1) * angular_velocity2(Rrel)])
        end

    else
        # Free motion joint does not have state
        if w_rel_a_fixed || z_rel_a_fixed
            append!(eqs,
                    w_rel_b .~ angular_velocity2(ori(frame_b)) - resolve2(ori(frame_b), angular_velocity1(ori(frame_a))))
        end
    end
    if state && !isroot
        compose(System(eqs, t; name), frame_a, frame_b, Rrel_f, Rrel_inv_f)
    elseif state
        compose(System(eqs, t; name), frame_a, frame_b, Rrel_f, )
    else
        compose(System(eqs, t; name), frame_a, frame_b)
    end
end

"""
    RevolutePlanarLoopConstraint(; name, n)

Revolute joint that is described by 2 positional constraints for usage in a planar loop (the ambiguous cut-force perpendicular to the loop and the ambiguous cut-torques are set arbitrarily to zero)

Joint where `frame_b` rotates around axis `n` which is fixed in `frame_a` and where this joint is used in a planar loop providing 2 constraint equations on position level.

If a planar loop is present, e.g., consisting of 4 revolute joints where the joint axes are all parallel to each other, then there is no unique mathematical solution if all revolute joints are modelled with [`Revolute`](@ref) and the symbolic algorithms will fail. The reason is that, e.g., the cut-forces in the revolute joints perpendicular to the planar loop are not uniquely defined when 3-dim. descriptions of revolute joints are used. In this case, one revolute joint in the loop has to be replaced by model `RevolutePlanarLoopConstraint`. The effect is that from the 5 constraints of a 3-dim. revolute joint, 3 constraints are removed and replaced by appropriate known variables (e.g., the force in the direction of the axis of rotation is treated as known with value equal to zero; for standard revolute joints, this force is an unknown quantity).
"""
@component function RevolutePlanarLoopConstraint(; name, n = Float64[0, 0, 1], radius = 0.05, length = radius, color = [0.5019608f0,0.0f0,0.5019608f0,1.0f0], render=true)
    norm(n) ≈ 1 || error("Axis of rotation must be a unit vector")
    @named frame_a = Frame()
    @named frame_b = Frame()

    # n isa Vector{Float64} || error("Parametric axis of rotation is currently not supported")

    # Activate this when symbolic parameters are a bit more robust
    # @parameters n[1:3]=n [description = "axis of rotation"]

    # # @parameters e[1:3] [description = "Unit vector in direction of rotation axis, resolved in frame_a (= same as in frame_b)"]
    # @parameters nnx_a[1:3] = [1,0,0]#ifelse(abs(n[1]) > 0.1, [0,1,0], ifelse(abs(n[2]) > 0.1, [0,0,1], [1,0,0])) [description = "Arbitrary vector that is not aligned with rotation axis n"]
    # @parameters ey_a[1:3] = [0,1,0]#inormalize(cross(n, nnx_a)) [description = "Unit vector orthogonal to axis n of revolute joint, resolved in frame_a"]
    # @parameters ex_a[1:3] = [1,0,0]#icross(ey_a, n) [description = "Unit vector orthogonal to axis n of revolute joint and to ey_a, resolved in frame_a"]
    # # @variables ey_b[1:3](t) [description = "ey_a, resolved in frame_b"]
    # # @variables ex_b[1:3](t) [description = "ex_a, resolved in frame_b"]


    pars = @parameters begin
        radius = radius, [description = "Radius of revolute cylinder in animations"]
        length = length, [description = "Length of revolute cylinder in animations"]
        color[1:4] = color, [description = "Color of revolute cylinder in animations"]
        render = render, [description = "Render the joint in animations"]
    end

    nnx_a = ifelse(abs(n[1]) > 0.1, [0,1,0], ifelse(abs(n[2]) > 0.1, [0,0,1], [1,0,0])) 
    ey_a = normalize(cross(n, nnx_a)) 
    ex_a = cross(ey_a, n) 
    
    
    @variables r_rel_a(t)[1:3] [description = "Position vector from origin of frame_a to origin of frame_b, resolved in frame_a"]
    @variables f_c(t)[1:2] [description = "Dummy or constraint forces in direction of ex_a, ey_a"]
    n0 = n
    @variables n(t)[1:3]
    

    # @named Rrel = NumRotationMatrix()

    Rrel0 = planar_rotation(n, 0, 0)
    varw = false
    @named Rrel = NumRotationMatrix(; R = Rrel0.R, w = Rrel0.w, varw, state_priority = -1)

    n = collect(n)
    ey_a = collect(ey_a)
    ex_a = collect(ex_a)
    r_rel_a = collect(r_rel_a)
    f_c = collect(f_c)

    Rb = ori(frame_b)

    eqs = [
        Rrel ~ relative_rotation(ori(frame_a), ori(frame_b))
        r_rel_a .~ resolve2(ori(frame_a), collect(frame_b.r_0 - frame_a.r_0))
        0 ~ (ex_a'r_rel_a)[]
        0 ~ (ey_a'r_rel_a)[]
        collect(frame_a.tau) .~ zeros(3)
        collect(frame_b.tau) .~ zeros(3)
        collect(frame_a.f) .~ vec([ex_a ey_a]*f_c)
        collect(frame_b.f) .~ -resolve2(Rrel, frame_a.f)
        collect(n) .~ n0
    ]
    sys = System(eqs, t; name=:nothing, systems=[frame_a, frame_b])
    add_params(sys, pars; name)
end

LinearAlgebra.normalize(a::Vector{Num}) = a / norm(a)


"""
    Planar(; n = [0,0,1], n_x = [1,0,0], cylinderlength = 0.1, cylinderdiameter = 0.05, cylindercolor = [1, 0, 1, 1], boxwidth = 0.3*cylinderdiameter, boxheight = boxwidth, boxcolor = [0, 0, 1, 1])

Joint where `frame_b` can move in a plane and can rotate around an
axis orthogonal to the plane. The plane is defined by
vector `n` which is perpendicular to the plane and by vector `n_x`,
which points in the direction of the x-axis of the plane.
`frame_a` and `frame_b` coincide when `s_x=prismatic_x.s=0,
s_y=prismatic_y.s=0` and `phi=revolute.phi=0`.

# Structural parameters
- `n`: Axis orthogonal to unconstrained plane, resolved in `frame_a` (= same as in `frame_b`)
- `n_x`: Vector in direction of x-axis of plane, resolved in `frame_a` (`n_x` shall be orthogonal to `n`)

# Connectors
- `frame_a`: Frame for the joint
- `frame_b`: Frame for the joint

# Variables
- `s_x`: Relative distance along first prismatic joint starting at `frame_a`
- `s_y`: Relative distance along second prismatic joint starting at first prismatic joint
- `phi`: Relative rotation angle from `frame_a` to `frame_b`
- `v_x`: Relative velocity along first prismatic joint
- `v_y`: Relative velocity along second prismatic joint
- `w`: Relative angular velocity around revolute joint
- `a_x`: Relative acceleration along first prismatic joint
- `a_y`: Relative acceleration along second prismatic joint
- `wd`: Relative angular acceleration around revolute joint

# Rendering parameters
- `cylinderlength`: Length of the revolute cylinder
- `cylinderdiameter`: Diameter of the revolute cylinder
- `cylindercolor`: (structural) Color of the revolute cylinder
- `boxwidth`: Width of the prismatic joint boxes
- `boxheight`: Height of the prismatic joint boxes
- `boxcolor`: (structural) Color of the prismatic joint boxes
- `radius`: (structural) Radius of the revolute cylinder
- `render`: Enable rendering of the joint in animations
"""
@component function Planar(; name, state_priority = 1, n, n_x, radius = 0.05,
                           cylindercolor = [1, 0, 1, 1], boxcolor = [0, 0, 1, 1],
                           cylinderlength = 0.1, cylinderdiameter = 0.05,
                           boxwidth = 0.3*cylinderdiameter, boxheight = boxwidth, render = true)
    pars = @parameters begin
        cylinderlength = cylinderlength, [description = "Length of revolute cylinder"]
        cylinderdiameter = cylinderdiameter, [description = "Diameter of revolute cylinder"]
        boxwidth = boxwidth, [description = "Width of prismatic joint boxes"]
        boxheight = boxheight, [description = "Height of prismatic joint boxes"]
        render = render, [description = "Enable rendering of the joint in animations"]
    end

    # n = collect(n)
    # n_x = collect(n_x)

    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
        prismatic_x = Prismatic(; state_priority=2.1, n=cross(cross(n, n_x), n), color=boxcolor, radius)
        prismatic_y = Prismatic(; state_priority=2.1, n=cross(n, n_x), color=boxcolor, radius)
        revolute = Revolute(; state_priority=2.1, n, isroot=false, color=cylindercolor, radius)
    end

    vars = @variables begin
        (s_x(t) = 0), [state_priority = 3.0, description = "Relative distance along first prismatic joint starting at frame_a"]
        (s_y(t) = 0), [state_priority = 3.0, description = "Relative distance along second prismatic joint starting at first prismatic joint"]
        (phi(t) = 0), [state_priority = 3.0, description = "Relative rotation angle from frame_a to frame_b"]
        (v_x(t) = 0), [state_priority = 3.0, description = "Relative velocity along first prismatic joint"]
        (v_y(t) = 0), [state_priority = 3.0, description = "Relative velocity along second prismatic joint"]
        (w(t) = 0), [state_priority = 3.0, description = "Relative angular velocity around revolute joint"]
        (a_x(t)), [description = "Relative acceleration along first prismatic joint"]
        (a_y(t)), [description = "Relative acceleration along second prismatic joint"]
        (wd(t)), [description = "Relative angular acceleration around revolute joint"]
    end

    equations = [
        s_x ~ prismatic_x.s
        s_y ~ prismatic_y.s
        phi ~ revolute.phi
        v_x ~ D(s_x)
        v_y ~ D(s_y)
        w   ~ D(phi)
        a_x ~ D(v_x)
        a_y ~ D(v_y)
        wd  ~ D(w)
        connect(frame_a, prismatic_x.frame_a)
        connect(prismatic_x.frame_b, prismatic_y.frame_a)
        connect(prismatic_y.frame_b, revolute.frame_a)
        connect(revolute.frame_b, frame_b)
    ]

    return System(equations, t; name, systems)
end

@component function Cylindrical(; name, n = [1, 0, 0], cylinder_color = [1, 0, 1, 1],
                                cylinder_diameter = 0.05, render = true)
    pars = @parameters begin
        n[1:3] = n, [description = "Cylinder axis resolved in frame_a (= same as in frame_b)"]
        cylinder_diameter = cylinder_diameter, [description = "Diameter of cylinder"]
        render = render, [description = "Enable rendering of the joint in animations"]
    end

    # n = collect(n)
    # cylinder_color = collect(cylinder_color)

    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
        prismatic = Prismatic(; n, state_priority=1, render = false)
        revolute = Revolute(; n, state_priority=1, color = cylinder_color, radius = cylinder_diameter/2)
    end

    vars = @variables begin
        (s(t) = 0), [state_priority = 200, description = "Relative distance between frame_a and frame_b"]
        (phi(t) = 0), [state_priority = 200, description = "Relative rotation angle from frame_a to frame_b"]
        (v(t) = 0), [state_priority = 200, description = "First derivative of s (relative velocity)"]
        (w(t) = 0), [state_priority = 200, description = "First derivative of angle phi (relative angular velocity)"]
        (a(t)), [description = "Second derivative of s (relative acceleration)"]
        (wd(t)), [description = "Second derivative of angle phi (relative angular acceleration)"]
    end

    equations = [
        phi ~ revolute.phi
        w ~ D(phi)
        wd ~ D(w)
        s ~ prismatic.s
        v ~ D(s)
        a ~ D(v)
        connect(frame_a, prismatic.frame_a)
        connect(prismatic.frame_b, revolute.frame_a)
        connect(revolute.frame_b, frame_b)
    ]

    return System(equations, t; name, systems)
end

@component function URDFRevolute(; name, r, R=I(3), axisflange = false, kwargs...)
    if R == I(3) && r == zeros(3)
        return Revolute(; axisflange, name, kwargs...)
    end

    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
        rev = Revolute(; axisflange, kwargs...)
    end
    if R == I(3)
        @named trans = FixedTranslation(; r, render=false)
    else
        R = RotMatrix{3}(R)
        n = rotation_axis(R)
        angle = rotation_angle(R)
        @named trans = FixedRotation(; r, n, angle, render=false)
    end
    push!(systems, trans)
    connections = [
        connect(frame_a, trans.frame_a)
        connect(trans.frame_b, rev.frame_a)
        connect(rev.frame_b, frame_b)
    ]   
    if axisflange
        more_systems = @named begin
            axis = Rotational.Flange()
            support = Rotational.Flange()
        end
        systems = [systems; more_systems]
        connections = [
            connections
            connect(axis, rev.axis)
            connect(support, rev.support)
        ]
    end
    System(connections, t; systems, name)
end

@component function URDFPrismatic(; name, r, R, axisflange = false, kwargs...)
    if R == I(3) && r == zeros(3)
        return Prismatic(; axisflange, name, kwargs...)
    end

    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
        rev = Prismatic(; axisflange, kwargs...)
    end
    if R == I(3)
        @named trans = FixedTranslation(; r, render=false)
    else
        R = RotMatrix{3}(R)
        n = rotation_axis(R)
        angle = rotation_angle(R)
        @named trans = FixedRotation(; r, n, angle, render=false)
    end
    push!(systems, trans)
    connections = [
        connect(frame_a, trans.frame_a)
        connect(trans.frame_b, rev.frame_a)
        connect(rev.frame_b, frame_b)
    ]   
    if axisflange
        more_systems = @named begin
            axis = Rotational.Flange()
            support = Rotational.Flange()
        end
        systems = [systems; more_systems]
        connections = [
            connections
            connect(axis, rev.axis)
            connect(support, rev.support)
        ]
    end
    System(connections, t; systems, name)
end

@component function NullJoint(; name)
    pars = @parameters begin
    end

    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
    end

    vars = @variables begin
    end

    equations = [
        connect(frame_a, frame_b)
    ]

    return System(equations, t; name, systems)
end