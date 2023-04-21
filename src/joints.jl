
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
                       ori(frame_b) ~ abs_rotation(ori(frame_a), Rrel)
                       collect(frame_a.f) .~ -resolve1(Rrel, frame_b.f)
                       collect(frame_a.tau) .~ -resolve1(Rrel, frame_b.tau)]
    else
        eqs = Equation[Rrel ~ planar_rotation(-n, phi, w)
                       ori(frame_a) ~ abs_rotation(ori(frame_b), Rrel)
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
        @named fixed = Translational.Fixed()
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

function Universal(; name, n_a = [1, 0, 0], n_b = [0, 1, 0])
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
        (phi_a(t) = 0),
        [
            state_priority = 10,
            description = "Relative rotation angle from frame_a to intermediate frame",
        ]
        (phi_b(t) = 0),
        [
            state_priority = 10,
            description = "Relative rotation angle from intermediate frame to frame_b",
        ]
        (w_a(t) = 0),
        [
            state_priority = 10,
            description = "First derivative of angle phi_a (relative angular velocity a)",
        ]
        (w_b(t) = 0),
        [
            state_priority = 10,
            description = "First derivative of angle phi_b (relative angular velocity b)",
        ]
        (a_a(t) = 0),
        [
            state_priority = 10,
            description = "Second derivative of angle phi_a (relative angular acceleration a)",
        ]
        (a_b(t) = 0),
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
