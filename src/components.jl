using LinearAlgebra

function isroot(sys)
    sys.metadata isa Dict || return false
    get(sys.metadata, :isroot, false)
end

function ori(sys, varw = false)
    if sys.metadata isa Dict && (O = get(sys.metadata, :orientation, nothing)) !== nothing
        R = collect(O.R)
        # Since we are using this function instead of sys.ori, we need to handle namespacing properly as well
        ns = nameof(sys)
        R = ModelingToolkit.renamespace.(ns, R) .|> Num
        if varw
            w = collect(O.w)
            w = ModelingToolkit.renamespace.(ns, w) .|> Num
        else
            w = get_w(R)
        end
        RotationMatrix(R, w)
    else
        error("System $(sys.name) does not have an orientation object.")
    end
end

function World(; name)
    # World should have
    # 3+3+9+3 // r_0+f+R.R+τ
    # - (3+3) // (f+t)
    # = 12 equations 
    @named frame_b = Frame()
    @parameters n[1:3]=[0, -1, 0] [description = "gravity direction of world"]
    @parameters g=9.81 [description = "gravitational acceleration of world"]
    O = ori(frame_b)
    eqs = Equation[collect(frame_b.r_0) .~ 0;
                   O ~ nullrotation()
                   # vec(D(O).R .~ 0); # QUESTION: not sure if I should have to add this, should only have 12 equations according to modelica paper
                   ]
    ODESystem(eqs, t, [], [n; g]; name, systems = [frame_b])
end

"""
The world component is the root of all multibody models. It is a fixed frame with a parallel gravitational field and a gravity vector specified by the unit direction `world.n` (defaults to [0, -1, 0]) and magnitude `world.g` (defaults to 9.81).
"""
const world = World(; name = :world)

"Compute the gravity acceleration, resolved in world frame"
gravity_acceleration(r) = world.g * world.n # NOTE: This is hard coded for now to use the the standard, parallel gravity model

function Fixed(; name, r = [0, 0, 0])
    systems = @named begin frame_b = Frame() end
    @parameters begin r[1:3] = r,
                               [
                                   description = "Position vector from world frame to frame_b, resolved in world frame",
                               ] end
    eqs = [collect(frame_b.r_0 .~ r)
           ori(frame_b) ~ nullRotation()]
    compose(ODESystem(eqs, t; name), systems...)
end

function Mounting1D(; name, n = [1, 0, 0], phi0 = 0)
    systems = @named begin
        flange_b = Rotational.Flange()
        frame_a = Frame()
        housing_frame_a = Frame()
    end
    @parameters begin
        phi0 = phi0, [
                   description = "Fixed offset angle of housing"]
        n[1:3] = n,
                 [
                     description = "Axis of rotation = axis of support torque (resolved in frame_a)",
                 ]
    end
    @variables begin (housing_tau(t)[1:3] = 0), [
                         description = "Torque",
                     ] end
    eqs = [collect(housing_tau) .~ collect(-n * flange_b.tau)
           flange_b.phi .~ phi0
           connect(housing_frame_a, frame_a)]
    compose(ODESystem(eqs, t; name), systems...)
end

function FixedTranslation(; name, r)
    @named frame_a = Frame()
    @named frame_b = Frame()
    @parameters r(t)[1:3]=r [
        description = "position vector from frame_a to frame_b, resolved in frame_a",
    ]
    fa = frame_a.f |> collect
    fb = frame_b.f |> collect
    taua = frame_a.tau |> collect
    taub = frame_b.tau |> collect
    eqs = Equation[collect(frame_b.r_0) .~ collect(frame_a.r_0) + resolve1(ori(frame_a), r)
                   ori(frame_b) ~ ori(frame_a)
                   0 .~ fa + fb
                   0 .~ taua + taub + cross(r, fb)]
    compose(ODESystem(eqs, t; name), frame_a, frame_b)
end

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
    Body(; name, m = 1, r_cm, I = collect(0.001 * LinearAlgebra.I(3)), isroot = false, phi0 = zeros(3), phid0 = zeros(3), r_0=zeros(3))

Representing a body with 3 translational and 3 rotational degrees-of-freedom.

- `m`: Mass
- `r_cm`: Vector from `frame_a` to center of mass, resolved in `frame_a`
- `I`: Inertia matrix of the body
- `isroot`: Indicate whether this component is the root of the system, useful when there are no joints in the model.
- `phi0`: Initial orientation, only applicable if `isroot = true`
- `phid0`: Initial angular velocity
"""
function Body(; name, m = 1, r_cm = [0, 0, 0],
              I_11 = 0.001,
              I_22 = 0.001,
              I_33 = 0.001,
              I_21 = 0,
              I_31 = 0,
              I_32 = 0,
              isroot = false, phi0 = zeros(3), phid0 = zeros(3), r_0 = 0)
    @variables r_0(t)[1:3]=r_0 [
        state_priority = 2,
        description = "Position vector from origin of world frame to origin of frame_a",
    ]
    @variables v_0(t)[1:3]=0 [
        state_priority = 2,
        description = "Absolute velocity of frame_a, resolved in world frame (= D(r_0))",
    ]
    @variables a_0(t)[1:3]=0 [
        description = "Absolute acceleration of frame_a resolved in world frame (= D(v_0))",
    ]
    @variables g_0(t)[1:3]=0 [description = "gravity acceleration"]
    @variables w_a(t)[1:3]=0 [
        description = "Absolute angular velocity of frame_a resolved in frame_a",
    ]
    @variables z_a(t)[1:3]=0 [
        description = "Absolute angular acceleration of frame_a resolved in frame_a",
    ]
    # 6*3 potential variables + Frame: 2*3 flow + 3 potential + 3 residual = 24 equations + 2*3 flow
    @parameters m=m [description = "mass"]
    @parameters r_cm[1:3]=r_cm [
        description = "Vector from frame_a to center of mass, resolved in frame_a",
    ]
    # @parameters I[1:3, 1:3]=I [description="inertia tensor"]

    @parameters I_11=I_11 [description = "Element (1,1) of inertia tensor"]
    @parameters I_22=I_22 [description = "Element (2,2) of inertia tensor"]
    @parameters I_33=I_33 [description = "Element (3,3) of inertia tensor"]
    @parameters I_21=I_21 [description = "Element (2,1) of inertia tensor"]
    @parameters I_31=I_31 [description = "Element (3,1) of inertia tensor"]
    @parameters I_32=I_32 [description = "Element (3,2) of inertia tensor"]

    I = [I_11 I_21 I_31; I_21 I_22 I_32; I_31 I_32 I_33]

    r_0, v_0, a_0, g_0, w_a, z_a, r_cm = collect.((r_0, v_0, a_0, g_0, w_a, z_a, r_cm))

    # DRa = D(Ra)

    eqs = if isroot # isRoot
        @variables phi(t)[1:3]=phi0 [state_priority = 10, description = "Euler angles"]
        @variables phid(t)[1:3]=phid0 [state_priority = 10]
        @variables phidd(t)[1:3]=zeros(3) [state_priority = 10]
        phi, phid, phidd = collect.((phi, phid, phidd))
        ar = axesRotations(phi, phid)

        @named frame_a = Frame(varw = true)
        Ra = ori(frame_a, true)

        Equation[
                 # 0 .~ orientation_constraint(Ra); 
                 phid .~ D.(phi)
                 phidd .~ D.(phid)
                 Ra ~ ar
                 Ra.w .~ ar.w
                 collect(w_a .~ Ra.w)]
    else
        @named frame_a = Frame()
        Ra = ori(frame_a)
        # This equation is defined here and not in the Rotation component since the branch above might use another equation
        Equation[
                 # collect(w_a .~ DRa.w); # angular_velocity2(R, D.(R)): skew(R.w) = R.T*der(transpose(R.T))
                 # vec(DRa.R .~ 0)
                 collect(w_a .~ angular_velocity2(Ra));]
    end

    eqs = [eqs;
           # collect(w_a .~ get_w(Ra));
           collect(r_0 .~ frame_a.r_0)
           collect(g_0 .~ gravity_acceleration(frame_a.r_0 .+ resolve1(Ra, r_cm)))
           collect(v_0 .~ D.(r_0))
           collect(a_0 .~ D.(v_0))
           collect(z_a .~ D.(w_a))
           collect(frame_a.f .~ m * (resolve2(Ra, a_0 - g_0) + cross(z_a, r_cm) +
                                 cross(w_a, cross(w_a, r_cm))))
           collect(frame_a.tau .~ I * z_a + cross(w_a, I * w_a) + cross(r_cm, frame_a.f))]

    ODESystem(eqs, t; name, metadata = Dict(:isroot => isroot), systems = [frame_a])
end

"""
    BodyShape(; name, m = 1, r, kwargs...)

The `BodyShape` component is similar to a [`Body`](@ref), but it has two frames and a vector `r` that describes the translation between them, while the body has a single frame only.

- `r`: Vector from `frame_a` to `frame_b` resolved in `frame_a`
- All `kwargs` are passed to the internal `Body` component.
"""
function BodyShape(; name, m = 1, r = [0, 0, 0], r_0 = 0, kwargs...)
    systems = @named begin
        frameTranslation = FixedTranslation(r = r)
        body = Body(; kwargs...)
        frame_a = Frame()
        frame_b = Frame()
    end

    @variables r_0(t)[1:3]=r_0 [
        state_priority = 2,
        description = "Position vector from origin of world frame to origin of frame_a",
    ]
    @variables v_0(t)[1:3]=0 [
        state_priority = 2,
        description = "Absolute velocity of frame_a, resolved in world frame (= D(r_0))",
    ]
    @variables a_0(t)[1:3]=0 [
        description = "Absolute acceleration of frame_a resolved in world frame (= D(v_0))",
    ]
    @parameters r[1:3]=r [
        description = "Vector from frame_a to frame_b resolved in frame_a",
    ]

    r_0, v_0, a_0 = collect.((r_0, v_0, a_0))

    eqs = [r_0 .~ collect(frame_a.r_0)
           v_0 .~ D.(r_0)
           a_0 .~ D.(v_0)
           connect(frame_a, frameTranslation.frame_a)
           connect(frame_b, frameTranslation.frame_b)
           connect(frame_a, body.frame_a)]
    ODESystem(eqs, t; name, systems)
end
