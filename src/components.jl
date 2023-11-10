using LinearAlgebra

function isroot(sys)
    sys.metadata isa Dict || return false
    get(sys.metadata, :isroot, false)
end

"""
    ori(frame, varw = false)

Get the orientation of `sys` as a `RotationMatrix` object.

For frames, the orientation is stored in the metadata field of the system as `sys.metadata[:orientation]`.

If `varw = true`, the angular velocity variables `w` of the frame is also included in the `RotationMatrix` object, otherwise `w` is derived as the time derivative of `R`. `varw = true` is primarily used when selecting a component as root.
"""
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

@component function World(; name)
    # World should have
    # 3+3+9+3 // r_0+f+R.R+τ
    # - (3+3) // (f+t)
    # = 12 equations 
    @named frame_b = Frame()
    @parameters n[1:3]=[0, -1, 0] [description = "gravity direction of world"]
    @parameters g=9.81 [description = "gravitational acceleration of world"]
    O = ori(frame_b)
    eqs = Equation[collect(frame_b.r_0) .~ 0;
                   O ~ nullRotation()
                   # vec(D(O).R .~ 0); # QUESTION: not sure if I should have to add this, should only have 12 equations according to modelica paper
                   ]
    ODESystem(eqs, t, [], [n; g]; name, systems = [frame_b])
end

"""
The world component is the root of all multibody models. It is a fixed frame with a parallel gravitational field and a gravity vector specified by the unit direction `world.n` (defaults to [0, -1, 0]) and magnitude `world.g` (defaults to 9.81).
"""
const world = World(; name = :world)

"Compute the gravity acceleration, resolved in world frame"
gravity_acceleration(r) = GlobalScope(world.g) * GlobalScope.(world.n) # NOTE: This is hard coded for now to use the the standard, parallel gravity model

@component function Fixed(; name, r = [0, 0, 0])
    systems = @named begin frame_b = Frame() end
    @parameters begin r[1:3] = r,
                               [
                                   description = "Position vector from world frame to frame_b, resolved in world frame",
                               ] end
    eqs = [collect(frame_b.r_0 .~ r)
           ori(frame_b) ~ nullRotation()]
    compose(ODESystem(eqs, t; name), systems...)
end

@component function Mounting1D(; name, n = [1, 0, 0], phi0 = 0)
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

"""
    FixedTranslation(; name, r)

Fixed translation of `frame_b` with respect to `frame_a` with position vector `r` resolved in `frame_a`.

Can be though of as a massless rod. For a massive rod, see [`BodyShape`](@ref) or [`BodyCylinder`](@ref).
"""
@component function FixedTranslation(; name, r)
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
    FixedRotation(; name, r, n, sequence, isroot = false, angle, n_x, n_y)

Fixed translation followed by a fixed rotation of `frame_b` with respect to `frame_a`


- `r`: Translation vector
- `n`: Axis of rotation, resolved in frame_a
- `sequence`: DESCRIPTION
- `angle`: Angle of rotation around `n`, given in radians
"""
@component function FixedRotation(; name, r, n = [1, 0, 0], sequence = [1, 2, 3], isroot = false,
                       angle, n_x = [1, 0, 0], n_y = [0, 1, 0])
    norm(n) ≈ 1 || error("n must be a unit vector")
    @named frame_a = Frame()
    @named frame_b = Frame()
    @parameters r(t)[1:3]=r [
        description = "position vector from frame_a to frame_b, resolved in frame_a",
    ]
    @parameters n(t)[1:3]=n [
        description = "axis of rotation, resolved in frame_a",
    ]
    # @parameters n_x(t)=n_x [
    #     description = "Vector along x-axis of frame_b resolved in frame_a",
    # ]
    # @parameters n_y(t)=n_y [
    #     description = "Vector along y-axis of frame_b resolved in frame_a",
    # ]
    @parameters angle(t)=angle [
        description = "angle of rotation in radians",
    ]
    fa = frame_a.f |> collect
    fb = frame_b.f |> collect
    taua = frame_a.tau |> collect
    taub = frame_b.tau |> collect

    # Relationships between quantities of frame_a and frame_b 

    if isroot
        R_rel = planar_rotation(n, angle, 0)
        eqs = [ori(frame_b) ~ absoluteRotation(frame_a, R_rel);
               zeros(3) .~ fa + resolve1(R_rel, fb);
               zeros(3) .~ taua + resolve1(R_rel, taub) - cross(r,
                                                                fa)]
    else
        R_rel_inv = planar_rotation(n, -angle, 0)
        eqs = [ori(frame_a) ~ absoluteRotation(frame_b, R_rel_inv);
               zeros(3) .~ fb + resolve1(R_rel_inv, fa);
               zeros(3) .~ taub + resolve1(R_rel_inv, taua) +
                           cross(resolve1(R_rel_inv, r), fb)]
    end
    eqs = collect(eqs)
    append!(eqs, collect(frame_b.r_0) .~ collect(frame_a.r_0) + resolve1(frame_a, r))

    compose(ODESystem(eqs, t; name), frame_a, frame_b)
end

"""
    Body(; name, m = 1, r_cm, I = collect(0.001 * LinearAlgebra.I(3)), isroot = false, phi0 = zeros(3), phid0 = zeros(3), r_0=zeros(3))

Representing a body with 3 translational and 3 rotational degrees-of-freedom.

# Parameters
- `m`: Mass
- `r_cm`: Vector from `frame_a` to center of mass, resolved in `frame_a`
- `I`: Inertia matrix of the body
- `isroot`: Indicate whether this component is the root of the system, useful when there are no joints in the model.
- `phi0`: Initial orientation, only applicable if `isroot = true`
- `phid0`: Initial angular velocity

# Variables
- `r_0`: Position vector from origin of world frame to origin of `frame_a`
- `v_0`: Absolute velocity of `frame_a`, resolved in world frame (= D(r_0))
- `a_0`: Absolute acceleration of `frame_a` resolved in world frame (= D(v_0))
"""
@component function Body(; name, m = 1, r_cm = [0, 0, 0],
              I_11 = 0.001,
              I_22 = 0.001,
              I_33 = 0.001,
              I_21 = 0,
              I_31 = 0,
              I_32 = 0,
              isroot = false,
              phi0 = zeros(3),
              phid0 = zeros(3),
              r_0 = 0,
              radius = 0.005,
              useQuaternions=false,)
    @variables r_0(t)[1:3]=r_0 [
        state_priority = 2,
        description = "Position vector from origin of world frame to origin of frame_a",
    ]
    @variables v_0(t)[1:3]=0 [
        state_priority = 2,
        description = "Absolute velocity of frame_a, resolved in world frame (= D(r_0))",
    ]
    @variables a_0(t)[1:3]=0 [
        state_priority = 2,
        description = "Absolute acceleration of frame_a resolved in world frame (= D(v_0))",
    ]
    @variables g_0(t)[1:3]=0 [description = "gravity acceleration"]
    @variables w_a(t)[1:3]=0 [
        state_priority = 2,
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
    @parameters radius=radius [
        description = "Radius of the body in animations",
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

    dvs = [r_0;v_0;a_0;g_0;w_a;z_a;]
    eqs = if isroot # isRoot
        
        if useQuaternions
            @named frame_a = Frame(varw = true)
            Ra = ori(frame_a, true)
            # @variables q(t)[1:4] = [0.0,0,0,1.0]
            # @variables qw(t)[1:3] = [0.0,0,0]
            # q = collect(q)
            # qw = collect(qw)
            # Q = Quaternion(q, qw)
            @named Q = NumQuaternion(varw=true) 
            ar = from_Q(Q, angularVelocity2(Q, D.(Q.Q)))
            Equation[
                0 ~ orientation_constraint(Q)
                Ra ~ ar
                Ra.w .~ ar.w
                Q.w .~ ar.w
                collect(w_a .~ Ra.w)
            ]
        else
            @named frame_a = Frame(varw = true)
            @variables phi(t)[1:3]=phi0 [state_priority = 10, description = "Euler angles"]
            @variables phid(t)[1:3]=phid0 [state_priority = 10]
            @variables phidd(t)[1:3]=zeros(3) [state_priority = 10]
            phi, phid, phidd = collect.((phi, phid, phidd))
            ar = axesRotations([1, 2, 3], phi, phid)

            Ra = ori(frame_a, true)

            Equation[
                    # 0 .~ orientation_constraint(Ra); 
                    phid .~ D.(phi)
                    phidd .~ D.(phid)
                    Ra ~ ar
                    Ra.w .~ ar.w
                    collect(w_a .~ Ra.w)]
        end
    else
        @named frame_a = Frame()
        Ra = ori(frame_a)
        # This equation is defined here and not in the Rotation component since the branch above might use another equation
        Equation[
                 # collect(w_a .~ DRa.w); # angularVelocity2(R, D.(R)): skew(R.w) = R.T*der(transpose(R.T))
                 # vec(DRa.R .~ 0)
                 collect(w_a .~ angularVelocity2(Ra));]
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

    # pars = [m;r_cm;radius;I_11;I_22;I_33;I_21;I_31;I_32;]
    
    ODESystem(eqs, t; name, metadata = Dict(:isroot => isroot), systems = [frame_a])
end


"""
    BodyShape(; name, m = 1, r, kwargs...)

The `BodyShape` component is similar to a [`Body`](@ref), but it has two frames and a vector `r` that describes the translation between them, while the body has a single frame only.

- `r`: Vector from `frame_a` to `frame_b` resolved in `frame_a`
- All `kwargs` are passed to the internal `Body` component.
"""
@component function BodyShape(; name, m = 1, r = [0, 0, 0], r_cm = 0.5*r, r_0 = 0, radius = 0.08, kwargs...)
    systems = @named begin
        frameTranslation = FixedTranslation(r = r)
        body = Body(; r_cm, r_0, kwargs...)
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
    @parameters begin
        r[1:3]=r, [
            description = "Vector from frame_a to frame_b resolved in frame_a",
        ]
        radius = radius, [description = "Radius of the body in animations"]
        # color = color, [description = "Color of the body in animations"]
    end

    pars = [r; radius]

    r_0, v_0, a_0 = collect.((r_0, v_0, a_0))

    eqs = [r_0 .~ collect(frame_a.r_0)
           v_0 .~ D.(r_0)
           a_0 .~ D.(v_0)
           connect(frame_a, frameTranslation.frame_a)
           connect(frame_b, frameTranslation.frame_b)
           connect(frame_a, body.frame_a)]
    ODESystem(eqs, t, [r_0; v_0; a_0], pars; name, systems)
end
