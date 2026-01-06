using LinearAlgebra
using ModelingToolkit: get_metadata
import ModelingToolkitStandardLibrary

struct IsRoot end

function isroot(sys)
    md = get_metadata(sys)
    md isa Dict || return false
    get(md, IsRoot, false)
end

purple = [0.5019608f0,0.0f0,0.5019608f0,1.0f0]

"""
    ori(frame, varw = false)

Get the orientation of `sys` as a `RotationMatrix` object. See also [`get_rot`](@ref). `ori(frame).R` is the rotation matrix that rotates a vector from the world coordinate system to the local frame.

For frames, the orientation is stored in the metadata field of the system as `get_metadata(sys)[:orientation]`.

If `varw = true`, the angular velocity variables `w` of the frame is also included in the `RotationMatrix` object, otherwise `w` is derived as the time derivative of `R`. `varw = true` is primarily used when selecting a component as root.
"""
function ori(sys, varw = false)
    md = get_metadata(sys)
    if md isa AbstractDict && (O = get(md, ModelingToolkit.FrameOrientation, nothing)) !== nothing
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
        error("System $(getfield(sys, :name)) does not have an orientation object.")
    end
end

"""
    World(; name, render=true, point_gravity=false, n = [0.0, -1.0, 0.0], g=9.80665, mu=3.986004418e14)

All multibody models must include exactly one world component defined at the top level. The `frame_b` of the world is fixed in the origin.

If a connection to the world is needed in a component model, use [`Fixed`](@ref) instead.

# Arguments
- `name`: Name of the world component
- `render`: Render the component in animations
- `point_gravity`: If `true`, the gravity is always opinting towards a single point in space. If `false`, the gravity is always pointing in the same direction `n`.
- `n`: Gravity direction unit vector, defaults to [0, -1, 0], only applicable if `point_gravity = false`
- `g`: Gravitational acceleration, defaults to 9.80665
- `mu`: Gravity field constant, defaults to 3.986004418e14, only applicable to point gravity
"""
@component function World(; name, render=true, point_gravity=false, n = [0.0, -1.0, 0.0], g=9.80665, mu=3.986004418e14)
    # World should have
    # 3+3+9+3 // r_0+f+R.R+τ
    # - (3+3) // (f+t)
    # = 12 equations 
    n0 = n
    g0 = g
    mu0 = mu
    @named frame_b = Frame()

    @parameters n[1:3] = n0 [description = "gravity direction"]
    @parameters g=g0 [description = "gravitational acceleration of world"]
    @parameters mu=mu0 [description = "Gravity field constant [m³/s²] (default = field constant of earth)"]
    @parameters render=render
    @parameters point_gravity = point_gravity

    @variables n_inner(t)[1:3]
    @variables g_inner(t)
    @variables mu_inner(t)
    @variables render_inner(t)
    @variables point_gravity_inner(t)

    n = Symbolics.scalarize(n)
    n_inner = GlobalScope.(Symbolics.scalarize(n_inner))
    g_inner = GlobalScope(g_inner)
    mu_inner = GlobalScope(mu_inner)
    render_inner = GlobalScope(render_inner)
    point_gravity_inner = GlobalScope(point_gravity_inner)

    O = ori(frame_b)
    eqs = Equation[
        frame_b.r_0 ~ zeros(3)
        O ~ nullrotation()
        n_inner .~ n
        g_inner ~ g
        mu_inner ~ mu
        render_inner ~ render
        point_gravity_inner ~ point_gravity
    ]
    System(eqs, t, [n_inner; g_inner; mu_inner; render_inner; point_gravity_inner], [n; g; mu; point_gravity; render]; name, systems = [frame_b])#, defaults=[n => n0; g => g0; mu => mu0])
end

"""
The world component is the root of all multibody models. It is a fixed frame with a parallel gravitational field and a gravity vector specified by the unit direction `world.n` (defaults to [0, -1, 0]) and magnitude `world.g` (defaults to 9.80665).
"""
const world = World(; name = :world)

"Compute the gravity acceleration, resolved in world frame"
function gravity_acceleration(r)
    inner_gravity(GlobalScope(world.point_gravity_inner), GlobalScope(world.mu_inner), GlobalScope(world.g_inner), GlobalScope.(collect(world.n_inner)), collect(r))
end

function inner_gravity(point_gravity, mu, g, n, r)
    # This is slightly inefficient, producing three if statements, one for each array entry. The function registration for array-valued does not work properly so this is a workaround for now. Hitting, among other problems, https://github.com/SciML/ModelingToolkit.jl/issues/2808
    gvp = -(mu/(r'r))*(r/_norm(r))
    gvu = g * n
    ifelse.(point_gravity==true, gvp, gvu)
end


"""
    Fixed(; name, r = [0, 0, 0], render = true)

A component rigidly attached to the world frame with translation `r` between the world and the `frame_b`. The position vector `r` is resolved in the world frame.
"""
@component function Fixed(; name, r = [0, 0, 0], render = true)
    systems = @named begin frame_b = Frame() end
    @parameters begin
        r[1:3] = r, [description = "Position vector from world frame to frame_b, resolved in world frame"]
        render = render, [description = "Render the component in animations"]
    end
    eqs = [frame_b.r_0 ~ r
           ori(frame_b) ~ nullrotation()]
    sys = compose(System(eqs, t; name=:nothing), systems...)
    add_params(sys, [render]; name)
end

"""
    Position(; name, r = [0, 0, 0], render = true, fixed_orientation = true)

Forced movement of a flange according to a reference position `r`. Similar to [`Fixed`](@ref), but `r` is not a parameter, and may thus be any time-varying symbolic expression. The reference position vector `r` is resolved in the world frame. Example: `r = [sin(t), 0, 0]`.

# Arguments:
- `r`: Position vector from world frame to frame_b, resolved in world frame
- `render`: Render the component in animations
- `fixed_orientation`: If `true`, the orientation of the frame is fixed to the world frame. If `false`, the orientation is free to change.

See also [`Pose`](@ref) for a component that allows for forced orientation as well.
"""
@component function Position(; name, r = [0, 0, 0], render = true, fixed_orientation = true, x_locked = true, y_locked = true, z_locked = true)
    systems = @named begin frame_b = Frame() end
    @parameters begin
        render = render, [description = "Render the component in animations"]
    end
    @variables begin
        p(t)[1:3], [description = "Position vector from world frame to frame_b, resolved in world frame"]
        v(t)[1:3], [description = "Absolute velocity of frame_b, resolved in world frame"]
        a(t)[1:3], [description = "Absolute acceleration of frame_b, resolved in world frame"]
    end
    eqs = [
        (frame_b.r_0 .~ r)[[x_locked, y_locked, z_locked]]
        frame_b.r_0 ~ p
        v ~ D(p)
        a ~ D(v)
        ]
    if fixed_orientation
        append!(eqs, ori(frame_b) ~ nullrotation())
    end
    sys = compose(System(eqs, t; name=:nothing), systems...)
    add_params(sys, [render]; name)
end

"""
    Pose(; name, r = [0, 0, 0], R, q, render = true)

Forced movement of a flange according to a reference position `r` and reference orientation `R`. The reference arrays `r` and `R` are resolved in the world frame, and may be any symbolic expression. As an alternative to specifying `R`, it is possible to specify a quaternion `q` (4-vector quaternion with real part first).

Example usage:
```
using Multibody.Rotations
R = RotXYZ(0, 0.5sin(t), 0)
@named rot = Multibody.Pose(; r=[sin(t), 0, 0], R)
```

# Connectors
- `frame_b`: The frame that is forced to move according to the reference position and orientation.

# Arguments 
- `r`: Position vector from world frame to frame_b, resolved in world frame
- `R`: Reference orientation matrix
- `q`: Reference quaternion (optional alternative to `R`)
- `render`: Render the component in animations
- `normalize`: If a quaternion is provided, normalize the quaternion (default = true)
- `x_locked`, `y_locked`, `z_locked`: Lock the translation in the x, y, and z directions, respectively. This allows for selective specification of the translation components, i.e., if `y_locked = false`, the y-component of the translation is not constrained to follow `r`.

See also [`Position`](@ref) for a component that allows for only forced translation.
"""
@component function Pose(; name, r = [0, 0, 0], R=nothing, q=nothing, render = true, normalize=true, x_locked = true, y_locked = true, z_locked = true)
    systems = @named begin frame_b = Frame() end
    @parameters begin
        render = render, [description = "Render the component in animations"]
    end
    @variables begin
        p(t)[1:3], [description = "Position vector from world frame to frame_b, resolved in world frame"]
        v(t)[1:3], [description = "Absolute velocity of frame_b, resolved in world frame"]
        a(t)[1:3], [description = "Absolute acceleration of frame_b, resolved in world frame"]
    end
    eqs = [
        (frame_b.r_0 .~ r)[[x_locked, y_locked, z_locked]]
        frame_b.r_0 ~ p
        v ~ D(p)
        a ~ D(v)
        if R !== nothing
            ori(frame_b).R ~ R
        elseif q !== nothing
            ori(frame_b) ~ from_Q(q, 0; normalize)
        else
            error("Either R or q must be provided. If you only want to specify the position, use the `Position` component instead.")
        end
    ]
    sys = compose(System(eqs, t; name=:nothing), systems...)
    add_params(sys, [render]; name)
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
    @variables begin (housing_tau(t)[1:3]), [
                         description = "Torque",
                     ] end
    eqs = [housing_tau ~ -n * flange_b.tau
           flange_b.phi ~ phi0
           connect(housing_frame_a, frame_a)]
    compose(System(eqs, t; name), systems...)
end

"""
    FixedTranslation(; name, r)

Fixed translation of `frame_b` with respect to `frame_a` with position vector `r` resolved in `frame_a`.

Can be thought of as a massless rod. For a massive rod, see [`BodyShape`](@ref) or [`BodyCylinder`](@ref).
"""
@component function FixedTranslation(; name, r, radius=0.02f0, color = purple, render = true)
    @named frame_a = Frame()
    @named frame_b = Frame()
    @parameters r[1:3]=r [
        description = "position vector from frame_a to frame_b, resolved in frame_a",
    ]
    @parameters begin
        radius = radius, [description = "Radius of the body in animations"]
        color[1:4] = color, [description = "Color of the body in animations (RGBA)"]
        render = render, [description = "Render the component in animations"]
    end
    fa = frame_a.f
    fb = frame_b.f
    taua = frame_a.tau
    taub = frame_b.tau
    eqs = Equation[frame_b.r_0 ~ frame_a.r_0 + resolve1(ori(frame_a), r)
                   ori(frame_b) ~ ori(frame_a)
                   zeros(3) ~ fa + fb
                   zeros(3) ~ taua + taub + cross(r, fb)]
    pars = [r; radius; color; render]
    vars = []
    compose(System(eqs, t, vars, pars; name), frame_a, frame_b)
end

"""
    FixedRotation(; name, r, n, sequence, isroot = false, angle)

Fixed translation followed by a fixed rotation of `frame_b` with respect to `frame_a`

- `r`: Translation vector
- `n`: Axis of rotation, resolved in frame_a
- `angle`: Angle of rotation around `n`, given in radians

To obtain an axis-angle representation of any rotation, see [Conversion between orientation formats](@ref)
"""
@component function FixedRotation(; name, r=[0, 0, 0], n = [1, 0, 0], isroot = false,
                       angle, render=true)
    norm(n) ≈ 1 || error("n must be a unit vector")
    @named frame_a = Frame()
    @named frame_b = Frame()
    @parameters r[1:3]=r [
        description = "position vector from frame_a to frame_b, resolved in frame_a",
    ]
    @parameters n[1:3]=n [
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
    @parameters begin
        render = render, [description = "Render the component in animations"]
    end

    pars = [r; n; angle; render]

    fa = frame_a.f
    fb = frame_b.f
    taua = frame_a.tau
    taub = frame_b.tau

    # Relationships between quantities of frame_a and frame_b 

    if isroot
        Rrel = planar_rotation(n, angle, 0)
        eqs = [ori(frame_b) ~ absolute_rotation(frame_a, Rrel);
               zeros(3) ~ fa + resolve1(Rrel, fb);
               zeros(3) ~ taua + resolve1(Rrel, taub) - cross(r,
                                                                fa)]
    else
        Rrel_inv = planar_rotation(n, -angle, 0)
        eqs = [ori(frame_a) ~ absolute_rotation(frame_b, Rrel_inv);
               zeros(3) ~ fb + resolve1(Rrel_inv, fa);
               zeros(3) ~ taub + resolve1(Rrel_inv, taua) +
                           cross(resolve1(Rrel_inv, r), fb)]
    end
    append!(eqs, frame_b.r_0 ~ frame_a.r_0 + resolve1(frame_a, r))

    compose(System(eqs, t, [], pars; name), frame_a, frame_b)
end

"""
    Body(; name, m = 1, r_cm, isroot = false, phi0 = zeros(3), phid0 = zeros(3), r_0 = zeros(3), state_priority = 2, quat = false, sparse_I = false)

Representing a body with 3 translational and 3 rotational degrees-of-freedom.

This component has a single frame, `frame_a`. To represent bodies with more than one frame, see [`BodyShape`](@ref), [`BodyCylinder`](@ref), [`BodyBox`](@ref). The inertia tensor is defined with respect to a coordinate system that is parallel to `frame_a` with the origin at the center of mass of the body.

# Performance optimization
- `sparse_I`: If `true`, the zero elements of the inerita matrix are considered "structurally zero", and this fact is used to optimize performance. When this option is enabled, the elements of the inertia matrix that were zero when the component was created cannot changed without reinstantiating the component. This performance optimization may be useful, e.g., when the inertia matrix is known to be diagonal.

# Parameters
- `m`: Mass
- `r_cm`: Vector from `frame_a` to center of mass, resolved in `frame_a`
- `I_11, I_22, I_33, I_21, I_31, I_32`: Inertia-matrix elements
- `isroot`: Indicate whether this component is the root of the system, useful when there are no joints in the model.
- `phi0`: Initial orientation, only applicable if `isroot = true` and `quat = false`
- `phid0`: Initial angular velocity

# Variables
- `r_0`: Position vector from origin of world frame to origin of `frame_a`
- `v_0`: Absolute velocity of `frame_a`, resolved in world frame (= D(r_0))
- `a_0`: Absolute acceleration of `frame_a` resolved in world frame (= D(v_0))

# Rendering options
- `radius`: Radius of the joint in animations
- `cylinder_radius`: Radius of the cylinder from frame to COM in animations (only drawn if `r_cm` is non-zero). Defaults to `radius/2`
- `color`: Color of the joint in animations, a vector of length 4 with values between [0, 1] providing RGBA values
"""
@component function Body(; name, m = 1, r_cm = [0, 0, 0],
              I_11 = 0.001,
              I_22 = 0.001,
              I_33 = 0.001,
              I_21 = 0,
              I_31 = 0,
              I_32 = 0,
              sparse_I = false,
              isroot = false,
              state = false,
              sequence = [1,2,3],
              quat = false,
              phi0 = state || isroot ? zeros(3) : nothing,
              phid0 = state || isroot ? zeros(3) : nothing,
              r_0 = state || isroot ? zeros(3) : nothing,
              v_0 = state || isroot ? zeros(3) : nothing,
              w_a = (state || isroot) && quat ? zeros(3) : nothing,
              radius = 0.05,
              cylinder_radius = radius/2,
              length_fraction = 1,
              air_resistance = 0.0,
              color = [1,0,0,1],
              state_priority = 2,
              render = true,
              )
    if state
        # @warn "Make the body have state variables by using isroot=true rather than state=true"
        isroot = true
    end
    @variables r_0(t)[1:3]=r_0 [
        state_priority = state_priority+isroot,
        description = "Position vector from origin of world frame to origin of frame_a",
    ]
    @variables v_0(t)[1:3]=v_0 [ 
        state_priority = state_priority+isroot,
        description = "Absolute velocity of frame_a, resolved in world frame (= D(r_0))",
    ]
    @variables a_0(t)[1:3] [ 
        description = "Absolute acceleration of frame_a resolved in world frame (= D(v_0))",
    ]
    @variables g_0(t)[1:3] [ description = "gravity acceleration"]
    @variables w_a(t)[1:3]=w_a [ 
        state_priority = isroot ? quat ? state_priority : -1 : 0,
        description = "Absolute angular velocity of frame_a resolved in frame_a",
    ]
    @variables z_a(t)[1:3] [ 
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
    @parameters cylinder_radius=cylinder_radius [
        description = "Radius of the cylinder from frame to COM in animations",
    ]
    @parameters color[1:4] = color [description = "Color of the body in animations (RGBA)"]
    @parameters length_fraction=length_fraction, [description = "Fraction of the length of the body that is the cylinder from frame to COM in animations"]
    @parameters render = render [description = "Render the component in animations"]
    # @parameters I[1:3, 1:3]=I [description="inertia tensor"]

    if sparse_I
        Isparsity = sparse(.!isequal.(0, [I_11 I_21 I_31; I_21 I_22 I_32; I_31 I_32 I_33]))
    end

    @parameters I_11=I_11 [description = "Element (1,1) of inertia tensor"]
    @parameters I_22=I_22 [description = "Element (2,2) of inertia tensor"]
    @parameters I_33=I_33 [description = "Element (3,3) of inertia tensor"]
    @parameters I_21=I_21 [description = "Element (2,1) of inertia tensor"]
    @parameters I_31=I_31 [description = "Element (3,1) of inertia tensor"]
    @parameters I_32=I_32 [description = "Element (3,2) of inertia tensor"]

    I = [I_11 I_21 I_31; I_21 I_22 I_32; I_31 I_32 I_33]
    if sparse_I
        I = I.*Isparsity
    end

    # r_0, v_0, a_0, g_0, w_a, z_a, r_cm = collect.((r_0, v_0, a_0, g_0, w_a, z_a, r_cm))

    # DRa = D(Ra)

    dvs = [r_0;v_0;a_0;g_0;w_a;z_a;]

    eqs = if isroot # isRoot
        
        if quat
            @named frame_a = Frame(varw=false)
            Ra = ori(frame_a, false)
            qeeqs = nonunit_quaternion_equations(Ra, w_a)
        else
            @named frame_a = Frame(varw=true)
            Ra = ori(frame_a, true)
            @variables phi(t)[1:3]=phi0 [state_priority = 10, description = "Euler angles"]
            @variables phid(t)[1:3]=phid0 [state_priority = 10]
            @variables phidd(t)[1:3] [state_priority = 0]
            # phi, phid, phidd = collect.((phi, phid, phidd))
            ar = axes_rotations(sequence, phi, phid)
            Equation[
                    phid .~ D.(phi)
                    phidd .~ D.(phid)
                    Ra.w .~ ar.w
                    collect(w_a .~ (angular_velocity2(ar)))
                        # w_a .~ ar.w # This one for most systems
                    Ra ~ ar
                    ]
        end
    else
        @named frame_a = Frame()
        Ra = ori(frame_a)
        # This branch has never proven to be incorrect
        # This equation is defined here and not in the Rotation component since the branch above might use another equation
        (w_a ~ angular_velocity2(Ra))
        # (w_a ~ get_w(Ra))
    end

    eqs = Equation[eqs;
           (r_0 ~ frame_a.r_0)
           (g_0 ~ gravity_acceleration(frame_a.r_0 .+ resolve1(Ra, r_cm)))
           (v_0 ~ D(r_0))
           (a_0 ~ D(v_0))
           (z_a ~ D(w_a))
           if air_resistance > 0
                (frame_a.f ~ m * (resolve2(Ra, a_0 - g_0 + air_resistance*_norm(v_0)*v_0) + cross(z_a, r_cm) +
                                        cross(w_a, cross(w_a, r_cm))))
           else
                (frame_a.f ~ m * (resolve2(Ra, a_0 - g_0) + cross(z_a, r_cm) +
                                        cross(w_a, cross(w_a, r_cm))))
           end
           (frame_a.tau ~ I * z_a + cross(w_a, collect(I * w_a)) + cross(r_cm, frame_a.f))]

    # pars = [m;r_cm;radius;I_11;I_22;I_33;I_21;I_31;I_32;color]
    
    sys = System(eqs, t; name=:nothing, metadata = Dict(IsRoot => isroot), systems = [frame_a])
    add_params(sys, [radius; cylinder_radius; color; length_fraction; render]; name)
end


"""
    BodyShape(; name, m = 1, r, kwargs...)

The `BodyShape` component is similar to a [`Body`](@ref), but it has two frames and a vector `r` that describes the translation between them, while the body has a single frame only.

- `r`: Vector from `frame_a` to `frame_b` resolved in `frame_a`
- All `kwargs` are passed to the internal `Body` component.
- `shapefile`: A path::String to a CAD model that can be imported by MeshIO for 3D rendering. If none is provided, a cylinder shape is rendered.

See also [`BodyCylinder`](@ref) and [`BodyBox`](@ref) for body components with predefined shapes and automatically computed inertial properties based on geometry and density.
"""
@component function BodyShape(; name, state=false, m = 1, r = [0, 0, 0], r_cm = 0.5*r,
    radius = 0.08, color=purple, shapefile="", shape_transform = I(4), shape_scale = 1,
    height = 0.1_norm(r), width = height, shape = "cylinder",
    I_11 = 0.001,
    I_22 = 0.001,
    I_33 = 0.001,
    I_21 = 0,
    I_31 = 0,
    I_32 = 0,
    kwargs...
    )
    pars = @parameters begin
        m = m, [description = "mass"]
        I_11=I_11, [description = "Element (1,1) of inertia tensor"]
        I_22=I_22, [description = "Element (2,2) of inertia tensor"]
        I_33=I_33, [description = "Element (3,3) of inertia tensor"]
        I_21=I_21, [description = "Element (2,1) of inertia tensor"]
        I_31=I_31, [description = "Element (3,1) of inertia tensor"]
        I_32=I_32, [description = "Element (3,2) of inertia tensor"]
    end
    systems = @named begin
        translation = FixedTranslation(r = r, render=false)
        translation_cm = FixedTranslation(r = r_cm, render=false)
        body = Body(; m, r_cm, I_11, I_22, I_33, I_21, I_31, I_32, render=false, kwargs...)
        frame_a = Frame()
        frame_b = Frame()
        frame_cm = Frame()
    end

    # NOTE: these parameters should be defined before the `systems` block above, but due to bugs in MTK/JSC with higher-order array parameters we cannot do that. We still define the parameters so that they are available to make animations
    @variables r_0(t)[1:3] [
        state_priority = 2,
        description = "Position vector from origin of world frame to origin of frame_a",
    ]
    @variables v_0(t)[1:3] [ guess=0,
        state_priority = 2,
        description = "Absolute velocity of frame_a, resolved in world frame (= D(r_0))",
    ]
    @variables a_0(t)[1:3] [ guess=0,
        description = "Absolute acceleration of frame_a resolved in world frame (= D(v_0))",
    ]

    shapecode = encode(shapefile)
    shape = encode(shape)
    more_pars = @parameters begin
        r[1:3]=r, [
            description = "Vector from frame_a to frame_b resolved in frame_a",
        ]
        radius = radius, [description = "Radius of the body in animations"]
        color[1:4] = color, [description = "Color of the body in animations"]
        shapefile[1:length(shapecode)] = shapecode
        shape_transform[1:16] = vec(shape_transform)
        shape_scale = shape_scale
        width = width, [description = """Width of the body in animations (if shape = "box")"""]
        height = height, [description = """Height of the body in animations (if shape = "box")"""]
        shape[1:length(shape)] = shape
    end

    pars = collect_all([pars; more_pars])

    eqs = [r_0 ~ frame_a.r_0
           v_0 ~ D(r_0)
           a_0 ~ D(v_0)
           connect(frame_a, translation.frame_a, translation_cm.frame_a)
           connect(frame_b, translation.frame_b)
           connect(frame_a, body.frame_a)
           connect(frame_cm, translation_cm.frame_b)
           ]
    System(eqs, t, [r_0; v_0; a_0], pars; name, systems)
end


"""
    Rope(; name, l = 1, n = 10, m = 1, c = 0, d = 0, kwargs)

Model a rope (string / cable) of length `l` and mass `m`.

The rope is modeled as a series of `n` links, each connected by a [`Spherical`](@ref) joint. The links are either fixed in length (default, modeled using [`BodyShape`](@ref)) or flexible, in which case they are modeled as a [`Translational.Spring`](@ref) and [`Translational.Damper`](@ref) in parallel with a [`Prismatic`](@ref) joint with a [`Body`](@ref) adding mass to the center of the link segment. The flexibility is controlled by the parameters `c` and `d`, which are the stiffness and damping coefficients of the spring and damper, respectively. The default values are `c = 0` and `d = 0`, which corresponds to a stiff rope.


- `l`: Unstretched length of rope
- `n`: Number of links used to model the rope. For accurate approximations to continuously flexible ropes, a large number may be required.
- `m`: The total mass of the rope. Each rope segment will have mass `m / n`.
- `c`: The equivalent stiffness of the rope, i.e., the rope will act like a spring with stiffness `c`. 
- `d`: The equivalent damping in the stretching direction of the rope, i.e., the taught rope will act like a damper with damping `d`.
- `d_joint`: Viscous damping in the joints between the links. A positive value makes the rope dissipate energy while flexing (as opposed to the damping `d` which dissipats energy due to stretching).
- `dir`: A vector of norm 1 indicating the initial direction of the rope.

## Damping
There are three different methods of adding damping to the rope:
- Damping in the stretching direction of the rope, controlled by the parameter `d`.
- Damping in flexing of the rope, modeled as viscous friction in the joints between the links, controlled by the parameter `d_joint`.
- Air resistance to the rope moving through the air, controlled by the parameter `air_resistance`. This damping is quadratic in the velocity (``f_d ~ -||v||v``) of each link relative to the world frame.

## Rendering
- `color = [255, 219, 120, 255]./255`
- `radius = 0.05f0`
- `jointradius=0`
- `jointcolor=color`
"""
function Rope(; name, l = 1, dir = [0,-1, 0], n = 10, m = 1, c = 0, d=0, air_resistance=0, d_joint = 0, cutspherical = false, cutprismatic=false, color = [255, 219, 120, 255]./255, radius = 0.05f0, jointradius=0, jointcolor=color, kwargs...)

    @assert n >= 1
    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
    end
    dir = dir / norm(dir)

    li = l / n # Segment length
    mi = m / n # Segment mass

    joints = [Spherical(name=Symbol("joint_$i"),
        # isroot=!(cutspherical && i == 1),
        isroot=true,
        iscut = cutspherical && i == 1,
        # state=!(cutspherical && i == 1),
        state=true,
        color=jointcolor,
        radius=jointradius,
        d = d_joint) for i = 1:n+1]

    eqs = [
        connect(frame_a, joints[1].frame_a)
        connect(frame_b, joints[end].frame_b)
    ]
    
    if c > 0
        ModelingToolkitStandardLibrary.@symcheck m > 0 || error("A rope with flexibility (c > 0) requires a non-zero mass (m > 0)")
        ci = n * c # Segment stiffness
        di = n * d # Segment damping
        springs = [Translational.Spring(c = ci, s_rel0=li, name=Symbol("link_$i")) for i = 1:n]
        dampers = [Translational.Damper(d = di, name=Symbol("damping_$i")) for i = 1:n]
        masses = [Body(; m = mi, name=Symbol("mass_$i"), isroot=false, r_cm = li/2*dir, air_resistance, color=0.9*color) for i = 1:n]
        links = [Prismatic(; n = dir, s0 = li, name=Symbol("flexibility_$i"), axisflange=true, color, radius, iscut = cutprismatic && i == 1) for i = 1:n]
        for i = 1:n
            push!(eqs, connect(links[i].support, springs[i].flange_a, dampers[i].flange_a))
            push!(eqs, connect(links[i].axis, springs[i].flange_b, dampers[i].flange_b))
            push!(eqs, connect(links[i].frame_a, masses[i].frame_a))
        end
        links = [links; springs; dampers; masses]
    else
        links = [BodyShape(; m = mi, r = li*dir, name=Symbol("link_$i"), isroot=false, air_resistance, color, radius) for i = 1:n]
    end


    for i = 1:n
        push!(eqs, connect(joints[i].frame_b, links[i].frame_a))
        push!(eqs, connect(links[i].frame_b, joints[i+1].frame_a))
    end

    System(eqs, t; name, systems = [systems; links; joints])
end

# @component function BodyCylinder(; name, m = 1, r = [0.1, 0, 0], r_0 = 0, r_shape=zeros(3), length = _norm(r - r_shape), kwargs...)
#     @parameters begin
#         # r[1:3]=r, [ # MTKs symbolic language is too weak to handle this as a symbolic parameter in from_nxy
#         #     description = "Vector from frame_a to frame_b resolved in frame_a",
#         # ]
#         # r_shape[1:3]=zeros(3), [
#         #     description = "Vector from frame_a to cylinder origin, resolved in frame_a",
#         # ]
#     end
#     r, r_shape = collect.((r, r_shape))
#     @parameters begin
#         dir[1:3] = r - r_shape, [
#             description = "Vector in length direction of cylinder, resolved in frame_a",
#         ]
#         length = _norm(r - r_shape), [
#             description = "Length of cylinder",
#         ]
#         length2 = _norm(r - r_shape), [ # NOTE: strange bug in JSCompiler when both I and r_cm that are parameters if Body depend on the same paramter length. Introducing a dummy parameter with the same value works around the issue. This is not ideal though, since the two parameters must have the same value.
#             description = "Length of cylinder",
#         ]
#         diameter = 1, [#length/5, [
#             description = "Diameter of cylinder",
#         ]
#         inner_diameter = 0, [
#             description = "Inner diameter of cylinder (0 <= inner_diameter <= Diameter)",
#         ]
#         density = 7700, [
#             description = "Density of cylinder (e.g., steel: 7700 .. 7900, wood : 400 .. 800)",
#         ]
#     end
#     # @variables length2(t)
#     # @assert isequal(length, length2)
#     dir = collect(dir) #.|> ParentScope # The ParentScope is required, otherwise JSCompiler thinks that these parameters belong to Body.
#     # length = ParentScope(length)
#     # diameter = ParentScope(diameter)
#     # inner_diameter = ParentScope(inner_diameter)
#     # density = ParentScope(density)

#     radius = diameter/2
#     innerRadius = inner_diameter/2
#     mo = density*pi*length*radius^2
#     mi = density*pi*length*innerRadius^2
#     I22 = (mo*(length^2 + 3*radius^2) - mi*(length^2 + 3*innerRadius^2))/12
#     m = mo - mi
#     R = from_nxy(r, [0, 1, 0]) 
#     r_cm = r_shape + _normalize(dir)*length2/2
#     I = resolve_dyade1(R, Diagonal([(mo*radius^2 - mi*innerRadius^2)/2, I22, I22])) 

#     # r_cm = ParentScope.(r_cm)
#     # I = ParentScope.(I)
#     # m = ParentScope(m)

#     @variables begin
#         r_0(t)[1:3]=r_0, [
#             state_priority = 2,
#             description = "Position vector from origin of world frame to origin of frame_a",
#         ]
#         v_0(t)[1:3]=0, [
#             state_priority = 2,
#             description = "Absolute velocity of frame_a, resolved in world frame (= D(r_0))",
#         ]
#         a_0(t)[1:3]=0, [
#             description = "Absolute acceleration of frame_a resolved in world frame (= D(v_0))",
#         ]
#     end

#     systems = @named begin
#         frame_a = Frame()
#         frame_b = Frame()
#         frameTranslation = FixedTranslation(r = r)
#         body = Body(; m, r_cm, I_11 = I[1,1], I_22 = I[2,2], I_33 = I[3,3], I_21 = I[2,1], I_31 = I[3,1], I_32 = I[3,2], kwargs...)

#     end
#     r_0, v_0, a_0 = collect.((r_0, v_0, a_0))

#     eqs = [r_0 .~ collect(frame_a.r_0)
#            v_0 .~ D.(r_0)
#            a_0 .~ D.(v_0)
#            connect(frame_a, frameTranslation.frame_a)
#            connect(frame_b, frameTranslation.frame_b)
#            connect(frame_a, body.frame_a)]

#     # pars = [
#     #     dir; length; diameter; inner_diameter; density
#     # ] 
#     # vars = [r_0; v_0; a_0]
#     # System(eqs, t, vars, pars; name, systems)
#     System(eqs, t; name, systems)
# end

"""
    BodyCylinder(; name, m = 1, r = [0.1, 0, 0], r_shape = [0, 0, 0], dir = r - r_shape, length = _norm(r - r_shape), diameter = 1, inner_diameter = 0, density = 7700, color = purple)

Rigid body with cylinder shape. The mass properties of the body (mass, center of mass, inertia tensor) are computed from the cylinder data. Optionally, the cylinder may be hollow. The two connector frames `frame_a` and `frame_b` are always parallel to each other.

# Parameters
- `r`: (Structural parameter) Vector from `frame_a` to `frame_b` resolved in `frame_a`
- `r_shape`: (Structural parameter) Vector from `frame_a` to cylinder origin, resolved in `frame_a`
- `dir`: Vector in length direction of cylinder, resolved in `frame_a`
- `length`: Length of cylinder
- `diameter`: Diameter of cylinder
- `inner_diameter`: Inner diameter of cylinder (0 <= inner_diameter <= diameter)
- `density`: Density of cylinder [kg/m³] (e.g., steel: 7700 .. 7900, wood : 400 .. 800)
- `color`: Color of cylinder in animations

# Variables
- `r_0`: Position vector from origin of world frame to origin of `frame_a`
- `v_0`: Absolute velocity of `frame_a`, resolved in world frame (= D(r_0))
- `a_0`: Absolute acceleration of `frame_a` resolved in world frame (= D(v_0))
"""
@component function BodyCylinder(; name, r = [1, 0, 0], r_shape = [0, 0, 0], isroot = false,
                                 state = false, quat = false, sequence = [1, 2, 3],
                                 diameter = 1, inner_diameter = 0, density = 7700, color = purple)
    pars = @parameters begin
        dir[1:3] = r - r_shape, [description = "Vector in length direction of cylinder, resolved in frame_a"]
        length = _norm(r - r_shape), [description = "Length of cylinder"]
        length2 = _norm(r - r_shape), [description = "Length of cylinder"]  # NOTE: strange bug in JSCompiler workaround
        diameter = diameter, [description = "Diameter of cylinder"]
        inner_diameter = inner_diameter, [description = "Inner diameter of cylinder (0 <= inner_diameter <= diameter)"]
        density = density, [description = "Density of cylinder (e.g., steel: 7700 .. 7900, wood : 400 .. 800) [kg/m³]"]
        color[1:4] = color, [description = "Color of cylinder in animations"]
    end

    # Calculations from begin blocks
    radius = diameter/2
    innerRadius = inner_diameter/2
    mo = density*pi*length*radius^2
    mi = density*pi*length*innerRadius^2
    I22 = (mo*(length^2 + 3*radius^2) - mi*(length^2 + 3*innerRadius^2))/12
    m = mo - mi
    R = from_nxy(r, [0, 1, 0])
    r_cm = r_shape .+ _normalize(dir)*length2/2
    I = resolve_dyade1(R, Diagonal([(mo*radius^2 - mi*innerRadius^2)/2, I22, I22]))

    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
        translation = FixedTranslation(r = r)
        body = Body(; m, r_cm, I_11 = I[1,1], I_22 = I[2,2], I_33 = I[3,3], I_21 = I[2,1], I_31 = I[3,1], I_32 = I[3,2], state, quat, isroot, sequence, sparse_I=true)
    end

    vars = @variables begin
        r_0(t)[1:3], [state_priority = 2, description = "Position vector from origin of world frame to origin of frame_a"]
        v_0(t)[1:3], [state_priority = 2, description = "Absolute velocity of frame_a, resolved in world frame (= D(r_0))"]
        a_0(t)[1:3], [description = "Absolute acceleration of frame_a resolved in world frame (= D(v_0))"]
    end

    # Additional calculations from second begin block

    equations = Equation[
        r_0 ~ frame_a.r_0
        v_0 ~ D(r_0)
        a_0 ~ D(v_0)
        connect(frame_a, translation.frame_a)
        connect(frame_b, translation.frame_b)
        connect(frame_a, body.frame_a)
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    BodyBox(; name, m = 1, r = [1, 0, 0], r_shape = [0, 0, 0], width_dir = [0,1,0])

Rigid body with box shape. The mass properties of the body (mass, center of mass, inertia tensor) are computed from the box data. Optionally, the box may be hollow. The (outer) box shape is used in the animation, the hollow part is not shown in the animation. The two connector frames `frame_a` and `frame_b` are always parallel to each other.

# Parameters
- `r`: (structural parameter) Vector from `frame_a` to `frame_b` resolved in `frame_a`
- `r_shape`: (structural parameter) Vector from `frame_a` to box origin, resolved in `frame_a`
- `width_dir`: (structural parameter) Vector in width direction of box, resolved in `frame_a`
- `length_dir`: (structural parameter) Vector in length direction of box, resolved in `frame_a`
- `length`: (structural parameter) Length of box
- `width = 0.3length`: Width of box
- `height = width`: Height of box
- `inner_width`: Width of inner box surface (0 <= inner_width <= width)
- `inner_height`: Height of inner box surface (0 <= inner_height <= height)
- `density = 7700`: Density of cylinder (e.g., steel: 7700 .. 7900, wood : 400 .. 800)
- `color`: Color of box in animations
"""
@component function BodyBox(; name, r = [1, 0, 0], r_shape = [0, 0, 0], width_dir = [0, 1, 0],
                            length_dir = _normalize(r - r_shape), length = _norm(r - r_shape),
                            isroot = false, state = false, quat = false,
                            width = 0.3*length, height = width, inner_width = 0, inner_height = inner_width,
                            density = 7700, color = purple)
    # Validation from first begin block
    iszero(r_shape) || error("non-zero r_shape not supported")

    pars = @parameters begin
        # NOTE: these are workarounds to allow rendering of this component. Unfortunately, MTK/JSCompiler cannot handle parameter arrays well enough to let these be actual parameters
        render_r[1:3] = r, [description = "For internal use only"]
        render_r_shape[1:3] = r_shape, [description = "For internal use only"]
        render_length = length, [description = "For internal use only"]
        render_length_dir[1:3] = length_dir, [description = "For internal use only"]
        render_width_dir[1:3] = width_dir, [description = "For internal use only"]

        width = width, [description = "Width of box"]
        height = height, [description = "Height of box"]
        inner_width = inner_width, [description = "Width of inner box surface (0 <= inner_width <= width)"]
        inner_height = inner_height, [description = "Height of inner box surface (0 <= inner_height <= height)"]
        density = density, [description = "Density of cylinder (e.g., steel: 7700 .. 7900, wood : 400 .. 800)"]
        color[1:4] = color, [description = "Color of box in animations"]
    end

    mo = density*length*width*height
    mi = density*length*inner_width*inner_height
    m = mo - mi
    R = from_nxy(r, width_dir)
    r_cm = r_shape + _normalize(length_dir)*length/2

    I11 = mo*(width^2 + height^2) - mi*(inner_width^2 + inner_height^2)
    I22 = mo*(length^2 + height^2) - mi*(length^2 + inner_height^2)
    I33 = mo*(length^2 + width^2) - mi*(length^2 + inner_width^2)
    I = resolve_dyade1(R, Diagonal([I11, I22, I33] ./ 12))

    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
        translation = FixedTranslation(r = r)
        body = Body(; m, r_cm, I_11 = I[1,1], I_22 = I[2,2], I_33 = I[3,3], I_21 = I[2,1], I_31 = I[3,1], I_32 = I[3,2], state, quat, isroot, sparse_I = true)
    end

    vars = @variables begin
        r_0(t)[1:3], [state_priority = 2, description = "Position vector from origin of world frame to origin of frame_a"]
        v_0(t)[1:3], [state_priority = 2, description = "Absolute velocity of frame_a, resolved in world frame (= D(r_0))"]
        a_0(t)[1:3], [description = "Absolute acceleration of frame_a resolved in world frame (= D(v_0))"]
    end

    equations = Equation[
        r_0 ~ frame_a.r_0
        v_0 ~ D(r_0)
        a_0 ~ D(v_0)
        connect(frame_a, translation.frame_a)
        connect(frame_b, translation.frame_b)
        connect(frame_a, body.frame_a)
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    CraigBampton(; name, boundary_positions, M_BB, M_BI, M_II, K_BB, K_II, ζ, ω)

A flexible body component using Craig-Bampton reduced-order modeling. This allows
including reduced flexible body models derived from external FEM tools.

The Craig-Bampton method reduces a FEM model by partitioning degrees of freedom into:
- **Boundary DOFs (q_B)**: Interface nodes that connect to other multibody components
- **Internal DOFs**: Reduced to modal coordinates η

The reduced equations of motion are:
```
[M_BB  M_BI] [q̈_B]   [C_BB  0  ] [q̇_B]   [K_BB  0  ] [q_B]   [f_B]
[M_IB  M_II] [η̈  ] + [0    C_II] [η̇  ] + [0    K_II] [η  ] = [f_I]
```

where `C_II = diag(2*ζ_i*ω_i*m_ii)` provides modal damping.

# Arguments
- `boundary_positions`: Vector of 3-element vectors, positions of boundary nodes in body reference frame.
  The first boundary corresponds to `frame_a`, additional boundaries create `frame_b`, `frame_c`, etc.
- `M_BB`: Boundary mass matrix `[6*n_boundaries × 6*n_boundaries]`
- `M_BI`: Boundary-internal coupling mass matrix `[6*n_boundaries × n_modes]`
- `M_II`: Internal/modal mass matrix `[n_modes × n_modes]` (often diagonal for mass-normalized modes)
- `K_BB`: Boundary stiffness matrix `[6*n_boundaries × 6*n_boundaries]`
- `K_II`: Internal/modal stiffness matrix `[n_modes × n_modes]` (diagonal = ω² for mass-normalized modes)
- `ζ`: Vector of modal damping ratios `[n_modes]`
- `ω`: Vector of modal frequencies in rad/s `[n_modes]`, used to compute damping matrix

Each boundary has 6 DOFs: 3 translations + 3 rotations (small angle approximation).

# Connectors
- `frame_a`: Reference boundary frame (first boundary position)
- `frame_b`, `frame_c`, ...: Additional boundary frames

# State Variables
- `η[1:n_modes]`: Modal amplitudes
- `η̇[1:n_modes]`: Modal velocities
- `q_B[7:6*n_boundaries]`: Boundary deformations for boundaries 2..n (relative to frame_a)
- `q̇_B[7:6*n_boundaries]`: Boundary deformation velocities

The first boundary (frame_a) serves as the reference with zero deformation. Body motion
comes from how frame_a is connected externally (e.g., to joints or world).

# Example
```julia
# Simple 2-boundary flexible beam (10 modes)
n_modes = 10
n_boundaries = 2
n_bdof = 12  # 6 DOFs per boundary

# Mass-normalized mode shapes give M_II = I, K_II = diag(ω²)
M_BB = [...] # From FEM export
M_BI = [...] # From FEM export
M_II = diagm(ones(n_modes))
K_II = diagm(ω.^2)
ζ = 0.02 * ones(n_modes)  # 2% damping all modes
ω = [10, 25, 40, ...]     # Modal frequencies (rad/s)

@named flexible_beam = CraigBampton(
    boundary_positions = [[0,0,0], [1,0,0]],  # Beam ends
    M_BB = M_BB, M_BI = M_BI, M_II = M_II,
    K_BB = zeros(n_bdof, n_bdof), K_II = K_II,
    ζ = ζ, ω = ω
)

# Connect to world at frame_a, attach payload at frame_b
connections = [
    connect(world.frame_b, flexible_beam.frame_a)
    connect(flexible_beam.frame_b, payload.frame_a)
]
```

See also [`Body`](@ref), [`BodyShape`](@ref) for rigid body components.
"""
@component function CraigBampton(;
    name,
    boundary_positions,  # Vector of 3-vectors: [[x1,y1,z1], [x2,y2,z2], ...]
    M_BB,               # Boundary mass matrix [6*n_b × 6*n_b]
    M_BI,               # Boundary-internal coupling [6*n_b × n_modes]
    M_II,               # Internal modal mass [n_modes × n_modes]
    K_BB,               # Boundary stiffness [6*n_b × 6*n_b]
    K_II,               # Internal modal stiffness [n_modes × n_modes]
    ζ,                  # Modal damping ratios [n_modes]
    ω,                  # Modal frequencies for computing damping [n_modes]
    render = true,
    color = [0.6, 0.6, 0.8, 1.0],
)
    n_boundaries = length(boundary_positions)
    n_modes = size(M_II, 1)
    n_bdof = 6 * n_boundaries  # 6 DOFs per boundary (3 trans + 3 rot)

    @assert size(M_BB) == (n_bdof, n_bdof) "M_BB must be $(n_bdof)×$(n_bdof)"
    @assert size(M_BI) == (n_bdof, n_modes) "M_BI must be $(n_bdof)×$(n_modes)"
    @assert size(M_II) == (n_modes, n_modes) "M_II must be $(n_modes)×$(n_modes)"
    @assert size(K_BB) == (n_bdof, n_bdof) "K_BB must be $(n_bdof)×$(n_bdof)"
    @assert size(K_II) == (n_modes, n_modes) "K_II must be $(n_modes)×$(n_modes)"
    @assert length(ζ) == n_modes "ζ must have length $(n_modes)"
    @assert length(ω) == n_modes "ω must have length $(n_modes)"

    # Compute modal damping matrix C_II = diag(2*ζ*ω*m_ii)
    C_II_diag = 2 .* ζ .* ω .* diag(M_II)

    # Create frames - frame_a is reference, then frame_b, frame_c, etc.
    @named frame_a = Frame()
    frames = Any[frame_a]
    if n_boundaries >= 2
        @named frame_b = Frame()
        push!(frames, frame_b)
    end
    if n_boundaries >= 3
        @named frame_c = Frame()
        push!(frames, frame_c)
    end
    if n_boundaries >= 4
        @named frame_d = Frame()
        push!(frames, frame_d)
    end
    if n_boundaries >= 5
        @named frame_e = Frame()
        push!(frames, frame_e)
    end
    if n_boundaries >= 6
        @named frame_f = Frame()
        push!(frames, frame_f)
    end
    if n_boundaries > 6
        error("CraigBampton currently supports at most 6 boundary frames")
    end

    # Store matrices as parameters
    pars = @parameters begin
        render = render, [description = "Render the component in animations"]
        color[1:4] = color, [description = "Color in animations (RGBA)"]
    end

    # State variables - modal coordinates
    @variables begin
        η(t)[1:n_modes] = zeros(n_modes), [description = "Modal amplitudes"]
        η̇(t)[1:n_modes] = zeros(n_modes), [description = "Modal velocities"]
    end

    # Boundary deformation variables for boundaries 2..n (boundary 1 = frame_a is reference)
    # Each boundary has 6 DOFs: [ux, uy, uz, θx, θy, θz]
    n_flex_dof = 6 * (n_boundaries - 1)  # DOFs for flexible deformations (excluding reference)

    @variables begin
        q_flex(t)[1:n_flex_dof] = zeros(n_flex_dof), [description = "Boundary deformations (boundaries 2..n)"]
        q̇_flex(t)[1:n_flex_dof] = zeros(n_flex_dof), [description = "Boundary deformation velocities"]
    end

    # Reference frame orientation
    R_a = ori(frame_a)

    # Build equation arrays
    eqs = Equation[]

    # Kinematic equations for states
    append!(eqs, D.(collect(η)) .~ collect(η̇))
    if n_flex_dof > 0
        append!(eqs, D.(collect(q_flex)) .~ collect(q̇_flex))
    end

    # Build full q_B vector: [zeros(6); q_flex] (frame_a deformation is zero)
    # q_B represents deformations in the body frame (frame_a coordinates)
    q_B_full = vcat(zeros(Num, 6), collect(q_flex))
    q̇_B_full = vcat(zeros(Num, 6), collect(q̇_flex))

    # Second derivatives for dynamics
    η_ddot = D.(collect(η̇))
    q̈_B_full = vcat(zeros(Num, 6), D.(collect(q̇_flex)))

    # Modal dynamics equation:
    # M_II * η̈ + C_II * η̇ + K_II * η = -M_BI' * q̈_B
    # Note: M_IB = M_BI' (symmetric coupling)
    modal_forcing = -M_BI' * q̈_B_full
    for i in 1:n_modes
        lhs = sum(M_II[i,j] * η_ddot[j] for j in 1:n_modes) +
              C_II_diag[i] * η̇[i] +
              sum(K_II[i,j] * η[j] for j in 1:n_modes)
        push!(eqs, lhs ~ modal_forcing[i])
    end

    # Boundary force equations:
    # The generalized forces at boundaries from CB dynamics:
    # f_B = K_BB * q_B + M_BB * q̈_B + M_BI * η̈
    f_B = K_BB * q_B_full + M_BB * q̈_B_full + M_BI * η_ddot

    # Frame kinematics and force balance for each boundary
    for i in 1:n_boundaries
        frame_i = frames[i]
        p_i = collect(boundary_positions[i])  # Undeformed position in body frame

        # Index into q_B for this boundary's DOFs
        idx_start = 6*(i-1) + 1
        idx_end = 6*i

        if i == 1
            # frame_a is the reference - its deformation is zero
            # Position equation: frame_a position is determined by external connections
            # No additional position constraint needed for frame_a

            # Force balance for frame_a
            f_i = f_B[idx_start:idx_start+2]      # Forces in body frame
            τ_i = f_B[idx_start+3:idx_end]        # Torques in body frame

            # The reaction force at frame_a equals the internal CB forces (resolved to world frame)
            append!(eqs, collect(frame_a.f) .~ resolve1(R_a, f_i))
            append!(eqs, collect(frame_a.tau) .~ resolve1(R_a, τ_i))
        else
            # Other boundaries have position determined by deformation
            flex_idx = 6*(i-2) + 1  # Index into q_flex
            u_i = q_flex[flex_idx:flex_idx+2]      # Translation deformation
            θ_i = q_flex[flex_idx+3:flex_idx+5]    # Rotation deformation (small angles)

            # Position: r_i = r_a + R_a * (p_i + u_i)
            r_rel = p_i + collect(u_i)  # Relative position in body frame
            append!(eqs, collect(frame_i.r_0) .~ collect(frame_a.r_0) + resolve1(R_a, r_rel))

            # Orientation: For small rotations, use same orientation as frame_a
            # (Full small rotation: R_i = R_a * (I + skew(θ_i)), but for simplicity use R_a)
            append!(eqs, ori(frame_i) ~ R_a)

            # Force balance for this boundary
            f_i = f_B[idx_start:idx_start+2]      # Forces in body frame
            τ_i = f_B[idx_start+3:idx_end]        # Torques in body frame

            # Include moment from offset for force balance
            append!(eqs, collect(frame_i.f) .~ resolve1(R_a, f_i))
            append!(eqs, collect(frame_i.tau) .~ resolve1(R_a, τ_i))
        end
    end

    System(eqs, t, [η; η̇; q_flex; q̇_flex], pars; name, systems=frames)
end

#=
CB TODO
The "simplified" version has these limitations:

  1. frame_a is the fixed reference - body motion comes entirely from external connections to frame_a
  2. No rigid body dynamics - the CB component doesn't have its own mass/inertia for overall body motion
  3. No coupling between rigid motion and flexible modes - Coriolis/centrifugal effects are ignored

  To lift these restrictions and create a floating frame of reference implementation:

  Key Changes Needed

  1. Add Rigid Body States

  Like the Body component with isroot=true:

  @variables begin
      # Body reference frame position/velocity in world
      r_0(t)[1:3]      # Position of body frame origin
      v_0(t)[1:3]      # Translational velocity

      # Orientation (quaternions or Euler angles)
      Q(t)[1:4]        # Quaternion orientation
      # OR phi(t)[1:3] for Euler angles

      w_a(t)[1:3]      # Angular velocity in body frame
  end

  2. Extended Mass Matrix

  The full floating-frame CB equations are:

  [M_rr  M_rf  M_rB] [q̈_r]   [h_r]   [f_r]
  [M_fr  M_ff  M_fB] [η̈  ] + [h_f] = [f_f]  
  [M_Br  M_Bf  M_BB] [q̈_B]   [h_B]   [f_B]

  Where:
  - q_r = rigid body DOFs (6: 3 translation + 3 rotation)
  - M_rr = rigid body inertia (6×6) - total mass and inertia tensor
  - M_rf = rigid-flexible coupling - how modes affect rigid body motion
  - h = quadratic velocity terms (Coriolis, centrifugal)

  3. New Input Parameters

  CraigBampton(;
      # ... existing parameters ...

      # Floating frame parameters
      floating = false,       # Enable floating frame
      m_total,               # Total mass of flexible body
      I_cm,                  # Inertia tensor at center of mass
      r_cm,                  # Position of CM in body frame
      M_rf,                  # Rigid-flexible coupling matrix [6 × n_modes]

      # Optional
      isroot = false,        # Whether this body provides root states
      quat = true,           # Use quaternions for orientation
  )

  4. Frame Kinematics Changes

  Currently frame_a is special (reference). With floating frame, all boundaries are equal:

  # Current (simplified):
  frame_i.r_0 = frame_a.r_0 + R_a * (p_i + u_i)

  # Floating frame version:
  frame_i.r_0 = r_0 + R_body * (p_i + u_i)  # All frames relative to body frame

  5. Quadratic Velocity Terms

  For rotating bodies, Coriolis and centrifugal forces couple rigid and flexible motion:

  # Coriolis terms for flexible modes
  h_f = -2 * M_rf' * [0; 0; 0; cross(w_a, w_a)]  # Simplified

  # Centrifugal stiffening (geometric stiffness)
  K_geo = f(w_a, stress)  # Stress-dependent stiffness

  6. Gravity Distribution

  Gravity must act on the distributed mass:

  g_0 = gravity_acceleration(r_0 + R_body * r_cm)
  f_gravity = m_total * g_0  # Applied at CM

=#