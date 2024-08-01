using LinearAlgebra
import ModelingToolkitStandardLibrary

function isroot(sys)
    sys.metadata isa Dict || return false
    get(sys.metadata, :isroot, false)
end

purple = [0.5019608f0,0.0f0,0.5019608f0,1.0f0]

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

"""
    World(; name, render=true)
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
    @parameters n[1:3]=n0 [description = "gravity direction of world"]
    @parameters g=g0 [description = "gravitational acceleration of world"]
    @parameters mu=mu0 [description = "Gravity field constant [m³/s²] (default = field constant of earth)"]
    @parameters render=render
    @parameters point_gravity = point_gravity
    n = Symbolics.scalarize(n)
    O = ori(frame_b)
    eqs = Equation[
        collect(frame_b.r_0) .~ 0;
        O ~ nullrotation()
    ]
    ODESystem(eqs, t, [], [n; g; mu; point_gravity; render]; name, systems = [frame_b])#, defaults=[n => n0; g => g0; mu => mu0])
end

"""
The world component is the root of all multibody models. It is a fixed frame with a parallel gravitational field and a gravity vector specified by the unit direction `world.n` (defaults to [0, -1, 0]) and magnitude `world.g` (defaults to 9.80665).
"""
const world = World(; name = :world)

"Compute the gravity acceleration, resolved in world frame"
function gravity_acceleration(r)
    inner_gravity(GlobalScope(world.point_gravity), GlobalScope(world.mu), GlobalScope(world.g), GlobalScope.(collect(world.n)), collect(r))
end

function inner_gravity(point_gravity, mu, g, n, r)
    # This is slightly inefficient, producing three if statements, one for each array entry. The function registration for array-valued does not work properly so this is a workaround for now. Hitting, among other problems, https://github.com/SciML/ModelingToolkit.jl/issues/2808
    gvp = -(mu/(r'r))*(r/_norm(r))
    gvu = g * n
    ifelse.(point_gravity==true, gvp, gvu)
end


@component function Fixed(; name, r = [0, 0, 0])
    systems = @named begin frame_b = Frame() end
    @parameters begin r[1:3] = r,
                               [
                                   description = "Position vector from world frame to frame_b, resolved in world frame",
                               ] end
    eqs = [collect(frame_b.r_0 .~ r)
           ori(frame_b) ~ nullrotation()]
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
    @variables begin (housing_tau(t)[1:3]), [
                         description = "Torque",
                     ] end
    eqs = [(collect(housing_tau) .~ collect(-n * flange_b.tau));
           (flange_b.phi ~ phi0);
           connect(housing_frame_a, frame_a)]
    compose(ODESystem(eqs, t; name), systems...)
end

"""
    FixedTranslation(; name, r)

Fixed translation of `frame_b` with respect to `frame_a` with position vector `r` resolved in `frame_a`.

Can be thought of as a massless rod. For a massive rod, see [`BodyShape`](@ref) or [`BodyCylinder`](@ref).
"""
@component function FixedTranslation(; name, r, radius=0.02f0, color = purple)
    @named frame_a = Frame()
    @named frame_b = Frame()
    @parameters r[1:3]=collect(r) [
        description = "position vector from frame_a to frame_b, resolved in frame_a",
    ]
    r = collect(r)
    @parameters begin
        radius = radius, [description = "Radius of the body in animations"]
        color[1:4] = color, [description = "Color of the body in animations (RGBA)"]
    end
    fa = frame_a.f |> collect
    fb = frame_b.f |> collect
    taua = frame_a.tau |> collect
    taub = frame_b.tau |> collect
    eqs = Equation[(collect(frame_b.r_0) .~ collect(frame_a.r_0) + resolve1(ori(frame_a), r))
                   (ori(frame_b) ~ ori(frame_a))
                   collect(0 .~ fa + fb)
                   (0 .~ taua + taub + cross(r, fb))]
    pars = [r; radius; color]
    vars = []
    compose(ODESystem(eqs, t, vars, pars; name), frame_a, frame_b)
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
                       angle)
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

    pars = [r; n; angle]

    fa = frame_a.f |> collect
    fb = frame_b.f |> collect
    taua = frame_a.tau |> collect
    taub = frame_b.tau |> collect

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
    eqs = collect(eqs)
    append!(eqs, collect(frame_b.r_0) .~ collect(frame_a.r_0) + resolve1(frame_a, r))

    compose(ODESystem(eqs, t, [], pars; name), frame_a, frame_b)
end

"""
    Body(; name, m = 1, r_cm, isroot = false, phi0 = zeros(3), phid0 = zeros(3), r_0=zeros(3), state_priority = 2, quat=false)

Representing a body with 3 translational and 3 rotational degrees-of-freedom.

This component has a single frame, `frame_a`. To represent bodies with more than one frame, see [`BodyShape`](@ref), [`BodyCylinder`](@ref), [`BodyBox`](@ref).

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
              isroot = false,
              state = false,
              sequence = [1,2,3],
              neg_w = true,
              phi0 = zeros(3),
              phid0 = zeros(3),
              r_0 = 0,
              v_0 = 0,
              w_a = 0,
              radius = 0.05,
              cylinder_radius = radius/2,
              length_fraction = 1,
              air_resistance = 0.0,
              color = [1,0,0,1],
              state_priority = 2,
              quat=false,)
    if state
        # @warn "Make the body have state variables by using isroot=true rather than state=true"
        isroot = true
    end
    @variables r_0(t)[1:3]=r_0 [
        state_priority = state_priority+isroot,
        description = "Position vector from origin of world frame to origin of frame_a",
    ]
    @variables v_0(t)[1:3]=v_0 [guess = 0, 
        state_priority = state_priority+isroot,
        description = "Absolute velocity of frame_a, resolved in world frame (= D(r_0))",
    ]
    @variables a_0(t)[1:3] [guess = 0, 
        state_priority = state_priority+isroot,
        description = "Absolute acceleration of frame_a resolved in world frame (= D(v_0))",
    ]
    @variables g_0(t)[1:3] [guess = 0, description = "gravity acceleration"]
    @variables w_a(t)[1:3]=w_a [guess = 0, 
        state_priority = state_priority-1+2quat*state,
        description = "Absolute angular velocity of frame_a resolved in frame_a",
    ]
    @variables z_a(t)[1:3] [guess = 0, 
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
        
        if quat
            @named frame_a = Frame(varw=false)
            Ra = ori(frame_a, false)
            qeeqs = nonunit_quaternion_equations(Ra, w_a; neg_w)
        else
            @named frame_a = Frame(varw=true)
            Ra = ori(frame_a, true)
            @variables phi(t)[1:3]=phi0 [state_priority = 10, description = "Euler angles"]
            @variables phid(t)[1:3]=phid0 [state_priority = 10]
            @variables phidd(t)[1:3] [state_priority = 0]
            phi, phid, phidd = collect.((phi, phid, phidd))
            ar = axes_rotations(sequence, phi, phid)
            Equation[
                    phid .~ D.(phi)
                    phidd .~ D.(phid)
                    Ra.w .~ ar.w
                    if neg_w
                        # w_a .~ -ar.w # This is required for FreeBody and ThreeSprings tests to pass, but the other one required for harmonic osciallator without joint to pass. FreeBody passes with quat=true so we use that instead
                        collect(w_a .~ -angular_velocity2(ar))
                    else
                        collect(w_a .~ (angular_velocity2(ar)))
                        # w_a .~ ar.w # This one for most systems
                    end
                    Ra ~ ar
                    ]
        end
    else
        @named frame_a = Frame()
        Ra = ori(frame_a)
        # This branch has never proven to be incorrect
        # This equation is defined here and not in the Rotation component since the branch above might use another equation
        collect(w_a .~ angular_velocity2(Ra))
        # collect(w_a .~ get_w(Ra))
    end

    eqs = [eqs;
           collect(r_0 .~ frame_a.r_0)
           collect(g_0 .~ gravity_acceleration(frame_a.r_0 .+ resolve1(Ra, r_cm)))
           collect(v_0 .~ D.(r_0))
           collect(a_0 .~ D.(v_0))
           collect(z_a .~ D.(w_a))
           if air_resistance > 0
                collect(frame_a.f .~ m * (resolve2(Ra, a_0 - g_0 + air_resistance*_norm(v_0)*v_0) + cross(z_a, r_cm) +
                                        cross(w_a, cross(w_a, r_cm))))
           else
                collect(frame_a.f .~ m * (resolve2(Ra, a_0 - g_0) + cross(z_a, r_cm) +
                                        cross(w_a, cross(w_a, r_cm))))
           end
           collect(frame_a.tau .~ I * z_a + cross(w_a, I * w_a) + cross(r_cm, frame_a.f))]

    # pars = [m;r_cm;radius;I_11;I_22;I_33;I_21;I_31;I_32;color]
    
    sys = ODESystem(eqs, t; name=:nothing, metadata = Dict(:isroot => isroot), systems = [frame_a])
    add_params(sys, [radius; cylinder_radius; color; length_fraction]; name)
end


"""
    BodyShape(; name, m = 1, r, kwargs...)

The `BodyShape` component is similar to a [`Body`](@ref), but it has two frames and a vector `r` that describes the translation between them, while the body has a single frame only.

- `r`: Vector from `frame_a` to `frame_b` resolved in `frame_a`
- All `kwargs` are passed to the internal `Body` component.
- `shapefile`: A path::String to a CAD model that can be imported by MeshIO for 3D rendering. If none is provided, a cylinder shape is rendered.

See also [`BodyCylinder`](@ref) and [`BodyBox`](@ref) for body components with predefined shapes and automatically computed inertial properties based on geometry and density.
"""
@component function BodyShape(; name, m = 1, r = [0, 0, 0], r_cm = 0.5*r, r_0 = 0, radius = 0.08, color=purple, shapefile="", kwargs...)
    systems = @named begin
        translation = FixedTranslation(r = r)
        body = Body(; r_cm, r_0, kwargs...)
        frame_a = Frame()
        frame_b = Frame()
    end

    @variables r_0(t)[1:3]=r_0 [
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
    @parameters begin
        r[1:3]=r, [
            description = "Vector from frame_a to frame_b resolved in frame_a",
        ]
        radius = radius, [description = "Radius of the body in animations"]
        color[1:4] = color, [description = "Color of the body in animations"]
        shapefile[1:length(shapecode)] = shapecode
    end


    pars = [r; radius; color; shapefile]

    r_0, v_0, a_0 = collect.((r_0, v_0, a_0))

    eqs = [r_0 .~ collect(frame_a.r_0)
           v_0 .~ D.(r_0)
           a_0 .~ D.(v_0)
           connect(frame_a, translation.frame_a)
           connect(frame_b, translation.frame_b)
           connect(frame_a, body.frame_a)]
    ODESystem(eqs, t, [r_0; v_0; a_0], pars; name, systems)
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

    ODESystem(eqs, t; name, systems = [systems; links; joints])
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
#     # ODESystem(eqs, t, vars, pars; name, systems)
#     ODESystem(eqs, t; name, systems)
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
@mtkmodel BodyCylinder begin

    @structural_parameters begin
        r = [1, 0, 0]
        r_shape = [0, 0, 0]
        isroot = false
        state = false
        quat = false
        sequence = [1,2,3]
        neg_w = true
    end

    @parameters begin
        # r[1:3]=r, [ # MTKs symbolic language is too weak to handle this as a symbolic parameter in from_nxy
        #     description = "Vector from frame_a to frame_b resolved in frame_a",
        # ]
        # r_shape[1:3]=zeros(3), [
        #     description = "Vector from frame_a to cylinder origin, resolved in frame_a",
        # ]
        dir[1:3] = r - r_shape, [
            description = "Vector in length direction of cylinder, resolved in frame_a",
        ]
        length = _norm(r - r_shape), [
            description = "Length of cylinder",
        ]
        length2 = _norm(r - r_shape), [ # NOTE: strange bug in JSCompiler when both I and r_cm that are parameters of Body depend on the same paramter length. Introducing a dummy parameter with the same value works around the issue. This is not ideal though, since the two parameters must have the same value.
            description = "Length of cylinder",
        ]
        diameter = 1, [#length/5, [
            description = "Diameter of cylinder",
        ]
        inner_diameter = 0, [
            description = "Inner diameter of cylinder (0 <= inner_diameter <= diameter)",
        ]
        density = 7700, [
            description = "Density of cylinder (e.g., steel: 7700 .. 7900, wood : 400 .. 800) [kg/m³]",
        ]
        color[1:4] = purple, [description = "Color of cylinder in animations"]
    end
    begin
        radius = diameter/2
        innerRadius = inner_diameter/2
        mo = density*pi*length*radius^2
        mi = density*pi*length*innerRadius^2
        I22 = (mo*(length^2 + 3*radius^2) - mi*(length^2 + 3*innerRadius^2))/12
        m = mo - mi
        R = from_nxy(r, [0, 1, 0]) 
        r_cm = r_shape + _normalize(dir)*length2/2
        I = resolve_dyade1(R, Diagonal([(mo*radius^2 - mi*innerRadius^2)/2, I22, I22])) 
    end

    @variables begin
        r_0(t)[1:3]=0, [
            state_priority = 2,
            description = "Position vector from origin of world frame to origin of frame_a",
        ]
        v_0(t)[1:3]=0, [
            state_priority = 2,
            description = "Absolute velocity of frame_a, resolved in world frame (= D(r_0))",
        ]
        a_0(t)[1:3]=0, [
            description = "Absolute acceleration of frame_a resolved in world frame (= D(v_0))",
        ]
    end
    begin
        r_0 = collect(r_0)
        r_cm = collect(r_cm)
    end
    @components begin
        frame_a = Frame()
        frame_b = Frame()
        translation = FixedTranslation(r = r)
        body = Body(; m, r_cm, I_11 = I[1,1], I_22 = I[2,2], I_33 = I[3,3], I_21 = I[2,1], I_31 = I[3,1], I_32 = I[3,2], state, quat, isroot, sequence, neg_w)
    end

    @equations begin
        r_0[1] ~ ((frame_a.r_0)[1])
        r_0[2] ~ ((frame_a.r_0)[2])
        r_0[3] ~ ((frame_a.r_0)[3])
        v_0[1] ~ D(r_0[1])
        v_0[2] ~ D(r_0[2])
        v_0[3] ~ D(r_0[3])
        a_0[1] ~ D(v_0[1])
        a_0[2] ~ D(v_0[2])
        a_0[3] ~ D(v_0[3])
        connect(frame_a, translation.frame_a)
        connect(frame_b, translation.frame_b)
        connect(frame_a, body.frame_a)
    end
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
@mtkmodel BodyBox begin
    @structural_parameters begin
        r = [1, 0, 0]
        r_shape = [0, 0, 0]
        width_dir = [0,1,0] # https://github.com/SciML/ModelingToolkit.jl/issues/2810
        length_dir = _normalize(r - r_shape)
        length = _norm(r - r_shape)
        isroot = false
        state = false
        quat = false
    end
    begin
        iszero(r_shape) || error("non-zero r_shape not supported")
        width_dir = collect(width_dir)
        length_dir = collect(length_dir)
    end

    @parameters begin
        # r[1:3]=r, [ # MTKs symbolic language is too weak to handle this as a symbolic parameter in from_nxy
        #     description = "Vector from frame_a to frame_b resolved in frame_a",
        # ]
        # r_shape[1:3]=zeros(3), [
        #     description = "Vector from frame_a to box origin, resolved in frame_a",
        # ]
        # length = _norm(r - r_shape), [
        #     description = "Length of box",
        # ]
        # length_dir[1:3] = _norm(r - r_shape), [
        #     description = "Vector in length direction of box, resolved in frame_a",
        # ]

        # width_dir[1:3] = width_dir0, [ 
        #     description = "Vector in width direction of box, resolved in frame_a",
        # ]

        # NOTE: these are workarounds to allow rendering of this component. Unfortunately, MTK/JSCompiler cannot handle parameter arrays well enough to let these be actual parameters
        render_r[1:3]=r, [description="For internal use only"]
        render_r_shape[1:3]=r_shape, [description="For internal use only"]
        render_length = length, [description="For internal use only"]
        render_length_dir[1:3] = length_dir, [description="For internal use only"]
        render_width_dir[1:3] = width_dir, [description="For internal use only"]


        width = 0.3*length, [
            description = "Width of box",
        ]
        height = width, [
            description = "Height of box",
        ]

        inner_width = 0, [
            description = "Width of inner box surface (0 <= inner_width <= width)",
        ]
        inner_height = inner_width, [
            description = "Height of inner box surface (0 <= inner_height <= height)",
        ]
        density = 7700, [
            description = "Density of cylinder (e.g., steel: 7700 .. 7900, wood : 400 .. 800)",
        ]
        color[1:4] = purple, [description = "Color of box in animations"]
    end
    begin
        mo = density*length*width*height
        mi = density*length*inner_width*inner_height
        m = mo - mi
        R = from_nxy(r, width_dir) 
        r_cm = r_shape + _normalize(length_dir)*length/2
        r_cm = collect(r_cm)

        I11 = mo*(width^2 + height^2) - mi*(inner_width^2 + inner_height^2)
        I22 = mo*(length^2 + height^2) - mi*(length^2 + inner_height^2)
        I33 = mo*(length^2 + width^2) - mi*(length^2 + inner_width^2)
        I = resolve_dyade1(R, Diagonal([I11, I22, I33] ./ 12)) 
    end

    @variables begin
        r_0(t)[1:3]=zeros(3), [
            state_priority = 2,
            description = "Position vector from origin of world frame to origin of frame_a",
        ]
        v_0(t)[1:3]=zeros(3), [
            state_priority = 2,
            description = "Absolute velocity of frame_a, resolved in world frame (= D(r_0))",
        ]
        a_0(t)[1:3]=zeros(3), [
            description = "Absolute acceleration of frame_a resolved in world frame (= D(v_0))",
        ]
    end

    @components begin
        frame_a = Frame()
        frame_b = Frame()
        translation = FixedTranslation(r = r)
        body = Body(; m, r_cm, I_11 = I[1,1], I_22 = I[2,2], I_33 = I[3,3], I_21 = I[2,1], I_31 = I[3,1], I_32 = I[3,2], state, quat, isroot)
    end

    @equations begin
        r_0[1] ~ ((frame_a.r_0)[1])
        r_0[2] ~ ((frame_a.r_0)[2])
        r_0[3] ~ ((frame_a.r_0)[3])
        v_0[1] ~ D(r_0[1])
        v_0[2] ~ D(r_0[2])
        v_0[3] ~ D(r_0[3])
        a_0[1] ~ D(v_0[1])
        a_0[2] ~ D(v_0[2])
        a_0[3] ~ D(v_0[3])
        connect(frame_a, translation.frame_a)
        connect(frame_b, translation.frame_b)
        connect(frame_a, body.frame_a)
    end
end

# function BodyBox2(;
#         r = [1, 0, 0],
#         # r_shape = [0, 0, 0],
#         width_dir = [0,1,0],
#         # length_dir = _normalize(r - r_shape),
#         # length = _norm(r - r_shape),

#         # r,
#         r_shape = nothing,
#         length = nothing,
#         length_dir = nothing,
#         # width_dir = nothing,
#         width = nothing,
#         height = nothing,
#         inner_width = nothing,
#         inner_height = nothing,
#         density = nothing,
#         color = nothing,
#         name,
# )

    
#     # @parameters r[1:3]=something(r, [1,0,0]), [
#     #         description = "Vector from frame_a to frame_b resolved in frame_a",
#     #     ]
#     r = collect(r)
#     @parameters r_shape[1:3]=something(r_shape, zeros(3)), [
#             description = "Vector from frame_a to box origin, resolved in frame_a",
#         ]
#     r_shape = collect(r_shape)
#     @parameters length_dir[1:3] = something(length_dir, _norm(r - r_shape)), [
#             description = "Vector in length direction of box, resolved in frame_a",
#         ]
#     length_dir = collect(length_dir)
#     # @parameters width_dir[1:3] = something(width_dir, [0,1,0]), [ 
#     #         description = "Vector in width direction of box, resolved in frame_a",
#     #     ]
#     width_dir = collect(width_dir)
#     @parameters color[1:4] = something(color, purple), [
#             description = "Color of box in animations"
#         ]
#     color = collect(color)
#     pars = @parameters begin
#         length = something(length, _norm(r - r_shape)), [
#             description = "Length of box",
#         ]
#         width = something(width, 0.3*length), [
#             description = "Width of box",
#         ]
#         height = something(height, width), [
#             description = "Height of box",
#         ]
#         inner_width = something(inner_width, 0), [
#             description = "Width of inner box surface (0 <= inner_width <= width)",
#         ]
#         inner_height = something(inner_height, inner_width), [
#             description = "Height of inner box surface (0 <= inner_height <= height)",
#         ]
#         density = something(density, 7700), [
#             description = "Density of box (e.g., steel: 7700 .. 7900, wood : 400 .. 800)",
#         ]
#     end
#     pars = [
#         pars; 
#         # collect(r);
#         collect(r_shape);
#         collect(length_dir);
#         # collect(width_dir);
#         collect(color);
#     ]
#     mo = density*length*width*height
#     mi = density*length*inner_width*inner_height
#     m = mo - mi
#     R = from_nxy(r, width_dir) 
#     r_cm = r_shape + _normalize(length_dir)*length/2
#     r_cm = collect(r_cm)

#     I11 = mo*(width^2 + height^2) - mi*(inner_width^2 + inner_height^2)
#     I22 = mo*(length^2 + height^2) - mi*(length^2 + inner_height^2)
#     I33 = mo*(length^2 + width^2) - mi*(length^2 + inner_width^2)
#     I = resolve_dyade1(R, Diagonal([I11, I22, I33] ./ 12)) 

#     @variables begin
#         r_0(t)[1:3]=zeros(3), [
#             state_priority = 2,
#             description = "Position vector from origin of world frame to origin of frame_a",
#         ]
#         v_0(t)[1:3]=zeros(3), [
#             state_priority = 2,
#             description = "Absolute velocity of frame_a, resolved in world frame (= D(r_0))",
#         ]
#         a_0(t)[1:3]=zeros(3), [
#             description = "Absolute acceleration of frame_a resolved in world frame (= D(v_0))",
#         ]
#     end
#     r_0, v_0, a_0 = collect.((r_0, v_0, a_0))
#     vars = [r_0; v_0; a_0]

#     systems = @named begin
#         frame_a = Frame()
#         frame_b = Frame()
#         translation = FixedTranslation(r = r)
#         body = Body(; m, r_cm, I_11 = I[1,1], I_22 = I[2,2], I_33 = I[3,3], I_21 = I[2,1], I_31 = I[3,1], I_32 = I[3,2])
#     end

#     equations = [
#         r_0[1] ~ ((frame_a.r_0)[1])
#         r_0[2] ~ ((frame_a.r_0)[2])
#         r_0[3] ~ ((frame_a.r_0)[3])
#         v_0[1] ~ D(r_0[1])
#         v_0[2] ~ D(r_0[2])
#         v_0[3] ~ D(r_0[3])
#         a_0[1] ~ D(v_0[1])
#         a_0[2] ~ D(v_0[2])
#         a_0[3] ~ D(v_0[3])
#         connect(frame_a, translation.frame_a)
#         connect(frame_b, translation.frame_b)
#         connect(frame_a, body.frame_a)
#     ]
#     ODESystem(equations, t, vars, pars; name, systems)
# end
