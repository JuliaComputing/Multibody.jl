import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica as TP

import ModelingToolkitStandardLibrary.Blocks

"""
    BasicTorque(; name, resolve_frame = :world)

Low-level torque component used to build [`Torque`](@ref)
"""
function BasicTorque(; name, resolve_frame = :world)
    @named ptf = PartialTwoFrames()
    @named torque = Blocks.RealInput(; nin = 3)
    @unpack frame_a, frame_b = ptf
    @variables begin
        (r_0(t)[1:3]),
        [
            description = "Position vector from origin of frame_a to origin of frame_b resolved in world frame",
        ]
        (t_b_0(t)[1:3]), [
            description = "frame_b.tau resolved in world frame"]
    end

    eqs = [r_0 ~ frame_b.r_0 - frame_a.r_0
           frame_a.f ~ zeros(3)
           frame_b.f ~ zeros(3)
           # torque balance
           zeros(3) ~ frame_a.tau + resolve2(frame_a, t_b_0)]

    if resolve_frame == :frame_a
        append!(eqs,
                [t_b_0 ~ -resolve1(frame_a, torque.u)
                 frame_b.tau ~ resolve2(frame_b, t_b_0)])
    elseif resolve_frame == :frame_b
        append!(eqs,
                [t_b_0 ~ -resolve1(frame_b, torque.u)
                 frame_b.tau ~ -torque.u])
    elseif resolve_frame == :world
        append!(eqs,
                [t_b_0 ~ -torque.u
                 frame_b.tau ~ resolve2(frame_b, t_b_0)])
    else
        error("Unknown value of argument resolve_frame")
        append!(eqs, [t_b_0 ~ zeros(3)
                      frame_b.tau ~ zeros(3)])
    end

    extend(System(eqs, t, name = name, systems = [torque]), ptf)
end

"""
    Torque(; name, resolve_frame = :frame_b)

Torque acting between two frames, defined by 3 input signals and resolved in frame `world`, `frame_a`, `frame_b` (default)

# Connectors:
- `frame_a`
- `frame_b`
- `torque`: Of type `Blocks.RealInput(3)`. x-, y-, z-coordinates of torque resolved in frame defined by `resolve_frame`.

# Keyword arguments:
- `resolve_frame`: The frame in which the cut force and cut torque are resolved. Default is `:frame_b`, options include `:frame_a` and `:world`.
"""
@component function Torque(; name, resolve_frame = :frame_b)
    @named ptf = PartialTwoFrames()
    @unpack frame_a, frame_b = ptf
    @named begin
        torque = Blocks.RealInput(; nin = 3)
        basicTorque = BasicTorque(; resolve_frame = resolve_frame)
    end

    eqs = [connect(basicTorque.frame_a, frame_a)
           connect(basicTorque.frame_b, frame_b)
           connect(basicTorque.torque, torque)]
    extend(System(eqs, t, name = name, systems = [torque, basicTorque]), ptf)
end

@component function BasicWorldTorque(; name, resolve_frame = :world)
    @named torque = Blocks.RealInput(; nin = 3)
    @named frame_b = Frame()
    eq1 = if resolve_frame == :world
        frame_b.tau ~ -resolve2(ori(frame_b), torque.u)
    elseif resolve_frame == :frame_b
        frame_b.tau ~ -torque.u
    else
        frame_b.tau ~ zeros(3)
    end
    eqs = [eq1; frame_b.f ~ zeros(3)]
    System(eqs, t; name, systems = [torque, frame_b])
end

"""
    WorldTorque(; name, resolve_frame = :world)

External torque acting at `frame_b`, defined by 3 input signals and resolved in frame `:world` or `:frame_b`.

# Connectors:
- `frame_b`: Frame at which the torque is acting
- `torque`: Of type `Blocks.RealInput(3)`. x-, y-, z-coordinates of torque resolved in frame defined by `resolve_frame`.

# Rendering options
- `scale = 0.1`: scaling factor for the force [m/N]
- `color = [0,1,0,0.5]`: color of the force arrow in rendering
- `radius = 0.05`: radius of the force arrow in rendering
"""
@component function WorldTorque(; name, resolve_frame = :world, scale = 0.1, color = [0, 1, 0, 0.5], radius = 0.05)
    @named begin
        torque = Blocks.RealInput(; nin = 3)
        frame_b = Frame()
        basicWorldTorque = BasicWorldTorque(; resolve_frame)
    end
    pars = @parameters begin
        scale = scale, [description = "scaling factor for the force [m/N]"]
        color[1:4] = color, [description = "color of the force arrow in rendering"]
        radius = radius, [description = "radius of the force arrow in rendering"]
    end
    eqs = [
        connect(basicWorldTorque.frame_b, frame_b)
        connect(basicWorldTorque.torque, torque)
    ]
    System(eqs, t, [], pars; name, systems = [torque, basicWorldTorque, frame_b])
end


@component function BasicForce(; name, resolve_frame = :frame_b)
    @named ptf = PartialTwoFrames()
    @named force = Blocks.RealInput(; nin = 3)
    @unpack frame_a, frame_b = ptf
    @variables begin
        (r_0(t)[1:3]),
        [
            description = "Position vector from origin of frame_a to origin of frame_b resolved in world frame",
        ]
        (f_b_0(t)[1:3]), [
            description = "frame_b.f resolved in world frame"]
    end

    eqs = [r_0 ~ frame_b.r_0 - frame_a.r_0
           frame_b.tau ~ 0
           0 ~ frame_a.f + resolve2(frame_a, f_b_0)
           0 ~ frame_a.tau + resolve2(frame_a, cross(r_0, f_b_0))]

    if resolve_frame == :frame_a
        append!(eqs,
                [f_b_0 ~ -resolve1(frame_a, force.u)
                 frame_b.f ~ resolve2(frame_b, f_b_0)])
    elseif resolve_frame == :frame_b
        append!(eqs,
                [f_b_0 ~ -resolve1(frame_b, force.u)
                 frame_b.f ~ -force.u])
    elseif resolve_frame == :world
        append!(eqs,
                [f_b_0 ~ -force.u
                 frame_b.f ~ resolve2(frame_b, f_b_0)])
    else
        error("Unknown value of argument resolve_frame")
    end

    extend(System(eqs, t, name = name, systems = [force]), ptf)
end

@component function BasicWorldForce(; name, resolve_frame = :world)
    @named force = Blocks.RealInput(; nin = 3)
    @named frame_b = Frame()
    eq1 = if resolve_frame == :world
        frame_b.f ~ -resolve2(ori(frame_b), force.u)
    elseif resolve_frame == :frame_b
        frame_b.f ~ -force.u
    else
        frame_b.f ~ zeros(3)
    end
    eqs = [eq1; frame_b.tau ~ zeros(3)]
    System(eqs, t; name, systems = [force, frame_b])
end

"""
    Force(; name, resolve_frame = :frame_b)

Force acting between two frames, defined by 3 input signals and resolved in frame `world`, `frame_a`, `frame_b` (default)

# Connectors:
- `frame_a`
- `frame_b`
- `force`: Of type `Blocks.RealInput(3)`. x-, y-, z-coordinates of force resolved in frame defined by `resolve_frame`.

# Keyword arguments:
- `resolve_frame`: The frame in which the cut force and cut torque are resolved. Default is `:frame_b`, options include `:frame_a` and `:world`.
"""
@component function Force(; name, resolve_frame = :frame_b)
    @named ptf = PartialTwoFrames()
    @unpack frame_a, frame_b = ptf
    @named begin
        force = Blocks.RealInput(; nin = 3) # x-, y-, z-coordinates of force resolved in frame defined by resolve_frame
        basicForce = BasicForce(; resolve_frame = resolve_frame)
    end

    eqs = [connect(basicForce.frame_a, frame_a)
           connect(basicForce.frame_b, frame_b)
           connect(basicForce.force, force)]
    extend(System(eqs, t, name = name, systems = [force, basicForce]), ptf)
end

"""
    WorldForce(; name, resolve_frame = :world)

External force acting at `frame_b`, defined by 3 input signals and resolved in frame `:world` or `:frame_b`.

# Connectors:
- `frame_b`: Frame at which the force is acting
- `force`: Of type `Blocks.RealInput(3)`. x-, y-, z-coordinates of force resolved in frame defined by `resolve_frame`.

# Rendering options
- `scale = 0.1`: scaling factor for the force [m/N]
- `color = [0,1,0,0.5]`: color of the force arrow in rendering
- `radius = 0.05`: radius of the force arrow in rendering
"""
@component function WorldForce(; name, resolve_frame = :world, scale = 0.1, color = [0, 1, 0, 0.5], radius = 0.05)
    @named begin
        frame_b = Frame()
        force = Blocks.RealInput(; nin = 3) # x-, y-, z-coordinates of force resolved in frame defined by resolve_frame
        basicWorldForce = BasicWorldForce(; resolve_frame)
    end
    pars = @parameters begin
        scale = scale, [description = "scaling factor for the force [m/N]"]
        color[1:4] = color, [description = "color of the force arrow in rendering"]
        radius = radius, [description = "radius of the force arrow in rendering"]
    end

    eqs = [
        connect(basicWorldForce.frame_b, frame_b)
        connect(basicWorldForce.force, force)
    ]
    System(eqs, t, [], pars; name, systems = [force, basicWorldForce, frame_b])
end

@component function LineForceBase(; name, length = 0, s_small = 1e-10, fixed_rotation_at_frame_a = false,
    fixed_rotation_at_frame_b = false, r_rel_0 = nothing, s0 = nothing)
    @named frame_a = Frame(varw = fixed_rotation_at_frame_a)
    @named frame_b = Frame(varw = fixed_rotation_at_frame_b)

    vars = @variables begin
        length(t), [
            description = "Distance between the origin of frame_a and the origin of frame_b",
        ]
        s(t)=s0, [
            description = "(Guarded) distance between the origin of frame_a and the origin of frame_b (>= s_small))",
        ]
        r_rel_0(t)[1:3]=r_rel_0, [
            description = "Position vector from frame_a to frame_b resolved in world frame",
        ]
        e_rel_0(t)[1:3], [
            description = "Unit vector in direction from frame_a to frame_b, resolved in world frame",
        ]
    end

    eqs = [r_rel_0 ~ frame_b.r_0 - frame_a.r_0
           length ~ _norm(r_rel_0)
           s ~ max(length, s_small)
           e_rel_0 ~ r_rel_0 / s]

    # Modelica stdlib has the option to inser special equations when two line forces are connected, this option does not yet exisst here https://github.com/modelica/ModelicaStandardLibrary/blob/10238e9927e2078571e41b53cda128c5207f69f7/Modelica/Mechanics/MultiBody/Interfaces/LineForceBase.mo#L49

    if fixed_rotation_at_frame_a
        # TODO: frame_a.R should be rooted here
        eqs = [eqs; vec(ori(frame_a).R .~ nullrotation().R)]
    else
        eqs = [eqs; frame_a.tau .~ 0]
    end

    if fixed_rotation_at_frame_b
        # TODO: frame_b.R should be rooted here
        eqs = [eqs; vec(ori(frame_b).R .~ nullrotation().R)]
    else
        eqs = [eqs; frame_b.tau .~ 0]
    end

    compose(System(eqs, t; name), frame_a, frame_b)
end

@component function LineForceWithMass(; name, length = 0, m = 1.0, lengthfraction = 0.5, kwargs...)
    m0 = m
    @named lfb = LineForceBase(; length, kwargs...)
    @unpack length, s, r_rel_0, e_rel_0, frame_a, frame_b = lfb
    @named flange_a = TP.Flange()
    @named flange_b = TP.Flange()
    @parameters m=m [description = "mass", bounds = (0, Inf)]
    @variables fa(t) [description = "scalar force from flange_a"]
    @variables fb(t) [description = "scalar force from flange_b"]
    @variables r_CM_0(t)[1:3] [
        description = "Position vector from world frame to point mass, resolved in world frame",
    ]
    @variables v_CM_0(t)[1:3] [description = "First derivative of r_CM_0"]
    @variables ag_CM_0(t)[1:3] [description = "D(v_CM_0) - gravityAcceleration"]
    @parameters lengthfraction=lengthfraction [
        description = "Location of point mass with respect to frame_a as a fraction of the distance from frame_a to frame_b",
        bounds = (0, 1),
    ]

    eqs = [flange_a.s ~ 0
           flange_b.s ~ length]

    # Determine translational flange forces
    eqs = [eqs
           fa ~ flange_a.f
           fb ~ flange_b.f]

    # NOTE, both frames are assumed to be connected, while modelica has special handling if they aren't
    if m0 > 0
        eqs = [eqs
               r_CM_0 ~ frame_a.r_0 + r_rel_0 * lengthfraction
               v_CM_0 ~ D(r_CM_0)
               ag_CM_0 ~ D(v_CM_0) - gravity_acceleration(r_CM_0)
               frame_a.f ~ resolve2(ori(frame_a),
                                     (m * (1 - lengthfraction)) * ag_CM_0 - e_rel_0 * fa)
               frame_b.f ~ resolve2(ori(frame_b),
                                     (m * lengthfraction) * ag_CM_0 - e_rel_0 * fb)]
    else
        eqs = [eqs
               r_CM_0 ~ zeros(3)
               v_CM_0 ~ zeros(3)
               ag_CM_0 ~ zeros(3)
               frame_a.f ~ -resolve2(ori(frame_a), fa*e_rel_0)
               frame_b.f ~ -resolve2(ori(frame_b), fb*e_rel_0)]
    end

    extend(System(eqs, t; name, systems = [flange_a, flange_b]), lfb)
end

@component function PartialLineForce(; name, kwargs...)
    @named lfb = LineForceBase(; kwargs...)
    @unpack length, s, r_rel_0, e_rel_0, frame_a, frame_b = lfb

    @variables begin
        (r_rel_a(t)[1:3]),
        [
            description = "Position vector from origin of frame_a to origin of frame_b, resolved in frame_a",
        ]
        (e_a(t)[1:3]),
        [
            description = "Unit vector on the line connecting the origin of frame_a with the origin of frame_b resolved in frame_a (directed from frame_a to frame_b)",
        ]
        (f(t)),
        [
            description = "Line force acting on frame_a and on frame_b (positive, if acting on frame_b and directed from frame_a to frame_b)",
        ]
    end
    equations = Equation[
                 # Determine relative position vector between the two frames
                 r_rel_a ~ resolve2(frame_a, r_rel_0)
                 e_a ~ r_rel_a / s

                 # Determine forces and torques at frame_a and frame_b
                 frame_a.f ~ -e_a * f
                 frame_b.f ~ -resolve2(relative_rotation(frame_a, frame_b),
                                                 frame_a.f)]
    extend(System(equations, t; name), lfb)
end

"""
    Spring(; c, name, m = 0, lengthfraction = 0.5, s_unstretched = 0, kwargs)

Linear spring acting as line force between `frame_a` and `frame_b`.
A force `f` is exerted on the origin of `frame_b` and with opposite sign
on the origin of `frame_a` along the line from the origin of `frame_a` to the origin
of `frame_b` according to the equation:
```math
f = c s
```
where `c` is the spring stiffness parameter, `s` is the
distance between the origin of `frame_a` and the origin of `frame_b`.

Optionally, the mass of the spring is taken into account by a
point mass located on the line between `frame_a` and `frame_b`
(default: middle of the line). If the spring mass is zero, the
additional equations to handle the mass are removed.

# Arguments:
- `c`: Spring stiffness
- `m`: Mass of the spring (can be zero)
- `lengthfraction`: Location of spring mass with respect to `frame_a` as a fraction of the distance from `frame_a` to `frame_b` (=0: at `frame_a`; =1: at `frame_b`)
- `s_unstretched`: Length of the spring when it is unstretched
- `kwargs`: are passed to `LineForceWithMass`

# Rendering
- `num_windings = 6`: Number of windings of the coil when rendered
- `color = [0,0,1,1]`: Color of the spring when rendered
- `radius = 0.1`: Radius of spring when rendered
- `N = 200`: Number of points in mesh when rendered. Rendering time can be reduced somewhat by reducing this number.
- `end_ratio = 0.0`: Ratio of the length of the spring [0, 0.5] that is rendered with decreasing radius at the ends. Set this to 0 to have uniform radius along the entire spring.

See also [`SpringDamperParallel`](@ref)
"""
@component function Spring(; c, name, m = 0, lengthfraction = 0.5, s_unstretched = 0, num_windings=6, color=[0,0,1,1], radius=0.1, N=200, end_ratio = 0.0, kwargs...)
    @named ptf = PartialTwoFrames()
    @unpack frame_a, frame_b = ptf
    pars = @parameters begin
        c=c, [description = "spring constant", bounds = (0, Inf)]
        s_unstretched=s_unstretched, [
            description = "unstretched length of spring",
            bounds = (0, Inf),
        ]
        num_windings = num_windings, [description = "Number of windings of the coil when rendered"]
        end_ratio = 0.0
        color[1:4] = color
        radius = radius, [description = "Radius of spring when rendered"]
        N = N, [description = "Number of points in mesh when rendered"]
    end
    @named lineforce = LineForceWithMass(; length = s_unstretched, m, lengthfraction,
                                         kwargs...)
    
    @named spring2d = TP.Spring(; c, s_rel0 = s_unstretched)

    @variables r_rel_a(t)[1:3] [
        description = "Position vector from origin of frame_a to origin of frame_b, resolved in frame_a",
    ]
    @variables e_a(t)[1:3] [
        description = "Unit vector on the line connecting the origin of frame_a with the origin of frame_b resolved in frame_a (directed from frame_a to frame_b)",
    ]
    @variables f(t) [
        description = "Line force acting on frame_a and on frame_b (positive, if acting on frame_b and directed from frame_a to frame_b)",
    ]
    @variables length(t) [
        description = "Distance between the origin of frame_a and the origin of frame_b",
    ]
    @variables s(t) [
        description = "(Guarded) distance between the origin of frame_a and the origin of frame_b (>= s_small))",
    ]
    @variables v(t) [
        description = "derivative of s",
    ]
    @variables r_rel_0(t)[1:3] [
        description = "Position vector from frame_a to frame_b resolved in world frame",
    ]
    @variables e_rel_0(t)[1:3] [
        description = "Unit vector in direction from frame_a to frame_b, resolved in world frame",
    ]

    eqs = [D(s) ~ v
           r_rel_a ~ resolve2(ori(frame_a), r_rel_0)
           e_a ~ r_rel_a / s
           f ~ spring2d.f
           length ~ lineforce.length
           s ~ lineforce.s
           r_rel_0 ~ lineforce.r_rel_0
           e_rel_0 ~ lineforce.e_rel_0
           connect(lineforce.frame_a, frame_a)
           connect(lineforce.frame_b, frame_b)
           connect(spring2d.flange_b, lineforce.flange_b)
           connect(spring2d.flange_a, lineforce.flange_a)]

    sys = extend(System(eqs, t; name=:nothing, systems = [lineforce, spring2d]), ptf)
    add_params(sys, pars; name)
end

"""
    Damper(; d, name, kwargs)

Linear damper acting as line force between `frame_a` and `frame_b`.
A force `f` is exerted on the origin of `frame_b` and with opposite sign
on the origin of `frame_a` along the line from the origin of `frame_a` to the origin
of `frame_b` according to the equation:
```math
f = d D(s)
```
where `d` is the (viscous) damping parameter, `s` is the
distance between the origin of `frame_a` and the origin of `frame_b`
and `D(s)` is the time derivative of `s`.

# Arguments:
- `d`: Damping coefficient

# Rendering
- `radius = 0.1`: Radius of damper when rendered
- `length_fraction = 0.2`: Fraction of the length of the damper that is rendered
- `color = [0.5, 0.5, 0.5, 1]`: Color of the damper when rendered

See also [`SpringDamperParallel`](@ref)
"""
@component function Damper(; d, name, radius = 0.1, length_fraction = 0.2, color = [0.5, 0.5, 0.5, 1],kwargs...)
    @named plf = PartialLineForce(; kwargs...)
    @unpack s, f = plf
    pars = @parameters begin
        d=d, [description = "damping constant", bounds = (0, Inf)]
        radius = radius, [description = "Radius of damper when rendered"]
        length_fraction = length_fraction, [description = "Fraction of the length of the damper that is rendered"]
        color[1:4] = color
    end
    eqs = [
        f ~ d * D(s),
    ]
    extend(System(eqs, t, [s, f], pars; name), plf)
end

"""
    SpringDamperParallel(; name, c, d, s_unstretched)

Linear spring and linear damper in parallel acting as line force between `frame_a` and `frame_b`. A force `f` is exerted on the origin of `frame_b` and with opposite sign on the origin of `frame_a` along the line from the origin of `frame_a` to the origin of `frame_b` according to the equation:
```math
f = c (s - s_{unstretched}) + d \\cdot D(s)
```
where `c`, `s_unstretched` and `d` are parameters, `s` is the distance between the origin of `frame_a` and the origin of `frame_b` and `D(s)` is the time derivative of `s`.
"""
@component function SpringDamperParallel(; name, c, d, s_unstretched=0,
    color = [0, 0, 1, 1], radius = 0.1, N = 200, num_windings = 6, end_ratio = 0.0, kwargs...)
    @named plf = PartialLineForce(; kwargs...)
    @unpack s, f = plf

    pars = @parameters begin
        c=c, [description = "spring constant", bounds = (0, Inf)]
        d=d, [description = "damping constant", bounds = (0, Inf)]
        s_unstretched=s_unstretched, [
            description = "unstretched length of spring",
            bounds = (0, Inf),
        ]
        color[1:4] = color, [description = "Color of the spring when rendered"]
        radius = radius, [description = "Radius of spring when rendered"]
        N = N, [description = "Number of points in mesh when rendered"]
        num_windings = num_windings, [description = "Number of windings of the coil when rendered"]
        end_ratio = end_ratio
    end
    # pars = collect_all(pars)
    

    f_d = d * D(s)
    eqs = [
           f ~ c * (s - s_unstretched) + f_d
           # lossPower ~ f_d*der(s)
           ]
    extend(System(eqs, t, [], pars; name), plf)
end
