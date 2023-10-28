import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica as TP

import ModelingToolkitStandardLibrary.Blocks

"""
    BasicTorque(; name, resolveInFrame = :world)

Low-level torque component used to build [`Torque`](@ref)
"""
function BasicTorque(; name, resolveInFrame = :world)
    @named ptf = PartialTwoFrames()
    @named torque = Blocks.RealInput(; nin = 3)
    @unpack frame_a, frame_b = ptf
    @variables begin
        (r_0(t)[1:3] = zeros(3)),
        [
            description = "Position vector from origin of frame_a to origin of frame_b resolved in world frame",
        ]
        (t_b_0(t)[1:3] = zeros(3)), [
            description = "frame_b.tau resolved in world frame"]
    end
    r_0, t_b_0 = collect.((r_0, t_b_0))

    eqs = [collect(r_0 .~ frame_b.r_0 - frame_a.r_0)
           collect(frame_a.f .~ zeros(3))
           collect(frame_b.f .~ zeros(3))
           # torque balance
           zeros(3) .~ collect(frame_a.tau) + resolve2(frame_a, t_b_0)]

    if resolveInFrame == :frame_a
        append!(eqs,
                [t_b_0 .~ -resolve1(frame_a, torque.u)
                 collect(frame_b.tau) .~ resolve2(frame_b, t_b_0)])
    elseif resolveInFrame == :frame_b
        append!(eqs,
                [t_b_0 .~ -resolve1(frame_b, torque.u)
                 collect(frame_b.tau) .~ collect(-torque.u)])
    elseif resolveInFrame == :world
        append!(eqs,
                [t_b_0 .~ collect(-torque.u)
                 collect(frame_b.tau) .~ resolve2(frame_b, t_b_0)])
    else
        error("Unknown value of argument resolveInFrame")
        append!(eqs, [t_b_0 .~ zeros(3)
                      collect(frame_b.tau) .~ zeros(3)])
    end

    extend(ODESystem(eqs, t, name = name, systems = [torque]), ptf)
end

"""
    Torque(; name, resolveInFrame = :frame_b)

Torque acting between two frames, defined by 3 input signals and resolved in frame `world`, `frame_a`, `frame_b` (default)

# Connectors:
- `frame_a`
- `frame_b`
- `torque`: Of type `Blocks.RealInput(3)`. x-, y-, z-coordinates of torque resolved in frame defined by `resolveInFrame`.

# Keyword arguments:
- `resolveInFrame`: The frame in which the cut force and cut torque are resolved. Default is `:frame_b`, options include `:frame_a` and `:world`.
"""
function Torque(; name, resolveInFrame = :frame_b)
    @named ptf = PartialTwoFrames()
    @unpack frame_a, frame_b = ptf
    @named begin
        torque = Blocks.RealInput(; nin = 3)
        basicTorque = BasicTorque(; resolveInFrame = resolveInFrame)
    end

    eqs = [connect(basicTorque.frame_a, frame_a)
           connect(basicTorque.frame_b, frame_b)
           connect(basicTorque.torque, torque)]
    extend(ODESystem(eqs, t, name = name, systems = [torque, basicTorque]), ptf)
end

function BasicForce(; name, resolveInFrame = :frame_b)
    @named ptf = PartialTwoFrames()
    @named force = Blocks.RealInput(; nin = 3)
    @unpack frame_a, frame_b = ptf
    @variables begin
        (r_0(t)[1:3] = zeros(3)),
        [
            description = "Position vector from origin of frame_a to origin of frame_b resolved in world frame",
        ]
        (f_b_0(t)[1:3] = zeros(3)), [
            description = "frame_b.f resolved in world frame"]
    end
    r_0, f_b_0 = collect.((r_0, f_b_0))

    eqs = [collect(r_0 .~ frame_b.r_0 - frame_a.r_0)
           0 .~ collect(frame_a.f) + resolve2(frame_a, f_b_0)
           0 .~ collect(frame_a.tau) + resolve2(frame_a, cross(r_0, f_b_0))]

    if resolveInFrame == :frame_a
        append!(eqs,
                [f_b_0 .~ -resolve1(frame_a, force.u)
                 collect(frame_b.tau) .~ resolve2(frame_b, f_b_0)])
    elseif resolveInFrame == :frame_b
        append!(eqs,
                [f_b_0 .~ -resolve1(frame_b, force.u)
                 collect(frame_b.tau) .~ collect(-force.u)])
    elseif resolveInFrame == :world
        append!(eqs,
                [f_b_0 .~ collect(-force.u)
                 collect(frame_b.tau) .~ resolve2(frame_b, f_b_0)])
    else
        error("Unknown value of argument resolveInFrame")
    end

    extend(ODESystem(eqs, t, name = name, systems = [force]), ptf)
end

"""
    Force(; name, resolveInFrame = :frame_b)

Force acting between two frames, defined by 3 input signals and resolved in frame `world`, `frame_a`, `frame_b` (default)

# Connectors:
- `frame_a`
- `frame_b`
- `force`: Of type `Blocks.RealInput(3)`. x-, y-, z-coordinates of force resolved in frame defined by `resolveInFrame`.

# Keyword arguments:
- `resolveInFrame`: The frame in which the cut force and cut torque are resolved. Default is `:frame_b`, options include `:frame_a` and `:world`.
"""
function Force(; name, resolveInFrame = :frame_b)
    @named ptf = PartialTwoFrames()
    @unpack frame_a, frame_b = ptf
    @named begin
        force = Blocks.RealInput(; nin = 3) # x-, y-, z-coordinates of force resolved in frame defined by resolveInFrame
        basicForce = BasicForce(; resolveInFrame = resolveInFrame)
    end

    eqs = [connect(basicForce.frame_a, frame_a)
           connect(basicForce.frame_b, frame_b)
           connect(basicForce.force, force)]
    extend(ODESystem(eqs, t, name = name, systems = [force, basicForce]), ptf)
end

function LineForceBase(; name, length = 0, s_small = 1e-10, fixedRotationAtFrame_a = false,
                       fixedRotationAtFrame_b = false, r_rel_0 = 0, s0 = 0)
    @named frame_a = Frame()
    @named frame_b = Frame()

    @variables length(t) [
        description = "Distance between the origin of frame_a and the origin of frame_b",
    ]
    @variables s(t)=s0 [
        description = "(Guarded) distance between the origin of frame_a and the origin of frame_b (>= s_small))",
    ]
    @variables r_rel_0(t)[1:3]=r_rel_0 [
        description = "Position vector from frame_a to frame_b resolved in world frame",
    ]
    @variables e_rel_0(t)[1:3] [
        description = "Unit vector in direction from frame_a to frame_b, resolved in world frame",
    ]

    eqs = [r_rel_0 .~ frame_b.r_0 .- frame_a.r_0;
           length ~ norm(r_rel_0)
           s ~ max(length, s_small)
           e_rel_0 .~ r_rel_0 ./ s]

    # Modelica stdlib has the option to inser special equations when two line forces are connected, this option does not yet exisst here https://github.com/modelica/ModelicaStandardLibrary/blob/10238e9927e2078571e41b53cda128c5207f69f7/Modelica/Mechanics/MultiBody/Interfaces/LineForceBase.mo#L49

    if fixedRotationAtFrame_a
        # TODO: frame_a.R should be rooted here
        eqs = [eqs; vec(ori(frame_a).R .~ nullRotation().R)]
    else
        eqs = [eqs; frame_a.tau .~ 0]
    end

    if fixedRotationAtFrame_b
        # TODO: frame_b.R should be rooted here
        eqs = [eqs; vec(ori(frame_b).R .~ nullRotation().R)]
    else
        eqs = [eqs; frame_b.tau .~ 0]
    end

    compose(ODESystem(eqs, t; name), frame_a, frame_b)
end

function LineForceWithMass(; name, length = 0, m = 1.0, lengthFraction = 0.5, kwargs...)
    m0 = m
    @named lfb = LineForceBase(; length, kwargs...)
    @unpack length, s, r_rel_0, e_rel_0, frame_a, frame_b = lfb
    @named flange_a = TP.Flange()
    @named flange_b = TP.Flange()
    @parameters m=m [description = "mass", bounds = (0, Inf)]
    @variables fa(t)=0 [description = "scalar force from flange_a"]
    @variables fb(t)=0 [description = "scalar force from flange_b"]
    @variables r_CM_0(t)[1:3]=zeros(3) [
        description = "Position vector from world frame to point mass, resolved in world frame",
    ]
    @variables v_CM_0(t)[1:3]=zeros(3) [description = "First derivative of r_CM_0"]
    @variables ag_CM_0(t)[1:3]=zeros(3) [description = "D(v_CM_0) - gravityAcceleration"]
    @parameters lengthFraction=lengthFraction [
        description = "Location of point mass with respect to frame_a as a fraction of the distance from frame_a to frame_b",
        bounds = (0, 1),
    ]

    r_rel_0, e_rel_0, r_CM_0, v_CM_0, ag_CM_0 = collect.((r_rel_0, e_rel_0, r_CM_0, v_CM_0,
                                                          ag_CM_0))

    eqs = [flange_a.s ~ 0
           flange_b.s ~ length]

    # Determine translational flange forces
    eqs = [eqs
           fa ~ flange_a.f
           fb ~ flange_b.f]

    # NOTE, both frames are assumed to be connected, while modelica has special handling if they aren't
    if m0 > 0
        eqs = [eqs
               r_CM_0 .~ frame_a.r_0 + r_rel_0 * lengthFraction
               v_CM_0 .~ D.(r_CM_0)
               ag_CM_0 .~ D.(v_CM_0) - gravity_acceleration(r_CM_0)
               frame_a.f .~ resolve2(ori(frame_a),
                                     (m * (1 - lengthFraction)) * ag_CM_0 - e_rel_0 * fa)
               frame_b.f .~ resolve2(ori(frame_b),
                                     (m * lengthFraction) * ag_CM_0 - e_rel_0 * fb)]
    else
        eqs = [eqs
               r_CM_0 .~ zeros(3)
               v_CM_0 .~ zeros(3)
               ag_CM_0 .~ zeros(3)
               frame_a.f .~ -resolve2(ori(frame_a), e_rel_0 * fa)
               frame_b.f .~ -resolve2(ori(frame_b), e_rel_0 * fb)]
    end

    extend(ODESystem(eqs, t; name, systems = [flange_a, flange_b]), lfb)
end

function PartialLineForce(; name, kwargs...)
    @named lfb = LineForceBase(; kwargs...)
    @unpack length, s, r_rel_0, e_rel_0, frame_a, frame_b = lfb

    @variables begin
        (r_rel_a(t)[1:3] = 0),
        [
            description = "Position vector from origin of frame_a to origin of frame_b, resolved in frame_a",
        ]
        (e_a(t)[1:3] = 0),
        [
            description = "Unit vector on the line connecting the origin of frame_a with the origin of frame_b resolved in frame_a (directed from frame_a to frame_b)",
        ]
        (f(t) = 0),
        [
            description = "Line force acting on frame_a and on frame_b (positive, if acting on frame_b and directed from frame_a to frame_b)",
        ]
    end
    equations = [
                 # Determine relative position vector between the two frames
                 collect(r_rel_a) .~ resolve2(frame_a, r_rel_0)
                 collect(e_a) .~ collect(r_rel_a ./ s)

                 # Determine forces and torques at frame_a and frame_b
                 collect(frame_a.f) .~ collect(-e_a * f)
                 collect(frame_b.f) .~ -resolve2(relativeRotation(frame_a, frame_b),
                                                 frame_a.f)]
    extend(ODESystem(equations, t; name), lfb)
end

"""
    Spring(; c, name, m = 0, lengthFraction = 0.5, s_unstretched = 0, kwargs)

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
- `lengthFraction`: Location of spring mass with respect to `frame_a` as a fraction of the distance from `frame_a` to `frame_b` (=0: at `frame_a`; =1: at `frame_b`)
- `s_unstretched`: Length of the spring when it is unstretched
- `kwargs`: are passed to `LineForceWithMass`

See also [`SpringDamperParallel`](@ref)
"""
@component function Spring(; c, name, m = 0, lengthFraction = 0.5, s_unstretched = 0, kwargs...)
    @named ptf = PartialTwoFrames()
    @unpack frame_a, frame_b = ptf
    @named lineForce = LineForceWithMass(; length = s_unstretched, m, lengthFraction,
                                         kwargs...)
    @named spring2d = TP.Spring(; c, s_rel0 = s_unstretched)
    @parameters c=c [description = "spring constant", bounds = (0, Inf)]
    @parameters s_unstretched=s_unstretched [
        description = "unstretched length of spring",
        bounds = (0, Inf),
    ]

    @variables r_rel_a(t)[1:3]=0 [
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
           r_rel_a .~ resolve2(ori(frame_a), r_rel_0)
           e_a .~ r_rel_a / s
           f ~ spring2d.f
           length ~ lineForce.length
           s ~ lineForce.s
           r_rel_0 .~ lineForce.r_rel_0
           e_rel_0 .~ lineForce.e_rel_0
           connect(lineForce.frame_a, frame_a)
           connect(lineForce.frame_b, frame_b)
           connect(spring2d.flange_b, lineForce.flange_b)
           connect(spring2d.flange_a, lineForce.flange_a)]

    extend(ODESystem(eqs, t; name, systems = [lineForce, spring2d]), ptf)
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

See also [`SpringDamperParallel`](@ref)
"""
@component function Damper(; d, name, kwargs...)
    @named plf = PartialLineForce(; kwargs...)
    @unpack s, f = plf
    @parameters d=d [description = "damping constant", bounds = (0, Inf)]
    eqs = [
        f ~ d * D(s),
    ]
    extend(ODESystem(eqs, t; name), plf)
end

"""
    SpringDamperParallel(; name, c, d, s_unstretched)

Linear spring and linear damper in parallel acting as line force between `frame_a` and `frame_b`. A force `f` is exerted on the origin of `frame_b` and with opposite sign on the origin of `frame_a` along the line from the origin of `frame_a` to the origin of `frame_b` according to the equation:
```math
f = c (s - s_{unstretched}) + d\cdot D(s)
```
where `c`, `s_unstretched` and `d` are parameters, `s` is the distance between the origin of `frame_a` and the origin of `frame_b` and `D(s)` is the time derivative of `s`.
"""
@component function SpringDamperParallel(; name, c, d, s_unstretched)
    @named plf = PartialLineForce(; kwargs...)
    @unpack s, f = plf

    @parameters c=c [description = "spring constant", bounds = (0, Inf)]
    @parameters d=d [description = "damping constant", bounds = (0, Inf)]
    @parameters s_unstretched=s_unstretched [
        description = "unstretched length of spring",
        bounds = (0, Inf),
    ]

    eqs = [f_d ~ d * D(s)
           f ~ c * (s - s_unstretched) + f_d
           # lossPower ~ f_d*der(s)
           ]
    extend(ODESystem(eqs, t; name), plf)
end
