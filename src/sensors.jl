function PartialRelativeBaseSensor(; name)
    @named begin
        frame_a = Frame()
        frame_b = Frame()
    end

    equations = [frame_a.f .~ zeros(3) |> collect
                 frame_a.tau .~ zeros(3) |> collect
                 frame_b.f .~ zeros(3) |> collect
                 frame_b.tau .~ zeros(3) |> collect]
    compose(System(equations, t; name), frame_a, frame_b)
end

function BasicRelativePosition(; name, resolve_frame)
    @named prb = PartialRelativeBaseSensor()
    @unpack frame_a, frame_b = prb
    systems = @named begin
        r_rel = RealOutput(nout = 3)
    end

    d = collect(frame_b.r_0 - frame_a.r_0)
    eqs = if resolve_frame === :frame_a
        collect(r_rel.u) .~ resolve2(ori(frame_a), d)
    elseif resolve_frame === :frame_b
        collect(r_rel.u) .~ resolve2(ori(frame_b), d)
    elseif resolve_frame === :world
        collect(r_rel.u) .~ d
    else
        error("resolve_frame must be :world, :frame_a or :frame_b, you provided $resolve_frame, which makes me sad.")
    end

    extend(compose(System(eqs, t; name), r_rel), prb)
end

function RelativePosition(; name, resolve_frame = :frame_a)

    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
        relativePosition = BasicRelativePosition(; resolve_frame)
        r_rel = RealOutput(nout = 3)
    end

    eqs = [
        connect(relativePosition.frame_a, frame_a)
        connect(relativePosition.frame_b, frame_b)
        connect(relativePosition.r_rel, r_rel)
    ]
    compose(System(eqs, t; name), frame_a, frame_b, r_rel, relativePosition)
end

function PartialAbsoluteSensor(; name, n_out)
    @named begin
        frame_a = Frame()
        y = Blocks.RealOutput(nout = n_out)
    end
    equations = []
    compose(System(equations, t; name), frame_a, y)
end

"""
    PartialCutForceBaseSensor(; name, resolve_frame = :frame_a)

- `resolve_frame`: The frame in which the cut force and cut torque are resolved. Default is `:frame_a`, options include `:frame_a` and `:world`.
"""
function PartialCutForceBaseSensor(; name, resolve_frame = :frame_a)
    @named begin
        frame_a = Frame()
        frame_b = Frame()
    end

    equations = [frame_a.r_0 .~ frame_b.r_0 |> collect
                 ori(frame_a) ~ ori(frame_b)
                 zeros(3) .~ frame_a.f + frame_b.f |> collect
                 zeros(3) .~ frame_a.tau + frame_b.tau |> collect]
    compose(System(equations, t; name), frame_a, frame_b)
end

"""
    CutTorque(; name, resolve_frame)

Basic sensor to measure cut torque vector. Contains a connector of type `Blocks.RealOutput` with name `torque`.

- `resolve_frame`: The frame in which the cut force and cut torque are resolved. Default is `:frame_a`, options include `:frame_a` and `:world`.
"""
function CutTorque(; name, resolve_frame = :frame_a)
    @named pcfbs = PartialCutForceBaseSensor(; resolve_frame)
    @named torque = Blocks.RealOutput(nout = 3) # "Cut torque resolved in frame defined by resolve_frame"
    @unpack frame_a, frame_b = pcfbs
    eqs = if resolve_frame === :world
        collect(torque.u) .~ resolve1(ori(frame_a), frame_a.tau)
    elseif resolve_frame === :frame_a
        collect(torque.u) .~ collect(frame_a.tau)
    else
        error("resolve_frame must be :world or :frame_a")
    end
    extend(compose(System(eqs, t; name), torque), pcfbs)
end

"""
    BasicCutForce(; name, resolve_frame)

Basic sensor to measure cut force vector. Contains a connector of type `Blocks.RealOutput` with name `force`.

- `resolve_frame`: The frame in which the cut force and cut torque are resolved. Default is `:frame_a`, options include `:frame_a` and `:world`.
"""
function CutForce(; name, resolve_frame = :frame_a)
    @named pcfbs = PartialCutForceBaseSensor(; resolve_frame)
    @named force = Blocks.RealOutput(nout = 3) # "Cut force resolved in frame defined by resolve_frame"
    @unpack frame_a, frame_b = pcfbs
    eqs = if resolve_frame === :world
        collect(force.u) .~ resolve1(ori(frame_a), frame_a.f)
    elseif resolve_frame === :frame_a
        collect(force.u) .~ collect(frame_a.f)
    else
        error("resolve_frame must be :world or :frame_a")
    end
    extend(compose(System(eqs, t; name), force), pcfbs)
end

function RelativeAngles(; name, sequence = [1, 2, 3])
    @named begin
        frame_a = Frame()
        frame_b = Frame()
        angles = Blocks.RealOutput(nout = 3)
    end
    @named Rrel = NumRotationMatrix()
    eqs = [frame_a.f .~ zeros(3) |> collect
           frame_a.tau .~ zeros(3) |> collect
           frame_b.f .~ zeros(3) |> collect
           frame_b.tau .~ zeros(3) |> collect
           Rrel ~ relative_rotation(frame_a, frame_b)
           angles .~ axes_rotationangles(Rrel, sequence, guessAngle1)]
    compose(System(eqs, t; name), frame_a, frame_b, angles)
end

function AbsoluteAngles(; name, sequence = [1, 2, 3])
    @named begin
        pas = PartialAbsoluteSensor(; name)
        angles = Blocks.RealOutput(nout = 3)
    end
    @unpack frame_a = pas
    @named R_abs = NumRotationMatrix()
    eqs = [collect(frame_a.f .~ 0)
           collect(frame_a.tau .~ 0)
           angles.u .~ axes_rotationangles(ori(frame_a), [1, 2, 3])]
    extend(compose(System(eqs, t; name)), pas)
end


"""
    Power(; name)

A sensor measuring mechanical power transmitted from `frame_a` to `frame_b`.

# Connectors:
`power` of type `RealOutput`.
"""
function Power(; name)
    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
        power = RealOutput()
    end
    eqs = [
        frame_a.r_0 .~ frame_b.r_0 |> collect
        ori(frame_a) ~ ori(frame_b)
        0 .~ frame_a.f + frame_b.f |> collect
        0 .~ frame_a.tau + frame_b.tau |> collect
        power.u ~ collect(frame_a.f)'resolve2(frame_a, D.(frame_a.r_0)) +
                    collect(frame_a.tau)'angular_velocity2(ori(frame_a))
    ]
    System(eqs, t; name, systems)
end