function PartialRelativeBaseSensor(; name)
    @named begin
        frame_a = Frame()
        frame_b = Frame()
    end

    equations = [frame_a.f .~ zeros(3) |> collect
                 frame_a.tau .~ zeros(3) |> collect
                 frame_b.f .~ zeros(3) |> collect
                 frame_b.tau .~ zeros(3) |> collect
                 frame_resolve.f .~ zeros(3) |> collect
                 frame_resolve.tau .~ zeros(3) |> collect]
    compose(ODESystem(equations, t; name), frame_a, frame_b)
end

function PartialAbsoluteSensor(; name, n_out)
    @named begin
        frame_a = Frame()
        y = Blocks.RealOutput(nout = n_out)
    end
    equations = []
    compose(ODESystem(equations, t; name), frame_a, y)
end

"""
    PartialCutForceBaseSensor(; name, resolveInFrame = :frame_a)

- `resolveInFrame`: The frame in which the cut force and cut torque are resolved. Default is `:frame_a`, options include `:frame_a` and `:world`.
"""
function PartialCutForceBaseSensor(; name, resolveInFrame = :frame_a)
    @named begin
        frame_a = Frame()
        frame_b = Frame()
    end

    equations = [frame_a.r_0 .~ frame_b.r_0 |> collect
                 ori(frame_a) ~ ori(frame_b)
                 zeros(3) .~ frame_a.f + frame_b.f |> collect
                 zeros(3) .~ frame_a.tau + frame_b.tau |> collect]
    compose(ODESystem(equations, t; name), frame_a, frame_b)
end

"""
    CutTorque(; name, resolveInFrame)

Basic sensor to measure cut torque vector. Contains a connector of type `Blocks.RealOutput` with name `torque`.

- `resolveInFrame`: The frame in which the cut force and cut torque are resolved. Default is `:frame_a`, options include `:frame_a` and `:world`.
"""
function CutTorque(; name, resolveInFrame = :frame_a)
    @named pcfbs = PartialCutForceBaseSensor(; resolveInFrame)
    @named torque = Blocks.RealOutput(nout = 3) # "Cut torque resolved in frame defined by resolveInFrame"
    @unpack frame_a, frame_b = pcfbs
    eqs = if resolveInFrame === :world
        collect(torque.u) .~ resolve1(ori(frame_a), frame_a.tau)
    elseif resolveInFrame === :frame_a
        collect(torque.u) .~ collect(frame_a.tau)
    else
        error("resolveInFrame must be :world or :frame_a")
    end
    extend(compose(ODESystem(eqs, t; name), torque), pcfbs)
end

"""
    BasicCutForce(; name, resolveInFrame)

Basic sensor to measure cut force vector. Contains a connector of type `Blocks.RealOutput` with name `force`.

- `resolveInFrame`: The frame in which the cut force and cut torque are resolved. Default is `:frame_a`, options include `:frame_a` and `:world`.
"""
function CutForce(; name, resolveInFrame = :frame_a)
    @named pcfbs = PartialCutForceBaseSensor(; resolveInFrame)
    @named force = Blocks.RealOutput(nout = 3) # "Cut force resolved in frame defined by resolveInFrame"
    @unpack frame_a, frame_b = pcfbs
    eqs = if resolveInFrame === :world
        collect(force.u) .~ resolve1(ori(frame_a), frame_a.f)
    elseif resolveInFrame === :frame_a
        collect(force.u) .~ collect(frame_a.f)
    else
        error("resolveInFrame must be :world or :frame_a")
    end
    extend(compose(ODESystem(eqs, t; name), force), pcfbs)
end

function RelativePosition(; name, resolveInFrame = :frame_a)
    @named begin prs = PartialRelativeBaseSensor(; name) end

    @unpack frame_a, frame_b = prs

    equations = [frame_a.r_0 .~ frame_b.r_0 |> collect
                 ori(frame_a) ~ ori(frame_b)
                 zeros(3) .~ frame_a.r_0 - frame_b.r_0 |> collect]
    extend(compose(ODESystem(equations, t; name), frame_a, frame_b), prs)
end

function RelativeAngles(; name, sequence = [1, 2, 3])
    @named begin
        frame_a = Frame()
        frame_b = Frame()
        angles = Blocks.RealOutput(nout = 3)
    end
    @named R_rel = NumRotationMatrix()
    eqs = [frame_a.f .~ zeros(3) |> collect
           frame_a.tau .~ zeros(3) |> collect
           frame_b.f .~ zeros(3) |> collect
           frame_b.tau .~ zeros(3) |> collect
           R_rel ~ relativeRotation(frame_a, frame_b)
           angles .~ axesRotationsAngles(R_rel, sequence, guessAngle1)]
    compose(ODESystem(eqs, t; name), frame_a, frame_b, angles)
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
           angles.u .~ axesRotationsAngles(ori(frame_a), [1, 2, 3])]
    extend(compose(ODESystem(eqs, t; name)), pas)
end
