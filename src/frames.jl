function Oriantation(; name)
    # T: Transformation matrix from world frame to local frame
    # w: Absolute angular velocity of local frame, resolved in local frame

    @variables T(t)[1:3, 1:3] w(t)[1:3]

    ODESystem(Equation[], t, [vec(T); w], [], name = name)
end

@connector function Frame(; name)
    # r_0: Position vector from world frame to the connector frame origin, resolved in world frame
    # R: Orientation object to rotate the world frame into the connector frame
    # f: Cut-force resolved in connector frame
    # tau: Cut-torque resolved in connector frame
    @variables r_0(t)[1:3] f(t)[1:3] [connect = Flow] tau(t)[1:4] [connect = Flow]
    @named R = Oriantation()

    compose(ODESystem(Equation[], t, [r_0; vec(R); w; f; tau], [], name = name), R)
end

# TODO: set angular velocity to zero as well
# Should `~` be defined for `System ~ System` and `System ~ NamedTuple`?
nullrotation() = zeros(Int, 3, 3)

"""
    h2 = resolve2(R12, h1)

`R` is a 3x3 matrix that transforms a vector from frame 1 to frame 2. `h1` is a
vector resolved in frame 1. `h2` is the same vector in frame 2.
"""
resolve2(R, v2) = collect(R) * collect(v2)

"""
    h1 = resolve1(R12, h2)

`R` is a 3x3 matrix that transforms a vector from frame 1 to frame 2. `h2` is a
vector resolved in frame 2. `h1` is the same vector in frame 1.
"""
resolve1(R, v2) = collect(R)' * collect(v2)
