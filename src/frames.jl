function Oriantation(; name)
    # T: Transformation matrix from world frame to local frame
    # w: Absolute angular velocity of local frame, resolved in local frame

    @variables R(t)[1:3, 1:3]=collect(1.0I(3)) [description="Orientation rotation matrix ∈ SO(3)"]
    @variables w(t)[1:3] [description="angular velocity"]
    R,w = collect.((R,w))

    ODESystem(Equation[], t, [vec(R); w], [], name = name)
end

@connector function Frame(; name)
    @variables r_0(t)[1:3] [description = "Position vector directed from the origin of the world frame to the connector frame origin, resolved in world frame"]
    @variables f(t)[1:3] [connect = Flow, description = "Cut force resolved in connector frame"]
    @variables tau(t)[1:3] [connect = Flow, description = "Cut torque resolved in connector frame"]
    r_0, f, tau = collect.((r_0, f, tau))
    # R: Orientation object to rotate the world frame into the connector frame
    @named R = Oriantation()

    compose(ODESystem(Equation[], t, [r_0; f; tau], [], name = name), R)
end

# TODO: set angular velocity to zero as well
# Should `~` be defined for `System ~ System` and `System ~ NamedTuple`?
nullrotation() = Matrix(1.0I(3))

"""
    h2 = resolve2(R12, h1)

`R` is a 3x3 matrix that transforms a vector from frame 1 to frame 2. `h1` is a
vector resolved in frame 1. `h2` is the same vector in frame 2.
"""
resolve2(R::AbstractMatrix, v2) = collect(R) * collect(v2)
resolve2(R, v2) = resolve2(R.R, v2)

"""
    h1 = resolve1(R12, h2)

`R` is a 3x3 matrix that transforms a vector from frame 1 to frame 2. `h2` is a
vector resolved in frame 2. `h1` is the same vector in frame 1.
"""
resolve1(R, v2) = collect(R)' * collect(v2)

orientation_constraint(q::AbstractVector) = q'q - 1

function angular_velocity2(q::AbstractVector, q̇) 
    Q = [q[4] q[3] -q[2] -q[1]; -q[3] q[4] q[1] -q[2]; q[2] -q[1] q[4] -q[3]]
    2*Q*q̇
end
