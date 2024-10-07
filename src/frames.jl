@connector function Frame(; name, varw = false, r_0 = zeros(3), render=false, length=1.0, radius=0.1)
    @variables r_0(t)[1:3]=r_0 [
        description = "Position vector directed from the origin of the world frame to the connector frame origin, resolved in world frame",
        state_priority=-1,
    ]
    @variables f(t)[1:3] [
        connect = Flow,
        description = "Cut force resolved in connector frame",
    ]
    @variables tau(t)[1:3] [
        connect = Flow,
        description = "Cut torque resolved in connector frame",
    ]
    pars = @parameters begin
        render = render
        length = length
        radius = radius
    end
    r_0, f, tau = collect.((r_0, f, tau))
    # R: Orientation object to rotate the world frame into the connector frame

    R = NumRotationMatrix(; name, varw, state_priority=-1)

    ODESystem(Equation[], t, [r_0; f; tau], pars; name,
              metadata = Dict(:orientation => R, :frame => true))
end

"""
    Frame(; name)

`Frame` is the fundamental 3D connector in the multibody library. Most components have one or several `Frame` connectors that can be connected together.

The `Frame` connector has internal variables for
- `r_0`: The position vector from the world frame to the frame origin, resolved in the world frame
- `f`: The cut force resolved in the connector frame
- `tau`: The cut torque resolved in the connector frame
- Depending on usage, also rotation and rotational velocity variables, accessed using `ori(frame).R` and `ori(frame).w`. The rotation matrix represents the rotation to rotate the world frame into the connector frame

# Parameters
- `render = false`: If true, the connector is rendered in animations
- `length = 1.0`: Length of the frame when rendered
- `radius = 0.1`: Radius of the frame when rendered
"""
Frame
