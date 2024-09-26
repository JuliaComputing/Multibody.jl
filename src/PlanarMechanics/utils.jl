@connector function Frame(; name, render=false, length=1.0, radius=0.1)
    vars = @variables begin
        x(t), [state_priority = -1, description = "x position"]
        y(t), [state_priority = -1, description = "y position"]
        phi(t), [state_priority = 0, description = "rotation angle (counterclockwise)"]
        fx(t), [connect = Flow, description = "force in the x direction"]
        fy(t), [connect = Flow, description = "force in the y direction"]
        tau(t), [connect = Flow, description = "torque (clockwise)"]
    end
    pars = @parameters begin
        render = render, [description = "Render the joint in animations"]
        length = length, [description = "Length of each axis in animations"]
        radius = radius, [description = "Radius of each axis in animations"]
    end

    ODESystem(Equation[], t, vars, pars; name, metadata = Dict(:frame_2d => true))
end
Base.@doc """
    Frame(;name)

Coordinate system (2-dim.) fixed to the component with one cut-force and cut-torque

# Variables:
- `x`: [m] x position
- `y`: [m] y position
- `phi`: [rad] rotation angle (counterclockwise)
- `fx`: [N] force in the x direction
- `fy`: [N] force in the y direction
- `tau`: [N.m] torque (clockwise)

# Parameters:
- `render`: [Bool] Render the joint in animations
- `length`: [m] Length of each axis in animations
- `radius`: [m] Radius of each axis in animations
""" Frame

function ori_2d(frame)
    phi = frame.phi
    return [cos(phi) -sin(phi); sin(phi) cos(phi)]
end

FrameResolve = Frame

@mtkmodel PartialTwoFrames begin
    @components begin
        frame_a = Frame()
        frame_b = Frame()
    end
end

Base.@doc """
    PartialTwoFrames(;name)
Partial model with two frames

# Connectors:
- `frame_a` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
- `frame_b` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
""" PartialTwoFrames

"""
    ZeroPosition(;name)
Set zero position vector and orientation object of frame_resolve

# Connectors:
- `frame_resolve` [FrameResolve](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
"""
@mtkmodel ZeroPosition begin
    @components begin
        frame_resolve = FrameResolve()
    end

    @equations begin
        frame_resolve.x ~ 0
        frame_resolve.y ~ 0
        frame_resolve.phi ~ 0
    end
end
