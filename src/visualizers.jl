@component function Visaluzable(; name, color = [1.0, 0.0, 0.0, 1.0], render = true)
    pars = @parameters begin
        color[1:4] = color, [description = "Color of the body in animations"]
        render = render, [description = "Render the joint in animations"]
    end

    systems = @named begin
        frame_a = Frame()
    end

    vars = @variables begin
    end

    equations = [
        frame_a.f ~ zeros(3)
        frame_a.tau ~ zeros(3)
    ]

    return System(equations, t; name, systems)
end

@component function SphereVisualizer(; name, color = [1.0, 0.0, 0.0, 1.0], render = true, radius = 0.1)
    @named v = Visaluzable(; color, render)

    pars = @parameters begin
        radius = radius, [description = "Radius of the sphere in animations"]
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = [
    ]

    return extend(System(equations, t; name, systems), v)
end

@component function CylinderVisualizer(; name, color = [1.0, 0.0, 0.0, 1.0], render = true,
                                       radius = 0.1, length = 1.0, length_direction = [1.0, 0.0, 0.0])
    @named v = Visaluzable(; color, render)

    pars = @parameters begin
        radius = radius, [description = "Radius of the cylinder in animations"]
        length = length, [description = "Length of the cylinder"]
        length_direction[1:3] = length_direction, [description = "Length direction of the cylinder, resolved in frame_a"]
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = [
    ]

    return extend(System(equations, t; name, systems), v)
end

@component function BoxVisualizer(; name, color = [1.0, 0.0, 0.0, 1.0], render = true,
                                  length_direction = [1.0, 0.0, 0.0], width_direction = [0.0, 1.0, 0.0],
                                  length = 1.0, width = 1.0, height = 1.0, r_shape = [0, 0, 0])
    @named v = Visaluzable(; color, render)

    pars = @parameters begin
        length_direction[1:3] = length_direction, [description = "Length direction of the box, resolved in frame_a"]
        width_direction[1:3] = width_direction, [description = "Width direction of the box, resolved in frame_a"]
        length = length, [description = "Length of the box"]
        width = width, [description = "Width of the box"]
        height = height, [description = "Height of the box"]
        r_shape[1:3] = r_shape, [description = "Position vector from frame_a to box center"]
    end

    systems = @named begin
    end

    vars = @variables begin
    end

    equations = [
    ]

    return extend(System(equations, t; name, systems), v)
end

