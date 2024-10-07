@mtkmodel Visaluzable begin
    @components begin
        frame_a = Frame()
    end
    @parameters begin
        color[1:4] = [1.0, 0.0, 0.0, 1.0], [description = "Color of the body in animations"]
        render = true, [description = "Render the joint in animations"]
    end
    @equations begin
        frame_a.f ~ zeros(3)
        frame_a.tau ~ zeros(3)
    end
end

@mtkmodel SphereVisualizer begin
    @extend () = v = Visaluzable()
    @parameters begin
        radius = 0.1, [description = "Radius of the sphere in animations"]
    end
end

@mtkmodel CylinderVisualizer begin
    @extend () = v = Visaluzable()
    @parameters begin
        radius = 0.1, [description = "Radius of the sphere in animations"]
        length = 1.0, [description = "length of the cylinder"]
        length_direction[1:3] = [1.0, 0.0, 0.0], [description = "length direction of the cylinder, resolved in frame_a"]
    end
end

@mtkmodel BoxVisualizer begin
    @extend () = v = Visaluzable()
    @parameters begin
        length_direction[1:3] = [1.0, 0.0, 0.0], [description = "length direction of the box, resolved in frame_a"]
        width_direction[1:3] = [0.0, 1.0, 0.0], [description = "width direction of the box, resolved in frame_a"]
        length = 1.0, [description = "length of the box"]
        width = 1.0, [description = "width of the box"]
        height = 1.0, [description = "height of the box"]
        r_shape[1:3] = [0,0,0]
    end
end

