module Visualizers
using ModelingToolkit
using LinearAlgebra
using ..Multibody
using Multibody: _norm, _normalize, resolve1
t = Multibody.t

# export CylinderPrimitive
# export Cylinder, Frame

@component function CylinderPrimitive(;
    name,
    render = true,
    color = [1,0,0,1],
    specular_coefficient = 1,
    radius = 0.01,
)
    systems = @named begin
        frame_a = Frame()
    end
    pars = @parameters begin
        render = render
        color[1:4] = color
        specular_coefficient = specular_coefficient
        radius = radius
    end
    pars = reduce(vcat, collect.(pars))
    vars = @variables begin
        r(t)[1:3], [description="Vector from frame_a along the length direction of the cylinder, resolved in frame_a"]
    end
    vars = reduce(vcat, collect.(vars))
    eqs = Equation[
        frame_a.f .~ 0;
        frame_a.tau .~ 0;
    ]
    ODESystem(eqs, t, vars, pars; name, systems)
end

@component function Cylinder(;
    name,
    render = true,
    color = [1,0,0,1],
    specular_coefficient = 1,
    radius = 0.01,
)
    pars = @parameters begin
        # render = render
        # color[1:4] = color
        # specular_coefficient = specular_coefficient
        # radius = radius
    end
    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
        primitive = CylinderPrimitive(; render, color, specular_coefficient, radius)
    end
    vars = @variables begin
        r(t)[1:3], [description="Position vector directed from the origin of frame_a to the origin of frame_b, resolved in frame_a"]
    end
    vars = reduce(vcat, collect.(vars))
    eqs = Equation[
        collect(primitive.r) .~ collect(r);
        connect(frame_a, primitive.frame_a);

        frame_b.r_0 .~ collect(frame_a.r_0) + resolve1(ori(frame_a), r);
        ori(frame_b) ~ ori(frame_a);
      
        zeros(3) .~ collect(frame_a.f .+ frame_b.f);
        zeros(3) .~ collect(frame_a.tau .+ frame_b.tau .+ cross(r, frame_b.f));
    ]
    ODESystem(eqs, t, vars, pars; name, systems)
end

@component function FixedFrame(;
    name,
    render = false,
    color_x = [1,0,0,1],
    color_y = [0,1,0,1],
    color_z = [0,0,1,1],
    length = 1,
    specular_coefficient = 1,
    radius = 0.01,
    r_0 = nothing,
)

    pars = @parameters begin
        render = render
        # color_x[1:4] = color_x
        # color_y[1:4] = color_y
        # color_z[1:4] = color_z
        length = length
        specular_coefficient = specular_coefficient
        radius = radius
    end
    pars = reduce(vcat, collect.(pars))
    vars = @variables begin
    end
    
    systems = @named begin
        frame_a = Frame()
        arrow_x = CylinderPrimitive(; color=collect(color_x), radius, render, specular_coefficient)
        arrow_y = CylinderPrimitive(; color=collect(color_y), radius, render, specular_coefficient)
        arrow_z = CylinderPrimitive(; color=collect(color_z), radius, render, specular_coefficient)
    end
    eqs = Equation[
        connect(frame_a, arrow_x.frame_a, arrow_y.frame_a, arrow_z.frame_a);
        collect(arrow_x.r) .~ [length,0,0];
        collect(arrow_y.r) .~ [0,length,0];
        collect(arrow_z.r) .~ [0,0,length];
    ]
    ODESystem(eqs, t, vars, pars; name, systems)
end


end