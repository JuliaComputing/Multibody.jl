module URDF
using Multibody
using LightXML
using StaticArrays
using Rotations
using LinearAlgebra: ×, I, norm, normalize
using Graphs
using MetaGraphsNext
using JuliaFormatter

function LightXML.find_element(x::XMLElement, n::AbstractString, ns::AbstractString...)
    inner = find_element(x, n)
    inner === nothing && return nothing
    find_element(inner, ns...)
end

function parse_scalar(::Type{T}, e::XMLElement, name::String) where {T}
    parse(T, attribute(e, name))
end

function parse_scalar(::Type{T}, e, name::String, default::String) where {T}
    parse(T, e == nothing ? default : attribute(e, name))
end

function parse_vector(::Type{T}, e::Union{XMLElement, Nothing}, name::String, default::String) where {T}
    usedefault = e == nothing || attribute(e, name) == nothing
    [parse(T, str) for str in split(usedefault ? default : attribute(e, name))]
end

function parse_color(xml_link, default)
    color_elem = find_element(xml_link, "visual", "material", "color")
    parse_vector(Float64, color_elem, "rgba", "1 0 0 1")
end

function parse_inertia_mat(xml_inertia::XMLElement)
    ixx = parse_scalar(Float64, xml_inertia, "ixx", "0")
    ixy = parse_scalar(Float64, xml_inertia, "ixy", "0")
    ixz = parse_scalar(Float64, xml_inertia, "ixz", "0")
    iyy = parse_scalar(Float64, xml_inertia, "iyy", "0")
    iyz = parse_scalar(Float64, xml_inertia, "iyz", "0")
    izz = parse_scalar(Float64, xml_inertia, "izz", "0")
    [ixx ixy ixz; ixy iyy iyz; ixz iyz izz]
end

function parse_pose(xml_pose::Nothing)
    rot = one(RotMatrix3{Float64})
    trans = zeros(3)
    rot, trans
end

function parse_pose(xml_pose::XMLElement)
    rpy = parse_vector(Float64, xml_pose, "rpy", "0 0 0")
    rot = RotMatrix(RotZYX(rpy[3], rpy[2], rpy[1]))
    trans = parse_vector(Float64, xml_pose, "xyz", "0 0 0")
    rot, trans
end

function parse_joint(xml_joint::XMLElement, render_fixed, render_joints)
    urdf_joint_type = attribute(xml_joint, "type")
    name = getname(xml_joint)

    R, r = parse_pose(find_element(xml_joint, "origin"))

    connections = ""
    if urdf_joint_type == "revolute" || urdf_joint_type == "continuous"
        axis = parse_vector(Float64, find_element(xml_joint, "axis"), "xyz", "1 0 0")
        damping = parse_scalar(Float64, find_element(xml_joint, "dynamics"), "damping", "0")
        if iszero(damping)
            components = "$name = URDFRevolute(; r=$r, R=$R, n=$axis)"
        else
            components = """
            $name = URDFRevolute(; r=$r, R=$R, n=$axis, axisflange=true, render=$render_joints)
            $(name)_damper = Rotational.Damper(; d=$damping)
            """
            connections = """
            connect($name.support, $(name)_damper.flange_a)
            connect($(name)_damper.flange_b, $name.axis)
            """
        end
    elseif urdf_joint_type == "prismatic"
        axis = parse_vector(T, find_element(xml_joint, "axis"), "xyz", "1 0 0")
        damping = parse_scalar(Float64, find_element(xml_joint, "dynamics"), "damping", "0")
        if iszero(damping)
            components = "$name = URDFPrismatic(; r=$r, R=$R, n=$axis, render=$render_joints)"
        else
            components = """
            $name = URDFPrismatic(; r=$r, R=$R, n=$axis, axisflange=true)
            $(name)_damper = Translational.Damper(; d=$damping)
            """
            connections = """
            connect($name.support, $(name)_damper.flange_a)
            connect($(name)_damper.flange_b, $name.axis)
            """
        end
    elseif urdf_joint_type == "floating"
        components = "$name = FreeMotion()"
        connections = "connect(world.frame_b, $name.frame_a)"
    elseif urdf_joint_type == "fixed"
        if norm(r) == 0 && R == I
            components = "$name = NullJoint()" # Null joint
        elseif R == I(3)
            components = "$name = FixedTranslation(; r=$r, render=$render_fixed)"
        else
            components = "$name = FixedTranslation(; r=$r, render=$render_fixed)"
            @warn "Ignoring rotation of joint $name"
            # R = RotMatrix3(R)
            # components = "$name = FixedRotation(; r=$r, n = $(rotation_axis(R)), angle = $(rotation_angle(R)))"
        end
    elseif urdf_joint_type == "planar"
        urdf_axis = parse_vector(Float64, find_element(xml_joint, "axis"), "xyz", "1 0 0")
        # The URDF spec says that a planar joint allows motion in a
        # plane perpendicular to the axis.
        R = Rotations.rotation_between([0, 0, 1], urdf_axis)
        x_axis = R * [1, 0, 0]
        components = "$name = $(joint_type)(; n_x=$x_axis, n=$urdf_axis)"
    else
        error("joint type $(urdf_joint_type) not recognized")
    end
    components, connections
end


function parse_inertia(xml_inertial::XMLElement)
    moment = parse_inertia_mat(find_element(xml_inertial, "inertia"))
    mass = parse_scalar(Float64, find_element(xml_inertial, "mass"), "value", "0")
    r_cm = parse_vector(Float64, find_element(xml_inertial, "origin"), "xyz", "0 0 0")
    rpy = parse_vector(Float64, find_element(xml_inertial, "origin"), "rpy", "0 0 0")
    if !iszero(rpy)
        R = RotMatrix(RotXYZ(rpy[1], rpy[2], rpy[3]))
        moment = R * moment * R'
        # TODO: Double-check and test inertia transform, rotation convention RotXYZ? Transformation RIR' or R'IR?
    end

    mass, moment, r_cm
end

getname(xml_link::XMLElement) = attribute(xml_link, "name")

function parse_geometry(xml_link::XMLElement)
    xml_geometry = find_element(xml_link, "visual", "geometry")
    if xml_geometry === nothing
        @error "No geometry found for link $(getname(xml_link)), using sphere"
        return (; radius=0.1), :sphere
    elseif (cylinder = find_element(xml_geometry, "cylinder")) !== nothing
        radius = parse_scalar(Float64, cylinder, "radius", "0.1")
        length = parse_scalar(Float64, cylinder, "length", "1")
        geometry = (; radius, length)
        type = :cylinder
    elseif (box = find_element(xml_geometry, "box")) !== nothing
        size = parse_vector(Float64, box, "size", "1 1 1")
        geometry = (; size)
        type = :box

    elseif (sphere = find_element(xml_geometry, "sphere")) !== nothing
        radius = parse_scalar(Float64, sphere, "radius", "0.1")
        geometry = (; radius)
        type = :sphere
    elseif (mesh = find_element(xml_geometry, "mesh")) !== nothing
        @warn "Mesh geometry is not yet supported, using cylinder instead"
        # filename = attribute(mesh, "filename")
        # scale = parse_scalar(Float64, mesh, "scale", "1")
        radius = 0.1
        length = 0.1
        geometry = (; radius, length)
        type = :cylinder
    end
    
    return geometry, type
end

function parse_body(graph, xml_link::XMLElement; min_mass)
    xml_inertial = find_element(xml_link, "inertial")
    mass,inertia,r_cm = xml_inertial == nothing ? (0.0,0.0*I(3),zeros(3)) : parse_inertia(xml_inertial)
    if min_mass > 0
        mass = max(mass, min_mass)
        for i = 1:3
            inertia[i,i] = max(inertia[i,i], min_mass*0.01)
        end
    end
    linkname = getname(xml_link)

    R, r = parse_pose(find_element(xml_link, "visual", "origin"))
    color = parse_color(xml_link, "1 0 0 1")
    geometry, type = parse_geometry(xml_link)
    if geometry === nothing
        error("No geometry found for link $linkname")
    end
    if R != I
        @warn "Ignoring rotation of link $linkname"
    end
    if mass == 0
        # We special case this since the dynamics becomes strange with zero mass
        if type === :sphere
            return "$(Symbol(linkname)) = SphereVisualizer(; radius=$(geometry.radius))" # color=$(color),
        elseif type === :cylinder
            if iszero(r)
                r = [1, 0, 0]
            end
            return "$(Symbol(linkname)) = CylinderVisualizer(; length_direction=$(r), radius=$(geometry.radius), length=$(geometry.length))" # color=$(color),
        elseif type === :box
            return "$(Symbol(linkname)) = BoxVisualizer(; length_direction=[1,0,0], width_direction=[0,1,0], length=$(geometry.size[1]), width=$(geometry.size[2]), height=$(geometry.size[3]))" # color=$(color),
        end

    elseif type === :sphere
        radius = geometry.radius
        "$(Symbol(linkname)) = Body(; m=$(mass), r_cm=$(r_cm), I_11 = $(inertia[1,1]), I_22 = $(inertia[2,2]), I_33 = $(inertia[3,3]), I_21 = $(inertia[2,1]), I_31 = $(inertia[3,1]), I_32 = $(inertia[3,2]), color=$(color), radius=$(radius), sparse_I=true)"
    elseif type === :cylinder
        (; radius, length) = geometry
        if iszero(r)
            r = [1, 0, 0]
        end
        r = normalize(r)*length
        "$(Symbol(linkname)) = BodyShape(; r=$(r), m=$(mass), r_cm=$(r_cm), I_11 = $(inertia[1,1]), I_22 = $(inertia[2,2]), I_33 = $(inertia[3,3]), I_21 = $(inertia[2,1]), I_31 = $(inertia[3,1]), I_32 = $(inertia[3,2]), color=$(color), radius=$(radius), sparse_I=true)"
    elseif type === :box
        (; size) = geometry
        length = size[1]
        width = size[2]
        height = size[3]
        if iszero(r)
            r = [1, 0, 0]
        end
        "$(Symbol(linkname)) = BodyBox(; r=$(r), body.m=$(mass), body.r_cm=$(r_cm), body.I_11 = $(inertia[1,1]), body.I_22 = $(inertia[2,2]), body.I_33 = $(inertia[3,3]), body.I_21 = $(inertia[2,1]), body.I_31 = $(inertia[3,1]), body.I_32 = $(inertia[3,2]), color=$(color), length=$(length), width=$(width), height=$(height))"
    else
        error("Geometry type $type not recognized")
    end
end

default_modelname(filename) = uppercasefirst(splitext(basename(filename))[1])

"""
    urdf2multibody(filename::AbstractString; extras=false, out=nothing, worldconnection = :rigid)


Example usage:
```
urdf2multibody(joinpath(dirname(pathof(Multibody)), "..", "test", "urdf", "doublependulum.urdf"), extras=true, out="/tmp/urdf_import.jl")
```

## Keyword arguments
- `extras=false`: If `true`, the generated code will include package imports, a simulation of the model and a rendering of the model.
- `out=nothing`: If provided, the generated code will be written to this file, otherwise the translated model is returned as a string.
- `worldconnection=:rigid`: If `:rigid`, the world frame will be connected to the root link with a rigid connection. If a joint constructor is provided, this component will be instantiated and the root link is connected to the world through this, e.g., `worldconnection = FreeMotion`, `()->Prismatic(n=[0, 1, 0])` etc.
- `modelname`: The name of the model, default is the name of the URDF file with the extension removed and the first letter upper case.
"""
function Multibody.urdf2multibody(filename::AbstractString; extras=false, out=nothing, worldconnection = :rigid, modelname = default_modelname(filename), solver = "FBDF", render_fixed=false, render_joints = true, min_mass=0)
    gravity = 9.81
    floating::Bool = false

    xdoc = parse_file(filename)
    xroot = LightXML.root(xdoc)
    @assert LightXML.name(xroot) == "robot"

    xml_links = get_elements_by_tagname(xroot, "link")
    xml_joints = get_elements_by_tagname(xroot, "joint")

    # create graph structure of XML elements
    graph = MetaGraph(DiGraph(), label_type=String, vertex_data_type=eltype(xml_links), edge_data_type=eltype(xml_joints))
    for (i,vertex) in enumerate(xml_links)
        # add_vertex!(graph)
        graph[attribute(vertex, "name")] = vertex
    end
    name_to_vertex = Dict(attribute(v, "name") => i for (i,v) in enumerate(xml_links))
    for xml_joint in xml_joints

        parent = attribute(find_element(xml_joint, "parent"), "link")
        child = attribute(find_element(xml_joint, "child"), "link")
        graph[parent, child] = xml_joint
    end

    cycles = simplecycles(graph)
    if !isempty(cycles)
        @info("The mechanism has kinematic loops. Kinematic loops must be broken manually, see https://juliacomputing.github.io/Multibody.jl/dev/examples/kinematic_loops/ for more info. The following loops were found", cycles)
    end


    # Parse all joints and possible extra connections due to axisflanges
    joints_extraconnections = map(xml_joints) do l
        parse_joint(l, render_fixed, render_joints)
    end
    joints = first.(joints_extraconnections)
    extra_connections = filter(!isempty, last.(joints_extraconnections))

    # Parse all bodies
    bodies = map(xml_links) do l
        parse_body(graph, l; min_mass)
    end

    roots = [v for v in vertices(graph) if indegree(graph, v) == 0]

    connections = map(edges(graph)) do e
        src_link = label_for(graph, e.src)
        dst_link = label_for(graph, e.dst)
        joint = attribute(graph[src_link, dst_link], "name")
        """connect($(src_link).frame_a, $(joint).frame_a)
        connect($(joint).frame_b, $(dst_link).frame_a)"""

    end

    for root in roots
        if worldconnection == :rigid
            pushfirst!(connections, "connect(world.frame_b, $(getname(xml_links[root])).frame_a)")
        elseif worldconnection isa Function
            push!(joints, "world_connection = $(nameof(worldconnection))()")
            push!(connections, "connect(world.frame_b, world_connection.frame_a)")
            push!(connections, "connect(world_connection.frame_b, $(getname(xml_links[root])).frame_a)")
        end
    end

    s = if extras
        """
        using ModelingToolkit, Multibody, JuliaSimCompiler, OrdinaryDiffEq, Plots
        import ModelingToolkit: t_nounits as t, D_nounits as D
        """
    else 
        ""
    end
    s = s * """
    @mtkmodel $(modelname) begin
        @components begin
            world = World()
            $(join(bodies, "\n"))
            $(join(joints, "\n"))
        end
        @equations begin
            $(join(connections, "\n"))
            $(join(extra_connections, "\n"))
        end
    end
    """
    if extras
        s = s * """
        @named model = $(modelname)()
        model = complete(model)
        ssys = structural_simplify(multibody(model))
        prob = ODEProblem(ssys, [], (0.0, 10.0))
        sol = solve(prob, $(solver)())
        plot(sol) |> display

        import GLMakie
        first(Multibody.render(model, sol, 0, show_axis=true))
        """
    end

    s = JuliaFormatter.format_text(s, align_assignment=true, align_struct_field=true, align_conditional=true, align_pair_arrow=true)

    out === nothing && return s
    write(out, s)
    s
end


end