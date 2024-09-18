using LightXML
using StaticArrays
using Rotations
using LinearAlgebra: ×, I
using Graphs
using MetaGraphsNext

@mtkmodel NullJoint begin
    @components begin
        frame_a = Frame()
        frame_b = Frame()
    end
    @equations begin
        connect(frame_a, frame_b)
    end
end

function LightXML.find_element(x::XMLElement, n::AbstractString, ns::AbstractString...)
    inner = find_element(x, n)
    inner === nothing && return nothing
    find_element(inner, ns...)
end

function parse_scalar(::Type{T}, e::XMLElement, name::String) where {T}
    parse(T, attribute(e, name))
end

function parse_scalar(::Type{T}, e::XMLElement, name::String, default::String) where {T}
    parse(T, e == nothing ? default : attribute(e, name))
end

function parse_vector(::Type{T}, e::Union{XMLElement, Nothing}, name::String, default::String) where {T}
    usedefault = e == nothing || attribute(e, name) == nothing # TODO: better handling of required attributes
    [parse(T, str) for str in split(usedefault ? default : attribute(e, name))]
end

function parse_color(xml_link, default)
    color_elem = find_element(xml_link, "visual", "material", "color")
    parse_vector(Float64, color_elem, "rgba", "1 0 0 1")
end

function parse_inertia_mat(::Type{T}, xml_inertia::XMLElement) where {T}
    ixx = parse_scalar(T, xml_inertia, "ixx", "0")
    ixy = parse_scalar(T, xml_inertia, "ixy", "0")
    ixz = parse_scalar(T, xml_inertia, "ixz", "0")
    iyy = parse_scalar(T, xml_inertia, "iyy", "0")
    iyz = parse_scalar(T, xml_inertia, "iyz", "0")
    izz = parse_scalar(T, xml_inertia, "izz", "0")
    [ixx ixy ixz; ixy iyy iyz; ixz iyz izz]
end

function parse_pose(::Type{T}, xml_pose::Nothing) where {T}
    rot = one(RotMatrix3{T})
    trans = zeros(3)
    rot, trans
end

function parse_pose(::Type{T}, xml_pose::XMLElement) where {T}
    rpy = parse_vector(T, xml_pose, "rpy", "0 0 0")
    rot = RotMatrix(RotZYX(rpy[3], rpy[2], rpy[1]))
    trans = parse_vector(T, xml_pose, "xyz", "0 0 0")
    rot, trans
end

function parse_joint_type(::Type{T}, xml_joint::XMLElement, name) where {T}
    urdf_joint_type = attribute(xml_joint, "type")
    joint_type = joint_types[urdf_joint_type]
    # TODO: add Damper from <dynamics damping="0.1" />
    # TODO: add color information
    connections = ""
    if urdf_joint_type == "revolute" || urdf_joint_type == "continuous"
        axis = parse_vector(T, find_element(xml_joint, "axis"), "xyz", "1 0 0")
        damping = parse_scalar(T, find_element(xml_joint, "dynamics"), "damping", "0")
        if damping == "0"
            components = "$name = Revolute(; n=$axis)"
        else
            components = """
            $name = Revolute(; n=$axis, axisflange=true)
            $(name)_damper = Rotational.Damper(; d=$damping)
            """
            connections = """
            connect($name.support, $(name)_damper.flange_a)
            connect($(name)_damper.flange_b, $name.axis)
            """
        end
    elseif urdf_joint_type == "prismatic"
        axis = parse_vector(T, find_element(xml_joint, "axis"), "xyz", "1 0 0")
        components = "$name = Prismatic(; n=$axis)"
    elseif urdf_joint_type == "floating"
        components = "$name = FreeMotion()"
        connections = "connect(world.frame_b, $name.frame_a)"
    elseif urdf_joint_type == "fixed"
        components = "$name = NullJoint()" # Null joint
    elseif urdf_joint_type == "planar"
        urdf_axis = parse_vector(T, find_element(xml_joint, "axis"), "xyz", "1 0 0")
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

function parse_joint(::Type{T}, xml_joint::XMLElement) where {T}
    name = attribute(xml_joint, "name")
    parse_joint_type(T, xml_joint, name)
    # position_bounds, velocity_bounds, effort_bounds = parse_joint_bounds(joint_type, xml_joint)
end

function parse_inertia(::Type{T}, xml_inertial::XMLElement) where {T}
    moment = parse_inertia_mat(T, find_element(xml_inertial, "inertia"))
    mass = parse_scalar(T, find_element(xml_inertial, "mass"), "value", "0")
    # TODO: handle transformation of inertia
    mass, moment
end

function parse_body(::Type{T}, xml_link::XMLElement) where {T}
    xml_inertial = find_element(xml_link, "inertial")
    mass,inertia = xml_inertial == nothing ? (0,0*I(3)) : parse_inertia(T, xml_inertial)
    linkname = attribute(xml_link, "name")
    R, r = parse_pose(T, find_element(xml_link, "visual", "origin"))

    @show color = parse_color(xml_link, "1 0 0 1")
    cylinder = find_element(xml_link, "visual", "geometry", "cylinder")
    radius = if cylinder === nothing
        0.1
    else
        parse_scalar(T, cylinder, "radius", "0.1")
    end
    if R != I
        @warn "Ignoring rotation of link $linkname"
    end
    "$(Symbol(linkname)) = BodyShape(r=$(r), m=$(mass), I_11 = $(inertia[1,1]), I_22 = $(inertia[2,2]), I_33 = $(inertia[3,3]), I_21 = $(inertia[2,1]), I_31 = $(inertia[3,1]), I_32 = $(inertia[3,2]), color=$(color), radius=$(radius))"
end


function parse_urdf(filename::AbstractString; extras=false, out=nothing)
    gravity = 9.81
    floating::Bool = false
    joint_types = default_urdf_joint_types()
    root_joint_type = joint_types[floating ? "floating" : "fixed"](name=:root)

    xdoc = parse_file(filename)
    xroot = LightXML.root(xdoc)
    @assert LightXML.name(xroot) == "robot"

    xml_links = get_elements_by_tagname(xroot, "link")
    xml_joints = get_elements_by_tagname(xroot, "joint")
    poses = map(xml_joints) do xml_joint
        pose = parse_pose(Float64, find_element(xml_joint, "origin"))
    end

    # create graph structure of XML elements
    graph = MetaGraph(Graph(), label_type=String, vertex_data_type=eltype(xml_links), edge_data_type=eltype(xml_joints))
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

    joints_extraconnections = map(xml_joints) do l
        parse_joint(Float64, l, joint_types)
    end
    joints = first.(joints_extraconnections)
    extra_connections = filter(!isempty, last.(joints_extraconnections))

    bodies = map(xml_links) do l
        parse_body(Float64, l)
    end

    connections = map(edges(graph)) do e
        src_link = label_for(graph, e.src)
        dst_link = label_for(graph, e.dst)
        joint = attribute(graph[src_link, dst_link], "name")
        """
        connect($(src_link).frame_b, $(joint).frame_a)
        connect($(joint).frame_b, $(dst_link).frame_a)
        """

    end

    s = if extras
        """
        using ModelingToolkit, Multibody, JuliaSimCompiler, OrdinaryDiffEq, Plots
        W(args...; kwargs...) = Multibody.world
        """
    else 
        ""
    end
    s = s * """
    @mtkmodel URDFModel begin
        @components begin
            world = W()
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
        @named model = URDFModel()
        model = complete(model)
        ssys = structural_simplify(IRSystem(model))
        prob = ODEProblem(ssys, [], (0.0, 10.0))
        sol = solve(prob, FBDF())
        plot(sol)
        """
    end
    out === nothing && return s

    write(out, s)

    @eval Main (using JuliaFormatter; JuliaFormatter.format_file($out, align_assignment=true, align_struct_field=true, align_conditional=true, align_pair_arrow=true))

    s

end

function outer_model(f, io, name="URDFModel")
    println(io, "@mtkmodel $name begin")
    f(io)
    println(io, "end")
end

# filename = "/home/fredrikb/.julia/dev/Multibody/test/doublependulum.urdf"

