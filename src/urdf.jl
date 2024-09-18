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

function parse_joint(xml_joint::XMLElement)
    name = getname(xml_joint)
    urdf_joint_type = attribute(xml_joint, "type")
    connections = ""
    if urdf_joint_type == "revolute" || urdf_joint_type == "continuous"
        axis = parse_vector(Float64, find_element(xml_joint, "axis"), "xyz", "1 0 0")
        damping = parse_scalar(Float64, find_element(xml_joint, "dynamics"), "damping", "0")
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
    # TODO: handle transformation of inertia
    mass, moment
end

getname(xml_link::XMLElement) = attribute(xml_link, "name")

function parse_body(graph, xml_link::XMLElement)
    xml_inertial = find_element(xml_link, "inertial")
    mass,inertia = xml_inertial == nothing ? (0,0*I(3)) : parse_inertia(xml_inertial)
    linkname = getname(xml_link)

    # Find link transformation by finding the any joint attached to the end of the link.
    # The transformation of the link is stored as the origin of the joint

    # Find all outgoing edges from this vertex
    connected_linknames = MetaGraphsNext.neighbor_labels(graph, linkname) |> collect
    length(connected_linknames) > 2 && error("Too many outgoing links from $linkname, this is not yet handled")

    joint = graph[linkname, connected_linknames[end]]
    R, r = parse_pose(find_element(joint, "origin"))

    # R, r = parse_pose(find_element(xml_link, "visual", "origin"))


    color = parse_color(xml_link, "1 0 0 1")
    cylinder = find_element(xml_link, "visual", "geometry", "cylinder")
    radius = if cylinder === nothing
        0.1
    else
        parse_scalar(Float64, cylinder, "radius", "0.1")
    end
    if R != I
        @warn "Ignoring rotation of link $linkname"
    end
    "$(Symbol(linkname)) = BodyShape(r=$(r), m=$(mass), I_11 = $(inertia[1,1]), I_22 = $(inertia[2,2]), I_33 = $(inertia[3,3]), I_21 = $(inertia[2,1]), I_31 = $(inertia[3,1]), I_32 = $(inertia[3,2]), color=$(color), radius=$(radius), sparse_I=true)"
end

"""

Example usage:
```
parse_urdf(joinpath(dirname(pathof(Multibody)), "..", "test/doublependulum.urdf"), extras=true, out="/tmp/urdf_import.jl")
```
"""
function parse_urdf(filename::AbstractString; extras=false, out=nothing)
    gravity = 9.81
    floating::Bool = false

    xdoc = parse_file(filename)
    xroot = LightXML.root(xdoc)
    @assert LightXML.name(xroot) == "robot"

    xml_links = get_elements_by_tagname(xroot, "link")
    xml_joints = get_elements_by_tagname(xroot, "joint")

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


    # Parse all joints and possible extra connections due to axisflanges
    joints_extraconnections = map(xml_joints) do l
        parse_joint(l)
    end
    joints = first.(joints_extraconnections)
    extra_connections = filter(!isempty, last.(joints_extraconnections))

    # Parse all bodies
    bodies = map(xml_links) do l
        parse_body(graph, l)
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

    s = @eval Main (using JuliaFormatter; JuliaFormatter.format_text($s, align_assignment=true, align_struct_field=true, align_conditional=true, align_pair_arrow=true))

    out === nothing && return s
    write(out, s)
    s
end

function outer_model(f, io, name="URDFModel")
    println(io, "@mtkmodel $name begin")
    f(io)
    println(io, "end")
end

# filename = "/home/fredrikb/.julia/dev/Multibody/test/doublependulum.urdf"

