module Render
using Makie
using Multibody
import Multibody: render, render!, loop_render, encode, decode, get_rot, get_trans, get_frame, get_systemtype
import Multibody.PlanarMechanics as P
using Rotations
using LinearAlgebra
using ModelingToolkit
using ModelingToolkit: get_metadata
export render, loop_render
using MeshIO, FileIO
using StaticArrays


import Multibody: get_rot_fun, get_fun, get_frame_fun, get_color, get_shapefile, get_shape

"""
This struct is used to mimic a solution object such that a model can be rendered without having an actual solution
"""
struct FakeSol
    model
    prob
    function FakeSol(model, prob)
        new(model, prob)
    end
end


import ModelingToolkit.getu

function recursive_extract(sol, s, depth=0)
    depth > 10 && error("Recursion depth exceeded with symbol $s")
    (; model, prob) = sol
    v = if ModelingToolkit.isparameter(s)
        prob.ps[s]
    else
        prob[s]
    end

    if ModelingToolkit.Symbolics.is_symbolic_or_array_of_symbolic(v)
        return recursive_extract(sol, v, depth+1)
    else
        return v
    end
end

function getu(sol::FakeSol, syms)
    t->[recursive_extract(sol, s) for s in syms]
end

function ModelingToolkit.parameter_values(sol::FakeSol)
    pars = Multibody.collect_all(parameters(sol.model))
    
    sol.prob.ps[pars]
end

function (sol::FakeSol)(t; idxs=nothing)
    if idxs === nothing
        sol.prob[unknowns(sol.model)]
    elseif idxs isa Real
        recursive_extract(sol, idxs)
    else
        # ret = zeros(length(idxs))
        # Threads.@threads for i in 1:length(idxs)
        #     ret[i] = recursive_extract(sol, idxs[i])
        # end
        # ret
        [recursive_extract(sol, i) for i in idxs]
    end
end

function Base.getproperty(sol::FakeSol, s::Symbol)
    s ∈ fieldnames(typeof(sol)) && return getfield(sol, s)
    if s === :t
        return [0.0]
    else
        throw(ArgumentError("$(typeof(sol)) has no property named $s"))
    end
end

mutable struct CacheSol
    model
    sol
    cache::Dict{Any, Any}
    vars
    last_t::Float64
    function CacheSol(model, sol)
        vars = get_all_vars(model, sol) |> unique
        vars = filter(v->!occursin("₊tau(t)", string(v)), vars)
        vars = filter(v->!occursin("₊a_0(t)", string(v)), vars)
        vars = filter(v->!occursin("₊w_a(t)", string(v)), vars)
        vars = filter(v->!occursin("₊v_0(t)", string(v)), vars)
        vars = filter(v->!occursin("₊u(t)", string(v)), vars)
        vars = filter(v->!occursin("₊z_a(t)", string(v)), vars)
        vars = filter(v->!occursin("₊phid(t)", string(v)), vars)
        # Main.vars = vars
        # @show length(vars)
        # filter!(v->ModelingToolkit.SymbolicIndexingInterface.is_variable(sol, v), vars) # To work around https://github.com/SciML/ModelingToolkit.jl/issues/3065
        # @show length(vars)
        values = sol(0.0, idxs=vars)
        new(model, sol, Dict(vars .=> values), vars, 0)
    end
end

function get_all_vars(model, sol, vars = Multibody.collect_all(unknowns(model)))
    @parameters render
    nsrender = Multibody.ModelingToolkit.renamespace(getfield(model, :name), render)
    dorender = try
        Bool(sol.prob.ps[nsrender])
    catch
        true
    end
    for sys in getfield(model, :systems)
        if ModelingToolkit.isframe(sys)
            dorender || continue
            newvars = Multibody.ModelingToolkit.renamespace.(getfield(model, :name), Multibody.Symbolics.unwrap.(vec(ori(sys).R)))
            append!(vars, newvars)
        else
            subsys_ns = getproperty(model, getfield(sys, :name))
            get_all_vars(subsys_ns, sol, vars)
        end
    end
    vars
end


get_cached(cs::CacheSol, t::AbstractArray, idxs) = cs.sol(t; idxs)
function get_cached(cs::CacheSol, t::Real, idxs)
    if ModelingToolkit.isparameter(idxs[1])
        return cs.prob.ps[idxs]
    end
    if idxs isa AbstractArray{Num}
        idxs = Multibody.Symbolics.unwrap.(idxs)
    end
    if !haskey(cs.cache, idxs[1])
        # Fallback for things not in cache
        return cs.sol(t; idxs)
    end
    if t != cs.last_t     
        values = cs.sol(t, idxs=cs.vars)
        cs.cache = Dict(cs.vars .=> values)
        cs.last_t = t
    end
    if idxs isa Real
        return cs.cache[idxs]
    else
        return [cs.cache[i] for i in idxs]
    end
end


function getu(cs::CacheSol, syms)
    t->get_cached(cs::CacheSol, t.t, syms)
end

function ModelingToolkit.parameter_values(cs::CacheSol)
    ModelingToolkit.parameter_values(cs.sol)
    # pars = Multibody.collect_all(parameters(cs.model))
    # cs.prob.ps[pars]
end

function (cs::CacheSol)(t; idxs=nothing)
    if idxs === nothing
        cs.sol(t)
    else
        get_cached(cs, t, idxs)
    end
end

function Base.getproperty(cs::CacheSol, s::Symbol)
    s ∈ fieldnames(typeof(cs)) && return getfield(cs, s)
    return getproperty(getfield(cs, :sol), s)
end

function get_cam()
    cam = GLMakie.current_axis().scene.camera
    (; eyeposition=cam.eyeposition[], lookat = cam.view_direction[])
end

"""
    get_rot_fun(sol, frame)

Return a function of `t` that returns the transpose of the rotation-matrix part of `frame`. The transpose is there to make the rotation matrix `R_W_F`.

See also [`get_rot`](@ref)
"""
function get_rot_fun(sol, frame)
    syms = vec(ori(frame).R')
    getter = getu(sol, syms)
    p = ModelingToolkit.parameter_values(sol)
    function (t)
        iv = sol(t)
        temp = ModelingToolkit.ProblemState(; u = iv, p, t) # bad name for something in SII :(
        reshape(getter(temp), 3, 3)
    end
end

"""
    get_fun(sol, syms)

Return a function of `t` that returns `syms` from the solution.
"""
function get_fun(sol, syms)
    getter = getu(sol, syms)
    p = ModelingToolkit.parameter_values(sol)
    function (t)
        iv = sol(t)
        temp = ModelingToolkit.ProblemState(; u = iv, p, t) # bad name for something in SII :(
        getter(temp)
    end
end


"""
    get_frame_fun(sol, frame)

Return a function of `t` that returns the transformation matrix from frame, where the rotation matrix is transposed. The transpose is there to make the rotation matrix `R_W_F`.

See also [`get_frame`](@ref)
"""
function get_frame_fun(sol, frame)
    R = get_rot_fun(sol, frame)
    tr = get_fun(sol, collect(frame.r_0))
    function (t)
        [R(t) tr(t); 0 0 0 1]
    end
end


function get_color(sys, sol, default, var_name = :color)
    try
        Makie.RGBA(sol(sol.t[1], idxs=collect(getproperty(sys, var_name)))...)
    catch
        if default isa AbstractVector
            Makie.RGBA(default...)
        else
            default
        end
    end
end

function get_shapefile(sys, sol)::String
    try
        sf = sol(sol.t[1], idxs=sys.shapefile)
        # decode(sf)
    catch
        ""
    end
end

function get_shape(sys, sol)::String
    try
        sf = sol(sol.t[1], idxs=sys.shape)
        # decode(sf)
    catch
        ""
    end
end


function default_scene(x,y,z; lookat=Vec3f(0,0,0),up=Vec3f(0,1,0),show_axis=false, size)
    # if string(Makie.current_backend()) == "CairoMakie"
    #     scene = Scene() # https://github.com/MakieOrg/Makie.jl/issues/3763
    #     fig = nothing
    # else
        fig = Figure(; size)
        # scene = LScene(fig[1, 1], scenekw = (lights = [DirectionalLight(RGBf(1, 1, 1), Vec3f(-1, 0, 0))],)).scene # This causes a black background for CairoMakie, issue link above
        scene = LScene(fig[1, 1])#.scene
    # end
    cam3d!(scene, center=false)
    # scene.scene.camera.view[] = [
    #     R [x,y,z]; 0 0 0 1
    # ]
    camc = cameracontrols(scene.scene)
    update_cam!(scene.scene, camc, Vec3f(x, y, z), Vec3f(lookat), Vec3f(up))
    fig.current_axis.x.show_axis[] = show_axis
    scene, fig
end

function default_framerate(filename)
    if parse(Bool, get(ENV, "DOCS_BUILD", "false"))
        return 20
    elseif lowercase(last(splitext(filename))) == ".gif"
        return 25
    else
        return 30
    end
end

function render(model, sol,
    timevec::Union{AbstractVector, Nothing} = nothing;
    filename = "multibody_$(getfield(model, :name)).mp4",
    framerate = default_framerate(filename),
    x = 2,
    y = 0.5,
    z = 2,
    lookat = Vec3f(0,0,0),
    up = Vec3f(0,1,0),
    show_axis = false,
    timescale = 1.0,
    traces = nothing,
    display = false,
    loop = 1,
    cache = true,
    size = (600,450),
    kwargs...
    )
    ModelingToolkit.iscomplete(model) || (model = complete(model, flatten=false))
    if sol isa ODEProblem
        sol = FakeSol(model, sol)
        return render(model, sol, 0; x, y, z, lookat, up, show_axis, kwargs...)[1]
    elseif cache
        sol = CacheSol(model, sol)
    end
    scene, fig = default_scene(x,y,z; lookat,up,show_axis, size)
    if timevec === nothing
        timevec = range(sol.t[1], sol.t[end]*timescale, step=1/framerate)
    end

    t = Observable(timevec[1])

    recursive_render!(scene, model, sol, t)

    if traces !== nothing
        tvec = range(sol.t[1], stop=sol.t[end], length=500)
        for frame in traces
            md = get_metadata(frame)
            (md !== nothing) || error("Only frames can be traced in animations.")
            if get(md, Multibody.ModelingToolkit.IsFrame, false)
                points = get_trans(sol, frame, tvec) |> Matrix
            elseif get(md, Multibody.PlanarMechanics.IsFrame2D, false)
                points = get_trans_2d(sol, frame, tvec) |> Matrix
            else
                error("Got fishy frame metadata", md)
            end
            Makie.lines!(scene, points)
        end
    end
    if loop > 1
        timevec = repeat(timevec, loop)
    end
    if display
        Base.display(fig)
        sleep(3)
        fnt = @async begin
            record(fig, filename, timevec; framerate) do time
                if time == timevec[1]
                    Base.display(fig)
                end
                t[] = time/timescale
                sleep(max(0, 1/framerate))
            end
        end
        fn = fetch(fnt)
    else
        fn = record(fig, filename, timevec; framerate) do time
            t[] = time/timescale
        end
    end

    fn, scene, fig
end

function render(model, sol, time::Real;
    traces = nothing,
    x = 2,
    y = 0.5,
    z = 2,
    cache = true,
    size = (1200,1000),
    slider = !(sol isa Union{ODEProblem, FakeSol}),
    kwargs...,
    )

    ModelingToolkit.iscomplete(model) || (model = complete(model, flatten=false))
    # model = ModelingToolkit.toggle_namespacing(model, false)

    if sol isa ODEProblem
        sol = FakeSol(model, sol)
    end
    if cache
        sol = CacheSol(model, sol)
    end

    # fig = Figure()
    # scene = LScene(fig[1, 1]).scene
    # cam3d!(scene)
    scene, fig = default_scene(x,y,z; size, kwargs...)
    # mesh!(scene, Rect3f(Vec3f(-5, -3.6, -5), Vec3f(10, 0.1, 10)), color=:gray) # Floor

    steps = range(sol.t[1], sol.t[end], length=3000)

    if slider
        t = Slider(fig[2, 1], range = steps, startvalue = time).value
    else
        t = Observable(time)
    end
    recursive_render!(scene, model, sol, t)

    if traces !== nothing
        tvec = range(sol.t[1], stop=sol.t[end], length=500)
        for frame in traces
            md = get_metadata(frame)
            (md !== nothing) || error("Only frames can be traced in animations.")
            if get(md, Multibody.ModelingToolkit.IsFrame, false)
                points = get_trans(sol, frame, tvec) |> Matrix
            elseif get(md, Multibody.PlanarMechanics.IsFrame2D, false)
                points = get_trans_2d(sol, frame, tvec) |> Matrix
            else
                error("Got fishy frame metadata", md)
            end
            Makie.lines!(scene, points)
        end
    end
    fig, t, scene
end

function Multibody.loop_render(model, sol; timescale = 1.0, framerate = 30, max_loop = 5, kwargs...)
    fig, t = render(model, sol, sol.t[1]; kwargs...)
    sleeptime = 1/framerate
    timevec = range(sol.t[1], sol.t[end]*timescale, step=sleeptime)
    display(fig)
    @async begin
        for i = 1:max_loop
            for ti in timevec
                execution_time = @elapsed t[] = ti
                sleep(max(0, sleeptime - execution_time))
            end
        end
    end
end

"""
Internal function: Recursively render all subsystem components of a multibody system. If a particular component returns `true` from its `render!` method, indicating that the component performaed rendering, the recursion stops.
"""
function recursive_render!(scene, model, sol, t, first=true)
    # Render the model itself
    if first
        did_render = render!(scene, get_systemtype(model), model, sol, t)
        if did_render === true
            return
        end
    end
    for subsys in getfield(model, :systems)
        system_type = get_systemtype(subsys)
        # did_render = render!(scene, system_type, subsys, sol, t)
        subsys_ns = getproperty(model, getfield(subsys, :name))
        did_render = render!(scene, system_type, subsys_ns, sol, t)
        if !something(did_render, false)
            recursive_render!(scene, subsys_ns, sol, t, false)
        end
    end
end

render!(scene, ::Any, args...) = false # Fallback for systems that have no rendering

function render!(scene, ::typeof(Body), sys, sol, t)
    render, radius, length_fraction, cylinder_radius = sol(sol.t[1], idxs=[sys.render, sys.radius, sys.length_fraction, sys.cylinder_radius]) .|> Float32
    render==true || return true # yes, == true
    color = get_color(sys, sol, :purple)
    r_cm = get_fun(sol, collect(sys.r_cm))
    framefun = get_frame_fun(sol, sys.frame_a)
    thing = @lift begin # Sphere
        Ta = framefun($t)
        coords = (Ta*[r_cm($t); 1])[1:3] # TODO: make use of a proper transformation library instead of rolling own?
        point = Point3f(coords)
        Sphere(point, Float32(radius))
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))

    iszero(r_cm(sol.t[1])) && (return false)

    thing2 = @lift begin # Cylinder
        Ta = framefun($t)
               
        r_cmt = r_cm($t) # r_cm is the center of the sphere in frame a
        iszero(r_cmt)
        coords = (Ta*[r_cmt; 1])[1:3]
        point = Point3f(coords) # Sphere center in world coords

        tip = Point3f(Ta[1:3, 4]) # Origin of frame a in world coords

        d = tip - point
        # d = Point3f(Ta[1:3, 1])
        tip = point + length_fraction*d
        extremity = tip
        origin = point

        Makie.GeometryBasics.Cylinder(origin, extremity, cylinder_radius)
    end
    mesh!(scene, thing2; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))
    false
end

function render!(scene, ::typeof(World), sys, sol, t)
    sol(sol.t[1], idxs=sys.render)==true || return true # yes, == true
    radius = 0.01f0
    r_0 = get_fun(sol, collect(sys.frame_b.r_0))

    thing = @lift begin
        O = Point3f(r_0($t)) # Assume world is never moving
        x = O .+ Point3f(1,0,0)
        Makie.GeometryBasics.Cylinder(O, x, radius)
    end
    mesh!(scene, thing, color=:red)

    thing = @lift begin
        O = Point3f(r_0($t))
        y = O .+ Point3f(0,1,0)
        Makie.GeometryBasics.Cylinder(O, y, radius)
    end
    mesh!(scene, thing, color=:green)

    thing = @lift begin
        O = Point3f(r_0($t))
        z = O .+ Point3f(0,0,1)
        Makie.GeometryBasics.Cylinder(O, z, radius)
    end
    mesh!(scene, thing, color=:blue)
    true
end

function render!(scene, ::typeof(Frame), sys, sol, t)
    sol(sol.t[1], idxs=sys.render)==true || return true # yes, == true
    radius = sol(sol.t[1], idxs=sys.radius) |> Float32
    length = sol(sol.t[1], idxs=sys.length) |> Float32
    T = get_frame_fun(sol, sys)

    thing = @lift begin
        Ti = T($t)
        Rx = Ti[:, 1]
        O = Point3f(Ti[1:3, 4]) # Assume world is never moving
        x = O .+ Point3f(length*Rx)
        Makie.GeometryBasics.Cylinder(O, x, radius)
    end
    mesh!(scene, thing, color=:red)

    thing = @lift begin
        Ti = T($t)
        Ry = Ti[:, 2]
        O = Point3f(Ti[1:3, 4]) # Assume world is never moving
        y = O .+ Point3f(length*Ry)
        Makie.GeometryBasics.Cylinder(O, y, radius)
    end
    mesh!(scene, thing, color=:green)

    thing = @lift begin
        Ti = T($t)
        Rz = Ti[:, 3]
        O = Point3f(Ti[1:3, 4]) # Assume world is never moving
        z = O .+ Point3f(length*Rz)
        Makie.GeometryBasics.Cylinder(O, z, radius)
    end
    mesh!(scene, thing, color=:blue)

    true
end

function render!(scene, ::Union{typeof(Revolute), typeof(RevolutePlanarLoopConstraint)}, sys, sol, t)
    sol(sol.t[1], idxs=sys.render)==true || return true # yes, == true
    r_0 = get_fun(sol, collect(sys.frame_a.r_0))
    n = get_fun(sol, collect(sys.n))
    color = get_color(sys, sol, :red)

    rotfun = get_rot_fun(sol, sys.frame_a)
    radius = try
        sol(sol.t[1], idxs=sys.radius)
    catch
        0.05f0
    end |> Float32
    length = try
        sol(sol.t[1], idxs=sys.length)
    catch
        radius
    end |> Float32
    thing = @lift begin
        O = r_0($t)
        n_a = n($t)
        R_w_a = rotfun($t)
        n_w = R_w_a*n_a # Rotate to the world frame
        p1 = Point3f(O + length*n_w)
        p2 = Point3f(O - length*n_w)
        Makie.GeometryBasics.Cylinder(p1, p2, radius)
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))
    true
end

render!(scene, ::typeof(Universal), sys, sol, t) = false # To recurse down to rendering the revolute joints

function render!(scene, ::typeof(Spherical), sys, sol, t)
    sol(sol.t[1], idxs=sys.render)==true || return true # yes, == true
    vars = get_fun(sol, collect(sys.frame_a.r_0))
    color = get_color(sys, sol, :yellow)
    radius = sol(sol.t[1], idxs=sys.radius)
    thing = @lift begin
        coords = vars($t)
        point = Point3f(coords)
        Sphere(point, Float32(radius))
    end
    mesh!(scene, thing; color)
    true
end

render!(scene, ::typeof(FreeMotion), sys, sol, t) = true


function render!(scene, ::typeof(FixedTranslation), sys, sol, t)
    sol(sol.t[1], idxs=sys.render)==true || return true # yes, == true
    r_0a = get_fun(sol, collect(sys.frame_a.r_0))
    r_0b = get_fun(sol, collect(sys.frame_b.r_0))
    color = get_color(sys, sol, :purple)
    thing = @lift begin
        r1 = Point3f(r_0a($t))
        r2 = Point3f(r_0b($t))
        origin = r1#(r1+r2) ./ 2
        extremity = r2#-r1 # Double pendulum is a good test for this
        radius = Float32(sol($t, idxs=sys.radius))
        Makie.GeometryBasics.Cylinder(origin, extremity, radius)
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1), transparency=true)
    true
end

function render!(scene, ::typeof(BodyShape), sys, sol, t)
    color = get_color(sys, sol, :purple)
    shapepath = get_shapefile(sys, sol)
    shape = get_shape(sys, sol)
    Tshape = reshape(sol(sol.t[1], idxs=sys.shape_transform), 4, 4)
    scale = Vec3f(Float32(sol(sol.t[1], idxs=sys.shape_scale))*ones(Float32, 3))
    if isempty(shapepath)
        r_0a = get_fun(sol, collect(sys.frame_a.r_0))
        r_0b = get_fun(sol, collect(sys.frame_b.r_0))
        shape = get_shape(sys, sol)
        if shape == "cylinder"
            radius = Float32(sol(sol.t[1], idxs=sys.radius))
            thing = @lift begin
                r1 = Point3f(r_0a($t))
                r2 = Point3f(r_0b($t))
                origin = r1
                extremity = r2
                Makie.GeometryBasics.Cylinder(origin, extremity, radius)
            end
            mesh!(scene, thing; color, specular = Vec3f(1.5), transparency=true)
        elseif shape == "box"
            Rfun = get_rot_fun(sol, sys.frame_a)

            width = 2*Float32(sol(sol.t[1], idxs=sys.radius))
            height = width
            r = sol(sol.t[1], idxs=collect(sys.r))
            length = norm(r)
            length_dir = normalize(r)

            width_dir = [1,0,0]'length_dir > 0.9 ? [0,1,0] : [1,0,0]
            height_dir = normalize(cross(length_dir, width_dir))
            width_dir = normalize(cross(height_dir, length_dir))

            r_0a = get_fun(sol, collect(sys.frame_a.r_0)) # Origin is translated by r_shape

            R0 = [length_dir width_dir height_dir]
            @assert isapprox(det(R0), 1.0, atol=1e-6)

            origin = Vec3f(0, -width/2, -height/2)
            extent = Vec3f(length, width, height) 
            thing = Makie.Rect3f(origin, extent)
            m = mesh!(scene, thing; color, specular = Vec3f(1.5), transparency=true)
            on(t) do t
                r1 = Point3f(r_0a(t))
                R = Rfun(t)
                q = Rotations.QuatRotation(R*R0).q
                Q = Makie.Quaternionf(q.v1, q.v2, q.v3, q.s)
                Makie.transform!(m, translation=r1, rotation=Q)
            end
        else
            error("Shape $shape not supported")
        end
    else
        T = get_frame_fun(sol, sys.frame_a)
        scale = Vec3f(Float32(sol(sol.t[1], idxs=sys.shape_scale))*ones(Float32, 3))
        Tshape = reshape(sol(sol.t[1], idxs=sys.shape_transform), 4, 4)
        
        @info "Loading shape mesh $shapepath"
        shapemesh = FileIO.load(shapepath)
        m = mesh!(scene, shapemesh; color, specular = Vec3f(1.5))

        @views on(t) do t
            Ta = T(t)*Tshape
            r1 = Point3f(Ta[1:3, 4])
            q = Rotations.QuatRotation(Ta[1:3, 1:3]).q
            Q = Makie.Quaternionf(q.v1, q.v2, q.v3, q.s)
            Makie.transform!(m; translation=r1, rotation=Q, scale)
        end
    end

    # thing = @lift begin
    #     r1 = Point3f(sol($t, idxs=collect(sys.frame_a.r_0)))
    #     r2 = Point3f(sol($t, idxs=collect(sys.frame_b.r_0)))
    #     origin = r1
    #     extremity = r2
    #     radius = 0.02f0
    #     Makie.GeometryBasics.Sphere((r1+r2) ./ 2, 0.1f0)
    # end
    # mesh!(scene, thing, color=:purple)
    true
end


function render!(scene, ::typeof(BodyCylinder), sys, sol, t)
    
    # NOTE: This draws a solid cylinder without the hole in the middle. Cannot figure out how to render a hollow cylinder
    color = get_color(sys, sol, [1, 0.2, 1, 0.9])
    radius = Float32(sol(sol.t[1], idxs=sys.diameter)/2)
    r_0a = get_fun(sol, collect(sys.frame_a.r_0))
    r_0b = get_fun(sol, collect(sys.frame_b.r_0))
    thing = @lift begin
        r1 = Point3f(r_0a($t))
        r2 = Point3f(r_0b($t))
        origin = r1
        extremity = r2
        Makie.GeometryBasics.Cylinder(origin, extremity, radius)
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5))

    true
end

function render!(scene, ::typeof(BodyBox), sys, sol, t)
    
    # NOTE: This draws a solid box without the hole in the middle. Cannot figure out how to render a hollow box
    color = get_color(sys, sol, [1, 0.2, 1, 0.9])
    width, height, length = Float32.(sol(sol.t[1], idxs=[sys.width, sys.height, sys.render_length]))

    length_dir = sol(sol.t[1], idxs=collect(sys.render_length_dir))
    width_dir = sol(sol.t[1], idxs=collect(sys.render_width_dir))
    height_dir = normalize(cross(normalize(length_dir), normalize(width_dir)))
    width_dir = normalize(cross(height_dir, length_dir))

    Rfun = get_rot_fun(sol, sys.frame_a)
    r_0a = get_fun(sol, collect(sys.frame_a.r_0)) # Origin is translated by r_shape
    r_shape = sol(sol.t[1], idxs=collect(sys.render_r_shape))
    r = sol(sol.t[1], idxs=collect(sys.render_r))

    R0 = [length_dir width_dir height_dir]
    # R0 = Multibody.from_nxy(r, width_dir).R'
    @assert isapprox(det(R0), 1.0, atol=1e-5) "Rotation matrix R0 is not a valid rotation matrix, got `R0 = $R0` with determinant `det(R0) = $(det(R0))`"
    # NOTE: The rotation by this R and the translation with r_shape needs to be double checked

    origin = Vec3f(-length/2, -width/2, -height/2) + r_shape
    extent = Vec3f(length, width, height) 
    thing = Makie.Rect3f(origin, extent)
    m = mesh!(scene, thing; color, specular = Vec3f(1.5))
    on(t) do t
        r1 = Point3f(r_0a(t))
        R = Rfun(t)
        q = Rotations.QuatRotation(R*R0).q
        Q = Makie.Quaternionf(q.v1, q.v2, q.v3, q.s)
        Makie.transform!(m, translation=r1, rotation=Q)
    end

    true
end

function render!(scene, ::typeof(UniversalSpherical), sys, sol, t)
    sphere_diameter = Float32(sol(sol.t[1], idxs=sys.sphere_diameter))
    sphere_color = get_color(sys, sol, [1, 0.2, 1, 0.9], :sphere_color)
    rod_width = Float32(sol(sol.t[1], idxs=sys.rod_width))
    rod_height = Float32(sol(sol.t[1], idxs=sys.rod_height))
    rod_color = get_color(sys, sol, [0, 0.1, 1, 0.9], :rod_color)
    cylinder_length = Float32(sol(sol.t[1], idxs=sys.cylinder_length))
    cylinder_diameter = Float32(sol(sol.t[1], idxs=sys.cylinder_diameter))
    cylinder_color = get_color(sys, sol, [1, 0.2, 0, 1], :cylinder_color)

    # NOTE: the rod is not currently drawn as a box and the revolute cylinders are not drawn at all
    r_0a = get_fun(sol, collect(sys.frame_a.r_0))
    r_0b = get_fun(sol, collect(sys.frame_b.r_0))
    thing = @lift begin
        r1 = Point3f(r_0a($t))
        r2 = Point3f(r_0b($t))
        origin = r1
        extremity = r2
        Makie.GeometryBasics.Cylinder(origin, extremity, rod_width/2)
    end
    mesh!(scene, thing; color=rod_color, specular = Vec3f(1.5))

    # render a sphere for the sperical joint at frame_b
    thing = @lift begin
        r2 = Point3f(r_0b($t))
        Sphere(r2, sphere_diameter/2)
    end
    mesh!(scene, thing; color=sphere_color, specular = Vec3f(1.5))
    true
end

function render!(scene, ::Union{typeof(OneDOFRollingWheelJoint), typeof(SlipWheelJoint)}, sys, sol, t)
    
    r_0 = get_fun(sol, collect(sys.frame_a.r_0))
    # framefun = get_frame_fun(sol, sys.frame_a)
    rotfun = get_rot_fun(sol, sys.frame_a)
    color = get_color(sys, sol, :red)

    radius = try
        sol(sol.t[1], idxs=sys.radius)
    catch
        0.05f0
    end |> Float32
    thing = @lift begin
        # radius = sol($t, idxs=sys.radius)
        O = r_0($t)
        # T_w_a = framefun($t)
        R_w_a = rotfun($t)
        n_w = R_w_a[:, 3] # Wheel rotates around z axis
        # n_w = R_w_a*n_a # Rotate to the world frame
        width = radius/10
        p1 = Point3f(O + width*n_w)
        p2 = Point3f(O - width*n_w)
        Makie.normal_mesh(Makie.Tesselation(Makie.GeometryBasics.Cylinder(p1, p2, radius), 64))
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))
    true
end

function render!(scene, ::typeof(Damper), sys, sol, t)
    r_0a = get_fun(sol, collect(sys.frame_a.r_0))
    r_0b = get_fun(sol, collect(sys.frame_b.r_0))
    color = get_color(sys, sol, :gray)
    radius = sol(sol.t[1], idxs=sys.radius) |> Float32
    length_fraction = sol(sol.t[1], idxs=sys.length_fraction) |> Float32
    thing = @lift begin
        r1 = Point3f(r_0a($t))
        r2 = Point3f(r_0b($t))
        origin = r1
        d = r2 - r1
        extremity = d / norm(d) * length_fraction + r1 
        Makie.GeometryBasics.Cylinder(origin, extremity, Float32(radius))
    end
    mesh!(scene, thing; color)
    true
end


function render!(scene, ::Union{typeof(Spring), typeof(SpringDamperParallel)}, sys, sol, t)
    r_0a = get_fun(sol, collect(sys.frame_a.r_0))
    r_0b = get_fun(sol, collect(sys.frame_b.r_0))
    color = get_color(sys, sol, :blue)
    n_wind = sol(sol.t[1], idxs=sys.num_windings)
    radius = sol(sol.t[1], idxs=sys.radius) |> Float32
    N = round(Int, sol(sol.t[1], idxs=sys.N))
    thing = @lift begin
        r1 = Point3f(r_0a($t))
        r2 = Point3f(r_0b($t))
        spring_mesh(r1,r2; n_wind, radius, N)
    end
    plot!(scene, thing; color)
    true
end

function render!(scene, ::typeof(Multibody.WorldForce), sys, sol, t)
    r_0b = get_fun(sol, collect(sys.frame_b.r_0))
    f = get_fun(sol, collect(sys.frame_b.f))
    R = get_rot_fun(sol, sys.frame_b)
    color = get_color(sys, sol, :green)
    radius = sol(sol.t[1], idxs=sys.radius) |> Float32
    scale = -sol(sol.t[1], idxs=sys.scale) |> Float32

    origin = @lift [Point3f(r_0b($t))]
    d = @lift [scale*Vec3f(R($t)*f($t))]

    Makie.arrows!(origin, d, linecolor = color, arrowcolor = color,
    linewidth = radius, arrowsize = Vec3f(1.3*radius, 1.3*radius, 1.4*radius))
    true
end

function render!(scene, ::Function, sys, sol, t, args...) # Fallback for systems that have at least two frames
    try
        sol(sol.t[1], idxs=sys.render)==true || return true # yes, == true
    catch
    end
    frameinds = findall(ModelingToolkit.isframe, collect(getfield(sys, :systems)))
    length(frameinds) == 2 || return false

    nameof(sys.systems[frameinds[1]]) ∈ (:frame_a, :frame_b) || return false
    nameof(sys.systems[frameinds[2]]) ∈ (:frame_a, :frame_b) || return false

    r_0a = get_fun(sol, collect(sys.frame_a.r_0))
    r_0b = get_fun(sol, collect(sys.frame_b.r_0))
    color = get_color(sys, sol, :green)
    radius = try
        sol(sol.t[1], idxs=sys.radius) |> Float32
    catch
        0.05f0
    end |> Float32
    try
        thing = @lift begin
            r1 = Point3f(r_0a($t))
            r2 = Point3f(r_0b($t))
            origin = r1
            extremity = r2
            Makie.GeometryBasics.Cylinder(origin, extremity, radius)
        end
        mesh!(scene, thing; color)
        return true
    catch
        return false
    end
    false
end

function spring_mesh(p1, p2; n_wind=6, radius=0.1f0, N=200)
    phis = range(0, n_wind*2π, length=N)
    d = p2 - p1
    z = range(0, norm(d), length=N) # Correct length
    dn = d ./ norm(d)
    R = rot_from_line(dn)

    points = map(enumerate(phis)) do (i,phi)
        x = radius*cos(phi)
        y = radius*sin(phi)
        pᵢ = Point3f(x, y, z[i])

        R * pᵢ + p1
    end

    Makie.GeometryBasics.LineString(points)
end

function rot_from_line(d)
    if d[1] == 0 && d[2] == 0
        return RotMatrix{3}(Matrix{Float32}(I, 3, 3))
    end
    d = d ./ norm(d)
    z = SA[0, 0, 1]
    x = cross(z, d)
    x = x ./ norm(x)
    y = cross(d, x)
    y = y ./ norm(y)
    RotMatrix{3}([x y d])
end

Multibody.render!(scene, ::typeof(Multibody.URDFRevolute), sys, sol, t) = false
Multibody.render!(scene, ::typeof(Multibody.URDFPrismatic), sys, sol, t) = false
Multibody.render!(scene, ::typeof(Multibody.NullJoint), sys, sol, t) = false

# ==============================================================================
## PlanarMechanics
# ==============================================================================
function get_rot_fun_2d(sol, frame)
    phifun = get_fun(sol, frame.phi)
    function (t)
        phi = phifun(t)
        [cos(phi) -sin(phi); sin(phi) cos(phi)]
    end
end

function get_frame_fun_2d(sol, frame)
    R = get_rot_fun_2d(sol, frame)
    tr = get_fun(sol, [frame.x, frame.y])
    function (t)
        [R(t) tr(t); 0 0 1]
    end
end

get_trans_2d(sol, frame, t) = SVector{2}(sol(t, idxs = [frame.x, frame.y]))
get_trans_2d(sol, frame, t::AbstractArray) = sol(t, idxs = [frame.x, frame.y])

function render!(scene, ::typeof(P.Frame), sys, sol, t)
    sol(sol.t[1], idxs=sys.render)==true || return true # yes, == true
    radius = sol(sol.t[1], idxs=sys.radius) |> Float32
    length = sol(sol.t[1], idxs=sys.length) |> Float32
    T = get_frame_fun_2d(sol, sys)

    thing = @lift begin
        Ti = T($t)
        Rx = Ti[:, 1]
        O = Point3f([Ti[1:2, 3]; 0]) # Assume world is never moving
        x = O .+ Point3f([length*Rx; 0])
        Makie.GeometryBasics.Cylinder(O, x, radius)
    end
    mesh!(scene, thing, color=:red)

    thing = @lift begin
        Ti = T($t)
        Ry = Ti[:, 2]
        O = Point3f([Ti[1:2, 3]; 0]) # Assume world is never moving
        y = O .+ Point3f([length*Ry; 0])
        Makie.GeometryBasics.Cylinder(O, y, radius)
    end
    mesh!(scene, thing, color=:green)


    true
end

function render!(scene, ::typeof(P.Body), sys, sol, t)
    sol(sol.t[1], idxs=sys.render)==true || return true # yes, == true
    color = get_color(sys, sol, :purple)
    r_cm = get_fun(sol, collect(sys.r))
    framefun = get_frame_fun_2d(sol, sys.frame_a)
    radius = sol(sol.t[1], idxs=sys.radius) |> Float32
    thing = @lift begin # Sphere
        # Ta = framefun($t)
        # coords = (Ta*[0; 0; 1])[1:2]  # r_cm($t)

        coords = r_cm($t)
        point = Point3f(coords..., 0)
        Sphere(point, Float32(radius))
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))
    false
end


function render!(scene, ::Union{typeof(P.FixedTranslation), typeof(P.BodyShape)}, sys, sol, t)
    sol(sol.t[1], idxs=sys.render)==true || return true # yes, == true
    r_0a = get_fun(sol, [sys.frame_a.x, sys.frame_a.y])
    r_0b = get_fun(sol, [sys.frame_b.x, sys.frame_b.y])
    color = get_color(sys, sol, :purple)
    thing = @lift begin
        r1 = Point3f(r_0a($t)..., 0)
        r2 = Point3f(r_0b($t)..., 0)
        origin = r1#(r1+r2) ./ 2
        extremity = r2#-r1 # Double pendulum is a good test for this
        radius = Float32(sol($t, idxs=sys.radius))
        Makie.GeometryBasics.Cylinder(origin, extremity, radius)
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))
    false
end

function render!(scene, ::typeof(P.Prismatic), sys, sol, t)
    sol(sol.t[1], idxs=sys.render)==true || return true # yes, == true
    r_0a = get_fun(sol, [sys.frame_a.x, sys.frame_a.y])
    r_0b = get_fun(sol, [sys.frame_b.x, sys.frame_b.y])
    color = get_color(sys, sol, :purple)
    thing = @lift begin
        r1 = Point3f(r_0a($t)..., 0)
        r2 = Point3f(r_0b($t)..., 0)
        origin = r1#(r1+r2) ./ 2
        extremity = r2#-r1 # Double pendulum is a good test for this
        radius = Float32(sol($t, idxs=sys.radius))
        Makie.GeometryBasics.Cylinder(origin, extremity, radius)
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))
    true
end

function render!(scene, ::typeof(P.Revolute), sys, sol, t)
    sol(sol.t[1], idxs=sys.render)==true || return true # yes, == true
    r_0 = get_fun(sol, [sys.frame_a.x, sys.frame_a.y])
    n = [0,0,1]
    color = get_color(sys, sol, :red)

    rotfun = get_rot_fun_2d(sol, sys.frame_a)
    radius = try
        sol(sol.t[1], idxs=sys.radius)
    catch
        0.05f0
    end |> Float32
    length = try
        sol(sol.t[1], idxs=sys.length)
    catch
        radius
    end |> Float32
    thing = @lift begin
        O = [r_0($t)..., 0]
        R_w_a = cat(rotfun($t), 1, dims=(1,2))
        n_w = R_w_a*n # Rotate to the world frame
        p1 = Point3f(O + length*n_w)
        p2 = Point3f(O - length*n_w)
        Makie.GeometryBasics.Cylinder(p1, p2, radius)
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))
    true
end

function render!(scene, ::Union{typeof(P.Spring), typeof(P.SpringDamper)}, sys, sol, t)
    sol(sol.t[1], idxs=sys.render)==true || return true # yes, == true
    r_0a = get_fun(sol, [sys.frame_a.x, sys.frame_a.y])
    r_0b = get_fun(sol, [sys.frame_b.x, sys.frame_b.y])
    color = get_color(sys, sol, :blue)
    n_wind = sol(sol.t[1], idxs=sys.num_windings)
    radius = sol(sol.t[1], idxs=sys.radius) |> Float32
    N = sol(sol.t[1], idxs=sys.N) |> Int
    thing = @lift begin
        r1 = Point3f(r_0a($t)..., 0)
        r2 = Point3f(r_0b($t)..., 0)
        spring_mesh(r1,r2; n_wind, radius, N)
    end
    plot!(scene, thing; color)
    true
end

function render!(scene, ::Union{typeof(P.SimpleWheel), typeof(P.SlipBasedWheelJoint)}, sys, sol, t)
    
    r_0 = get_fun(sol, [sys.frame_a.x, sys.frame_a.y])
    rotfun = get_rot_fun_2d(sol, sys.frame_a)
    color = get_color(sys, sol, :red)

    # TODO: add some form of assumetry to indicate that the wheel is rotating

    radius = try
        sol(sol.t[1], idxs=sys.radius)
    catch
        0.05f0
    end |> Float32
    z = try
        sol(sol.t[1], idxs=sys.z)
    catch
        0.0
    end |> Float32
    r = try
        sol(sol.t[1], idxs=collect(sys.r))
    catch
        [1, 0]
    end .|> Float32
    n_a = perp(normalize(r)) # Rotation axis
    thing = @lift begin
        O = [r_0($t)..., z]
        R_w_a = rotfun($t)
        n_w = [R_w_a*n_a; 0] # Rotate to the world frame
        width = radius/10
        p1 = Point3f(O + width*n_w)
        p2 = Point3f(O - width*n_w)
        Makie.GeometryBasics.Cylinder(p1, p2, radius)
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))
    true
end

function render!(scene, ::Union{typeof(P.OneDOFSlippingWheelJoint), typeof(P.OneDOFRollingWheelJoint)}, sys, sol, t)
    
    r_0 = get_fun(sol, [sys.frame_a.x, sys.frame_a.y])
    rotfun = get_rot_fun_2d(sol, sys.frame_a)
    color = get_color(sys, sol, :red)

    # TODO: add some form of assumetry to indicate that the wheel is rotating

    radius = try
        sol(sol.t[1], idxs=sys.radius)
    catch
        0.05f0
    end |> Float32
    z = try
        sol(sol.t[1], idxs=sys.z)
    catch
        0.0
    end |> Float32
    r = try
        sol(sol.t[1], idxs=collect(sys.r))
    catch
        [1, 0]
    end .|> Float32
    n_a = [0, 1.0] # Rotation axis
    thing = @lift begin
        O = [r_0($t)..., z]
        n_w = [0, 0, 1.0] # Rotate to the world frame
        width = radius/10
        p1 = Point3f(O + width*n_w)
        p2 = Point3f(O - width*n_w)
        Makie.GeometryBasics.Cylinder(p1, p2, radius)
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))
    false
end

function perp(r)
    if r[1] == 0
        return [1, 0]
    else
        return [-r[2], r[1]]
    end
end


# ==============================================================================
## Visualizers
# ==============================================================================

function render!(scene, ::typeof(Multibody.BoxVisualizer), sys, sol, t)
    
    # NOTE: This draws a solid box without the hole in the middle. Cannot figure out how to render a hollow box
    color = get_color(sys, sol, [1, 0.2, 1, 0.9])
    width = Float32(sol(sol.t[1], idxs=sys.width))
    height = Float32(sol(sol.t[1], idxs=sys.height))
    length = Float32(sol(sol.t[1], idxs=sys.length))

    length_dir = sol(sol.t[1], idxs=collect(sys.render_length_dir))
    width_dir = sol(sol.t[1], idxs=collect(sys.render_width_dir))
    height_dir = normalize(cross(normalize(length_dir), normalize(width_dir)))
    width_dir = normalize(cross(height_dir, length_dir))

    Rfun = get_rot_fun(sol, sys.frame_a)
    r_0a = get_fun(sol, collect(sys.frame_a.r_0)) # Origin is translated by r_shape
    r_shape = sol(sol.t[1], idxs=collect(sys.r_shape))
    # r = sol(sol.t[1], idxs=collect(sys.r))

    R0 = [length_dir width_dir height_dir]
    # R0 = Multibody.from_nxy(r, width_dir).R'
    @assert isapprox(det(R0), 1.0, atol=1e-5) "Rotation matrix R0 is not a valid rotation matrix, got `R0 = $R0` with determinant `det(R0) = $(det(R0))`"
    # NOTE: The rotation by this R and the translation with r_shape needs to be double checked

    origin = Vec3f(-length/2, -width/2, -height/2) + r_shape
    extent = Vec3f(length, width, height) 
    thing = Makie.Rect3f(origin, extent)
    m = mesh!(scene, thing; color, specular = Vec3f(1.5))
    on(t) do t
        r1 = Point3f(r_0a(t))
        R = Rfun(t)
        q = Rotations.QuatRotation(R*R0).q
        Q = Makie.Quaternionf(q.v1, q.v2, q.v3, q.s)
        Makie.transform!(m, translation=r1, rotation=Q)
    end

    true
end

function render!(scene, ::typeof(Multibody.SphereVisualizer), sys, sol, t)
    sol(sol.t[1], idxs=sys.render)==true || return true # yes, == true
    color = get_color(sys, sol, :purple)
    framefun = get_frame_fun(sol, sys.frame_a)
    radius = sol(sol.t[1], idxs=sys.radius) |> Float32
    thing = @lift begin # Sphere
        Ta = framefun($t)
        coords = Ta[1:3, 4]
        point = Point3f(coords)
        Sphere(point, Float32(radius))
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))
end

function render!(scene, ::typeof(Multibody.CylinderVisualizer), sys, sol, t)
    color = get_color(sys, sol, :purple)
    radius = Float32(sol(sol.t[1], idxs=sys.radius))
    r_0a = get_fun(sol, collect(sys.frame_a.r_0))




    length_dir = sol(sol.t[1], idxs=collect(sys.render_length_dir))
    width_dir = randn(3,3)
    height_dir = normalize(cross(normalize(length_dir), normalize(width_dir)))
    width_dir = normalize(cross(height_dir, length_dir))

    Rfun = get_rot_fun(sol, sys.frame_a)

    R0 = [length_dir width_dir height_dir]

    r1 = Point3f(0,0,0)
    r2 = Point3f((length*length_direction)...)
    origin = r1
    extremity = r2
    thing = Makie.GeometryBasics.Cylinder(origin, extremity, radius)
    m = mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1), transparency=true)


    on(t) do t
        r1 = Point3f(r_0a(t))
        R = Rfun(t)
        q = Rotations.QuatRotation(R*R0).q
        Q = Makie.Quaternionf(q.v1, q.v2, q.v3, q.s)
        Makie.transform!(m, translation=r1, rotation=Q)
    end

    true
end

end
