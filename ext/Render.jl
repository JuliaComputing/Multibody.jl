module Render
using Makie
using Multibody
import Multibody: render, render!, encode, decode
using Rotations
using LinearAlgebra
using ModelingToolkit
export render
using MeshIO, FileIO


function get_rot(sol, frame, t)
    reshape(sol(t, idxs = vec(ori(frame).R.mat)), 3, 3)
end

function get_rot_fun(sol, frame)
    syms = vec(ori(frame).R.mat)
    getter = ModelingToolkit.getu(sol, syms)
    p = ModelingToolkit.parameter_values(sol)
    function (t)
        iv = sol(t)
        temp = ModelingToolkit.ProblemState(; u = iv, p, t) # bad name for something in SII :(
        reshape(getter(temp), 3, 3)
    end
end

function get_fun(sol, syms)
    getter = ModelingToolkit.getu(sol, syms)
    p = ModelingToolkit.parameter_values(sol)
    function (t)
        iv = sol(t)
        temp = ModelingToolkit.ProblemState(; u = iv, p, t) # bad name for something in SII :(
        getter(temp)
    end
end


"""
    get_frame(sol, frame, t)

Extract a 4×4 transformation matrix ∈ SE(3) from a solution at time `t`.
"""
function get_frame(sol, frame, t)
    R = get_rot(sol, frame, t)
    tr = sol(t, idxs = collect(frame.r_0))
    [R' tr; 0 0 0 1]
end

function get_frame_fun(sol, frame)
    R = get_rot_fun(sol, frame)
    tr = get_fun(sol, collect(frame.r_0))
    function (t)
        [R(t)' tr(t); 0 0 0 1]
    end
end


"get_systemtype(sys): Get the constructor of a component for dispatch purposes. This only supports components that have the `gui_metadata` property set. If no metadata is available, nothing is returned."
function get_systemtype(sys)
    meta = getfield(sys, :gui_metadata)
    meta === nothing && return nothing
    eval(meta.type)
end

function get_color(sys, sol, default)
    try
        Makie.RGBA(sol(sol.t[1], idxs=collect(sys.color))...)
    catch
        default
    end
end

function get_shape(sys, sol)::String
    try
        sf = sol(sol.t[1], idxs=collect(sys.shapefile))
        decode(sf)
    catch
        ""
    end
end


function default_scene(x,y,z; lookat=Vec3f(0,0,0),up=Vec3f(0,1,0),show_axis=false)
    # if string(Makie.current_backend()) == "CairoMakie"
    #     scene = Scene() # https://github.com/MakieOrg/Makie.jl/issues/3763
    #     fig = nothing
    # else
        fig = Figure()
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


function render(model, sol,
    timevec::Union{AbstractVector, Nothing} = nothing;
    framerate = 30,
    x = 2,
    y = 0.5,
    z = 2,
    lookat = Vec3f(0,0,0),
    up = Vec3f(0,1,0),
    show_axis = false,
    timescale = 1.0,
    filename = "multibody_$(model.name).mp4",
    kwargs...
    )
    scene, fig = default_scene(x,y,z; lookat,up,show_axis)
    if timevec === nothing
        timevec = range(sol.t[1], sol.t[end]*timescale, step=1/framerate)
    end

    t = Observable(timevec[1])

    recursive_render!(scene, complete(model), sol, t)
    fn = record(fig, filename, timevec; framerate) do time
        t[] = time/timescale
    end
    fn, scene, fig
end

function render(model, sol, time::Real;
    kwargs...,
    )

    # fig = Figure()
    # scene = LScene(fig[1, 1]).scene
    # cam3d!(scene)
    scene, fig = default_scene(0,0,10; kwargs...)
    # mesh!(scene, Rect3f(Vec3f(-5, -3.6, -5), Vec3f(10, 0.1, 10)), color=:gray) # Floor

    steps = range(sol.t[1], sol.t[end], length=3000)

    t = Slider(fig[2, 1], range = steps, startvalue = time).value
    
    recursive_render!(scene, complete(model), sol, t)
    fig, t
end

"""
Internal function: Recursively render all subsystem components of a multibody system. If a particular component returns `true` from its `render!` method, indicating that the component performaed rendering, the recursion stops.
"""
function recursive_render!(scene, model, sol, t)
    for subsys in model.systems
        system_type = get_systemtype(subsys)
        # did_render = render!(scene, system_type, subsys, sol, t)
        subsys_ns = getproperty(model, subsys.name)
        did_render = render!(scene, system_type, subsys_ns, sol, t)
        if !something(did_render, false)
            recursive_render!(scene, subsys_ns, sol, t)
        end
    end
end

render!(scene, ::Any, args...) = false # Fallback for systems that have no rendering

function render!(scene, ::typeof(Body), sys, sol, t)
    color = get_color(sys, sol, :purple)
    r_cm = get_fun(sol, collect(sys.r_cm))
    framefun = get_frame_fun(sol, sys.frame_a)
    radius = try
        sol(sol.t[1], idxs=sys.radius)
    catch
        0.05f0
    end
    thing = @lift begin # Sphere
        Ta = framefun($t)
        coords = (Ta*[r_cm($t); 1])[1:3] # TODO: make use of a proper transformation library instead of rolling own?
        point = Point3f(coords)
        Sphere(point, Float32(radius))
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))

    iszero(r_cm(sol.t[1])) && (return true)

    thing = @lift begin # Cylinder
        Ta = framefun($t)
               
        r_cmt = r_cm($t) # r_cm is the center of the sphere in frame a
        iszero(r_cmt)
        coords = (Ta*[r_cmt; 1])[1:3]
        point = Point3f(coords) # Sphere center in world coords

        tip = Point3f(Ta[1:3, 4]) # Origin of frame a in world coords

        d = tip - point
        d = d ./ norm(d)
        # d = Point3f(Ta[1:3, 1])
        tip = point + 0.2f0d
        extremity = tip
        origin = point

        radius = 0.05f0
        Makie.GeometryBasics.Cylinder(origin, extremity, radius)
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))
    true
end

function render!(scene, ::typeof(World), sys, sol, t)
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

function render!(scene, ::Union{typeof(Revolute), typeof(RevolutePlanarLoopConstraint)}, sys, sol, t)
    # TODO: change to cylinder
    r_0 = get_fun(sol, collect(sys.frame_a.r_0))
    n = get_fun(sol, collect(sys.n))
    color = get_color(sys, sol, :red)

    rotfun = get_rot_fun(sol, sys.frame_a)
    radius = try
        sol(sol.t[1], idxs=sys.radius)
    catch
        0.05f0
    end |> Float32
    thing = @lift begin
        # radius = sol($t, idxs=sys.radius)
        O = r_0($t)
        n_a = n($t)
        R_w_a = rotfun($t)
        n_w = R_w_a'*n_a # Rotate to the world frame
        p1 = Point3f(O + radius*n_w)
        p2 = Point3f(O - radius*n_w)
        Makie.GeometryBasics.Cylinder(p1, p2, radius)
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))
    true
end

function render!(scene, ::typeof(Spherical), sys, sol, t)
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
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))
    true
end

function render!(scene, ::typeof(BodyShape), sys, sol, t)
    color = get_color(sys, sol, :purple)
    shapepath = get_shape(sys, sol)
    if isempty(shapepath)
        radius = Float32(sol(sol.t[1], idxs=sys.radius))
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
    else
        T = get_frame_fun(sol, sys.frame_a)

        @info "Loading shape mesh $shapepath"
        shapemesh = FileIO.load(shapepath)
        m = mesh!(scene, shapemesh; color, specular = Vec3f(1.5))

        on(t) do t
            Ta = T(t)
            r1 = Point3f(Ta[1:3, 4])
            q = Rotations.QuatRotation(Ta[1:3, 1:3]).q
            Q = Makie.Quaternionf(q.v1, q.v2, q.v3, q.s)
            Makie.transform!(m, translation=r1, rotation=Q)
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


function render!(scene, ::typeof(Damper), sys, sol, t)
    r_0a = get_fun(sol, collect(sys.frame_a.r_0))
    r_0b = get_fun(sol, collect(sys.frame_b.r_0))
    color = get_color(sys, sol, :gray)
    thing = @lift begin
        r1 = Point3f(r_0a($t))
        r2 = Point3f(r_0b($t))
        origin = r1
        d = r2 - r1
        extremity = d / norm(d) * 0.2f0 + r1 
        radius = 0.1f0
        Makie.GeometryBasics.Cylinder(origin, extremity, Float32(radius))
    end
    mesh!(scene, thing; color)
    true
end


function render!(scene, ::typeof(Spring), sys, sol, t)
    r_0a = get_fun(sol, collect(sys.frame_a.r_0))
    r_0b = get_fun(sol, collect(sys.frame_b.r_0))
    color = get_color(sys, sol, :blue)
    thing = @lift begin
        r1 = Point3f(r_0a($t))
        r2 = Point3f(r_0b($t))
        spring_mesh(r1,r2)
    end
    plot!(scene, thing; color)
    true
end


function render!(scene, ::Function, sys, sol, t, args...) # Fallback for systems that have at least two frames
    count(ModelingToolkit.isframe, sys.systems) == 2 || return false
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
    phi = range(0, n_wind*2π, length=N)
    
    d = p2 - p1
    
    x = radius*cos.(phi)
    y = radius*sin.(phi)
    z = range(0, norm(d), length=N) # Correct length
    points = Point3f.(x, y, z)

    d = d ./ norm(d)

    # Rotate
    R = rot_from_line(d)
    points = Ref(R) .* points

    # Translate
    points = points .+ p1


    Makie.GeometryBasics.LineString(points)
end

function rot_from_line(d)
    d = d ./ norm(d)
    if d[1] == 0 && d[2] == 0
        return RotMatrix{3}(Matrix{Float32}(I, 3, 3))
    end
    d = d ./ norm(d)
    z = [0, 0, 1]
    x = cross(z, d)
    x = x ./ norm(x)
    y = cross(d, x)
    y = y ./ norm(y)
    RotMatrix{3}([x y d])
end
end