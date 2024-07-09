module Render
using Makie
using Multibody
import Multibody: render, render!, loop_render, encode, decode, get_rot, get_trans, get_frame
using Rotations
using LinearAlgebra
using ModelingToolkit
export render, loop_render
using MeshIO, FileIO
using StaticArrays


"""
    get_rot_fun(sol, frame)

Return a function of `t` that returns the transpose of the rotation-matrix part of `frame`. The transpose is there to make the rotation matrix `R_W_F`.

See also [`get_rot`](@ref)
"""
function get_rot_fun(sol, frame)
    syms = vec(ori(frame).R.mat')
    getter = ModelingToolkit.getu(sol, syms)
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
    getter = ModelingToolkit.getu(sol, syms)
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
        if default isa AbstractVector
            Makie.RGBA(default...)
        else
            default
        end
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
    filename = "multibody_$(model.name).mp4",
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
    kwargs...
    )
    scene, fig = default_scene(x,y,z; lookat,up,show_axis)
    if timevec === nothing
        timevec = range(sol.t[1], sol.t[end]*timescale, step=1/framerate)
    end

    t = Observable(timevec[1])

    recursive_render!(scene, complete(model), sol, t)

    if traces !== nothing
        tvec = range(sol.t[1], stop=sol.t[end], length=500)
        for frame in traces
            (frame.metadata !== nothing && get(frame.metadata, :frame, false)) || error("Only frames can be traced in animations.")
            points = get_trans(sol, frame, tvec) |> Matrix
            Makie.lines!(scene, points)
        end
    end
    if loop > 1
        timevec = repeat(timevec, loop)
    end
    if display
        Base.display(fig)
        sleep(2)
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
    kwargs...,
    )

    # fig = Figure()
    # scene = LScene(fig[1, 1]).scene
    # cam3d!(scene)
    scene, fig = default_scene(x,y,z; kwargs...)
    # mesh!(scene, Rect3f(Vec3f(-5, -3.6, -5), Vec3f(10, 0.1, 10)), color=:gray) # Floor

    steps = range(sol.t[1], sol.t[end], length=3000)

    t = Slider(fig[2, 1], range = steps, startvalue = time).value
    recursive_render!(scene, complete(model), sol, t)

    if traces !== nothing
        tvec = range(sol.t[1], stop=sol.t[end], length=500)
        for frame in traces
            (frame.metadata !== nothing && get(frame.metadata, :frame, false)) || error("Only frames can be traced in animations.")
            points = get_trans(sol, frame, tvec) |> Matrix
            Makie.lines!(scene, points)
        end
    end
    fig, t
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
    radius = sol(sol.t[1], idxs=sys.radius) |> Float32
    cylinder_radius = sol(sol.t[1], idxs=sys.cylinder_radius) |> Float32
    thing = @lift begin # Sphere
        Ta = framefun($t)
        coords = (Ta*[r_cm($t); 1])[1:3] # TODO: make use of a proper transformation library instead of rolling own?
        point = Point3f(coords)
        Sphere(point, Float32(radius))
    end
    mesh!(scene, thing; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))

    iszero(r_cm(sol.t[1])) && (return true)

    thing2 = @lift begin # Cylinder
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

        Makie.GeometryBasics.Cylinder(origin, extremity, cylinder_radius)
    end
    mesh!(scene, thing2; color, specular = Vec3f(1.5), shininess=20f0, diffuse=Vec3f(1))
    true
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

function render!(scene, T::Union{typeof(Revolute), typeof(RevolutePlanarLoopConstraint)}, sys, sol, t)
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

        @views on(t) do t
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
    width = Float32(sol(sol.t[1], idxs=sys.width))
    height = Float32(sol(sol.t[1], idxs=sys.height))
    length = Float32(sol(sol.t[1], idxs=sys.render_length))

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
    @assert isapprox(det(R0), 1.0, atol=1e-6)
    # NOTE: The rotation by this R and the translation with r_shape needs to be double checked

    origin = Vec3f(0, -width/2, -height/2) + r_shape
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
    n_wind = sol(sol.t[1], idxs=sys.num_windings)
    radius = sol(sol.t[1], idxs=sys.radius) |> Float32
    N = sol(sol.t[1], idxs=sys.N) |> Int
    thing = @lift begin
        r1 = Point3f(r_0a($t))
        r2 = Point3f(r_0b($t))
        spring_mesh(r1,r2; n_wind, radius, N)
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
end