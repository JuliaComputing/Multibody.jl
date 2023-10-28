module Render
using Makie
using Multibody
import Multibody: render, render!
using Rotations
using LinearAlgebra
using ModelingToolkit
export render

"""
    get_frame(sol, frame, t)

Extract a 4×4 transformation matrix ∈ SE(3) from a solution at time `t`.
"""
function get_frame(sol, frame, t)
    R = reshape(sol(t, idxs = vec(ori(frame).R.mat)), 3, 3)
    t = sol(t, idxs = collect(frame.r_0))
    [R' t; 0 0 0 1]
end


"get_systemtype(sys): Get the constructor of a component for dispatch purposes. This only supports components that have the `gui_metadata` property set. If no metadata is available, nothing is returned."
function get_systemtype(sys)
    meta = getfield(sys, :gui_metadata)
    meta === nothing && return nothing
    eval(meta.type)
end



function render(model, sol,
    timevec::Union{AbstractVector, Nothing} = nothing;
    framerate = 30,
    x = 0,
    y = 0,
    z = -10,
    R = I(3),
    filename = "multibody_$(model.name).mp4",
    )
    if timevec === nothing
        timevec = range(sol.t[1], sol.t[end], step=1/framerate)
    end
    # with_theme(theme_dark()) do
        scene = Scene()
        
        cam3d!(scene)
        scene.camera.view[] = [
            R [x,y,z]; 0 0 0 1
        ]
        t = 0.0
        t = Observable(timevec[1])

        recursive_render!(scene, complete(model), sol, t)

        record(scene, filename, timevec; framerate) do time
            t[] = time
        end
    # end
end

function render(model, sol, time::Real;
    kwargs...,
    )

    fig = Figure()
    scene = LScene(fig[1, 1]).scene

    steps = range(sol.t[1], sol.t[end], length=300)

    t = Slider(fig[2, 1], range = steps, startvalue = time).value
    
    cam3d!(scene)
    scene.camera.view[] = [
        1 0 0 0
        0 1 0 0
        0 0 1 -10
        0 0 0 1
    ]

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
        did_render = render!(scene, system_type, getproperty(model, subsys.name), sol, t)
        if !something(did_render, false)
            recursive_render!(scene, getproperty(model, subsys.name), sol, t)
        end
    end
end

render!(scene, ::Any, args...) = false # Fallback for systems that have no rendering

function render!(scene, ::typeof(Body), sys, sol, t)
    thing = @lift begin # Sphere
        r = sol($t, idxs=collect(sys.r_cm))
        Ta = get_frame(sol, sys.frame_a, $t)
        coords = (Ta*[r; 1])[1:3] # TODO: make use of a proper transformation library instead of rolling own?
        point = Point3f(coords)
        Sphere(point, 0.1)
    end
    mesh!(scene, thing, color=:purple)

    # thing = @lift begin # Cylinder
    #     Ta = get_frame(sol, sys.frame_a, $t)
               
    #     r = sol($t, idxs=collect(sys.r_cm))
    #     coords = (Ta*[r; 1])[1:3]
    #     point = Point3f(coords)

    #     rt = sol($t, idxs=collect(sys.r_0))
    #     coords = (Ta*[rt; 1])[1:3]
    #     tip = Point3f(coords)

    #     d = tip - point
    #     d = d ./ norm(d)
    #     # d = Point3f(Ta[1:3, 1])
    #     tip = point + 0.2f0d
    #     extremity = tip
    #     origin = point

    #     radius = 0.05f0
    #     Makie.GeometryBasics.Cylinder(origin, extremity, radius)
    # end
    # mesh!(scene, thing, color=:purple)
    true
end

function render!(scene, ::typeof(World), sys, sol, t)
    radius = 0.01f0
    thing = @lift begin
        O = Point3f(sol($t, idxs=collect(sys.frame_b.r_0)))
        x = O .+ Point3f(1,0,0)
        Makie.GeometryBasics.Cylinder(O, x, radius)
    end
    mesh!(scene, thing, color=:red)

    thing = @lift begin
        O = Point3f(sol($t, idxs=collect(sys.frame_b.r_0)))
        y = O .+ Point3f(0,1,0)
        Makie.GeometryBasics.Cylinder(O, y, radius)
    end
    mesh!(scene, thing, color=:green)

    thing = @lift begin
        O = Point3f(sol($t, idxs=collect(sys.frame_b.r_0)))
        z = O .+ Point3f(0,0,1)
        Makie.GeometryBasics.Cylinder(O, z, radius)
    end
    mesh!(scene, thing, color=:blue)
    true
end

function render!(scene, ::typeof(Revolute), sys, sol, t)
    # TODO: change to cylinder
    thing = @lift begin
        # radius = sol($t, idxs=sys.radius)
        O = sol($t, idxs=collect(sys.frame_a.r_0))
        n_a = sol($t, idxs=collect(sys.n))
        R_w_a = get_frame(sol, sys.frame_a, $t)[1:3, 1:3]
        n_w = R_w_a*n_a # Rotate to the world frame
        p1 = Point3f(O + 0.1*n_w)
        p2 = Point3f(O - 0.1*n_w)
        Makie.GeometryBasics.Cylinder(p1, p2, 0.1f0)
    end
    mesh!(scene, thing, color=:red)
    true
end

function render!(scene, ::typeof(Spherical), sys, sol, t)
    thing = @lift begin
        vars = collect(sys.frame_a.r_0)
        coords = sol($t, idxs=vars)
        point = Point3f(coords)
        Sphere(point, 0.1)
    end
    mesh!(scene, thing, color=:yellow)
    true
end

render!(scene, ::typeof(FreeMotion), sys, sol, t) = true


function render!(scene, ::typeof(FixedTranslation), sys, sol, t)
    thing = @lift begin
        r1 = Point3f(sol($t, idxs=collect(sys.frame_a.r_0)))
        r2 = Point3f(sol($t, idxs=collect(sys.frame_b.r_0)))
        origin = r1#(r1+r2) ./ 2
        extremity = r2#-r1 # Double pendulum is a good test for this
        radius = 0.08f0
        Makie.GeometryBasics.Cylinder(origin, extremity, radius)
    end
    mesh!(scene, thing, color=:purple)
    true
end

function render!(scene, ::typeof(BodyShape), sys, sol, t)
    thing = @lift begin
        r1 = Point3f(sol($t, idxs=collect(sys.frame_a.r_0)))
        r2 = Point3f(sol($t, idxs=collect(sys.frame_b.r_0)))
        origin = r1
        extremity = r2
        radius = Float32(sol($t, idxs=sys.radius))
        Makie.GeometryBasics.Cylinder(origin, extremity, radius)
    end
    mesh!(scene, thing, color=:purple)
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
    thing = @lift begin
        r1 = Point3f(sol($t, idxs=collect(sys.frame_a.r_0)))
        r2 = Point3f(sol($t, idxs=collect(sys.frame_b.r_0)))
        origin = r1
        d = r2 - r1
        extremity = d / norm(d) * 0.2f0 + r1 
        radius = 0.1f0
        Makie.GeometryBasics.Cylinder(origin, extremity, radius)
    end
    mesh!(scene, thing, color=:gray)
    true
end


function render!(scene, ::typeof(Spring), sys, sol, t)
    thing = @lift begin
        r1 = Point3f(sol($t, idxs=collect(sys.frame_a.r_0)))
        r2 = Point3f(sol($t, idxs=collect(sys.frame_b.r_0)))
        spring_mesh(r1,r2)
    end
    plot!(scene, thing, color=:blue)
    true
end


function render!(scene, ::Function, sys, sol, t, args...) # Fallback for systems that have at least two frames
    count(ModelingToolkit.isframe, sys.systems) == 2 || return false
    try
        thing = @lift begin
            r1 = Point3f(sol($t, idxs=collect(sys.frame_a.r_0)))
            r2 = Point3f(sol($t, idxs=collect(sys.frame_b.r_0)))
            origin = r1
            extremity = r2
            radius = 0.05f0
            Makie.GeometryBasics.Cylinder(origin, extremity, radius)
        end
        mesh!(scene, thing, color=:green, alpha=0.5)
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
        return Matrix{Float32}(I, 3, 3)
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