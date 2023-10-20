using GLMakie

"""
    get_frame(sol, frame, t)

Extract a 4×4 transformation matrix ∈ SE(3) from a solution at time `t`.
"""
function get_frame(sol, frame, t)
    R = reshape(sol(t, idxs = vec(ori(frame).R.mat)), 3, 3)
    t = sol(t, idxs = collect(frame.r_0))
    [R t; 0 0 0 1]
end


"""
    render(model, sol; framerate = 30, timevec = range(sol.t[1], sol.t[end], step = 1 / framerate))

Create a 3D animation of a multibody system

# Arguments:
- `model`: The unsimplified multibody model
- `sol`: The `ODESolution` produced by simulating the system using `solve`
- `framerate`: Number of frames per second.
- `timevec`: The times at which to render a frame.
"""
function render(model, sol;
    framerate = 30,
    timevec = range(sol.t[1], sol.t[end], step=1/framerate),
    )
    # with_theme(theme_dark()) do
        scene = Scene()
        
        cam3d!(scene)
        scene.camera.view[] = [
            1 0 0 0
            0 1 0 0
            0 0 1 -10
            0 0 0 1
        ]
        t = 0.0
        t = Observable(timevec[1])

        for subsys in model.systems
            system_type = get_systemtype(subsys)
            render!(scene, system_type, subsys, sol, t)
        end

        record(scene, "multibody_$(model.name).mp4", timevec; framerate) do time
            t[] = time
        end
    # end
end

function render(model, sol, time::Real;
    framerate = 30,
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
    # t = Observable(time)

    for subsys in model.systems
        system_type = get_systemtype(subsys)
        render!(scene, system_type, subsys, sol, t)
    end
    fig, t
end


"""
    render!(scene, ::typeof(ComponentConstructor), sys, sol, t)

Each component that can be rendered must have a `render!` method. This method is called by `render` for each component in the system.

This method is responsible for drawing the component onto the scene the way it's supposed to look at time `t` in the solution `sol`.
`t` is an Observable. It's recommended to follow the pattern
```julia
thing = @lift begin
    acces relevant coordinates from sol at time t
    create a geometric object that can be rendered
end
mesh!(scene, thing; style...)
```
"""
function render!(scene, ::typeof(Body), sys, sol, t)
    thing = @lift begin
        r = sol($t, idxs=collect(sys.r_cm))
        # Ra, ta = get_frame(sol, sys.frame_a, $t)
        # coords = ta + Ra * r

        Ta = get_frame(sol, sys.frame_a, $t)
        coords = (Ta\[r; 1])[1:3] # TODO: make use of a proper transformation library instead of rolling own?
        point = Point3f(coords)
        Sphere(point, 0.1)
    end
    mesh!(scene, thing, color=:purple)
    thing = @lift begin
        r = sol($t, idxs=collect(sys.r_cm))
        # Ra, ta = get_frame(sol, sys.frame_a, $t)
        # coords = ta + Ra * r

        Ta = get_frame(sol, sys.frame_a, $t)
        coords = (Ta\[r; 1])[1:3]
        ta = Ta[1:3, 4]

        point = Point3f(coords)

        origin = Point3f(ta)
        extremity = point
        radius = 0.02f0
        Makie.GeometryBasics.Cylinder(origin, extremity, radius)
    end
    mesh!(scene, thing, color=:purple)
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
end

function render!(scene, ::typeof(Revolute), sys, sol, t)
    thing = @lift begin
        vars = collect(sys.frame_a.r_0)
        coords = sol($t, idxs=vars)
        point = Point3f(coords)
        Sphere(point, 0.1)
    end
    mesh!(scene, thing, color=:red)
end


function render!(scene, ::typeof(FixedTranslation), sys, sol, t)
    thing = @lift begin
        r1 = Point3f(sol($t, idxs=collect(sys.frame_a.r_0)))
        r2 = Point3f(sol($t, idxs=collect(sys.frame_b.r_0)))
        origin = r1#(r1+r2) ./ 2
        extremity = r2
        radius = 0.1f0
        Makie.GeometryBasics.Cylinder(origin, extremity, radius)
    end
    mesh!(scene, thing, color=:purple)
end

function render!(scene, ::typeof(BodyShape), sys, sol, t)
    thing = @lift begin
        r1 = Point3f(sol($t, idxs=collect(sys.frame_a.r_0)))
        r2 = Point3f(sol($t, idxs=collect(sys.frame_b.r_0)))
        origin = r1
        extremity = r2
        radius = 0.08f0
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
end

"get_systemtype(sys): Get the constructor of a component for dispatch purposes. This only supports components that have the `gui_metadata` property set. If no metadata is available, nothing is returned."
function get_systemtype(sys)
    meta = getfield(sys, :gui_metadata)
    meta === nothing && return nothing
    eval(meta.type)
end

render!(scene, ::Any, args...) = () # Fallback for systems that have no rendering

function render!(scene, ::Function, sys, sol, t, args...) # Fallback for systems that have at least two frames
    count(ModelingToolkit.isframe, sys.systems) == 2 || return
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
    catch
    end
end

