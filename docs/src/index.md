```@meta
CurrentModule = Multibody
```

# Multibody

Documentation for [Multibody](https://github.com/JuliaComputing/Multibody.jl).

```@setup logo
using ModelingToolkit
using Multibody, JuliaSimCompiler
using OrdinaryDiffEq # Contains the ODE solver we will use
using Plots
t = Multibody.t
W(args...; kwargs...) = Multibody.world


JULIASIM_PURPLE = [87,87,219,255.0f0]./255 # RGBA
length_scale = 0.5 # This controls the frequency of the oscillations, smaller scale -> faster oscillations
radius_small = length_scale*0.2
radius_large = length_scale*0.3

@mtkmodel Logo begin
    @components begin
        world = World(render=false)
        revl  = Revolute(; radius = radius_large, color=JULIASIM_PURPLE, axisflange=true)
        revl2 = Revolute(; radius = radius_large, color=JULIASIM_PURPLE, axisflange=true)
        revr  = Revolute(; radius = radius_small, color=JULIASIM_PURPLE, axisflange=true)
        bodyl = Body(m=1, radius = radius_small, color=JULIASIM_PURPLE)
        bodyr = Body(m=1, radius = radius_large, color=JULIASIM_PURPLE)
        bar_top = FixedTranslation(r=length_scale*[1, 0.05, 0], radius=length_scale*0.025, color=JULIASIM_PURPLE)
        barl   = FixedTranslation(r=length_scale*[1, -1, 0], radius=length_scale*0.025, color=JULIASIM_PURPLE)
        barr   = FixedTranslation(r=1.1*length_scale*[-1, -1, 0], radius=length_scale*0.025, color=JULIASIM_PURPLE)

        damperl  = Rotational.Damper(d=0.1)
        damperl2 = Rotational.Damper(d=0.01)
        damperr  = Rotational.Damper(d=0.01)
    end
    @equations begin
        connect(revl.frame_a, world.frame_b)

        connect(revl.frame_b, barl.frame_a)
        connect(barl.frame_b, bodyl.frame_a)

        connect(world.frame_b, revl2.frame_a)
        connect(revl2.frame_b, bar_top.frame_a)
        connect(bar_top.frame_b, revr.frame_a)
        connect(revr.frame_b, barr.frame_a)
        connect(barr.frame_b, bodyr.frame_a)

        connect(revl.axis, damperl.flange_a)
        connect(revl2.axis, damperl2.flange_a)
        connect(revr.axis, damperr.flange_a)

        connect(revl.support, damperl.flange_b)
        connect(revl2.support, damperl2.flange_b)
        connect(revr.support, damperr.flange_b)
    end
end

@named logo = Logo()
logo = complete(logo)
ssys = structural_simplify(IRSystem(logo))
prob = ODEProblem(ssys, [], (0.0, 3.51))
sol = solve(prob, Rodas5P())
Plots.plot(sol)

import GLMakie
framerate = 30
timevec = [zeros(30); range(sol.t[1], sol.t[end], step=1/framerate)] |> reverse
render(logo, sol, timevec; z=-2.5, x=0, y=-0.5, lookat=[0.2,-0.5, 0], filename="JuliaSim_logo.gif", framerate)
```
![animated logo](JuliaSim_logo.gif)


Welcome to the world of Multibody.jl, a powerful and flexible component of JuliaSim designed to model, analyze, and simulate multibody systems in Julia. As a state-of-the-art tool, Multibody.jl enables users to efficiently study the dynamics of complex mechanical systems in various fields, such as robotics, biomechanics, aerospace, and vehicle dynamics.

Built on top of the Julia language and the JuliaSim suite of tools for modeling, simulation, optimization and control, Multibody.jl harnesses the power of Julia's high-performance computing capabilities, making it a go-to choice for both researchers and engineers who require fast simulations and real-time performance. With an intuitive syntax and a comprehensive set of features, this package seamlessly integrates with other Julia and JuliaSim libraries, enabling users to tackle diverse and sophisticated problems in multibody dynamics.

In this documentation, you will find everything you need to get started with Multibody.jl, from basic component descriptions to detailed examples showcasing the package's capabilities. As you explore this documentation, you'll learn how to create complex models, work with forces and torques, simulate various types of motions, and visualize your results in both 2D and 3D. Whether you are a seasoned researcher or a newcomer to the field, Multibody.jl will empower you to bring your ideas to life and unlock new possibilities in the fascinating world of multibody dynamics.

## Example overview
The following animations give a quick overview of simple mechanisms that can be modeled using Multibody.jl. The examples are ordered from simple at the top, to more advanced at the bottom. Please browse the examples for even more examples!
```@raw html
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GIF Grid</title>
    <style>
        .grid-container {
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            grid-template-rows: repeat(4, auto);
            gap: 10px;
            padding: 20px;
        }
        .grid-item {
            width: 100%;
            height: auto;
        }
        .grid-item img {
            width: 100%;
            height: auto;
            display: block;
        }
    </style>
</head>
<body>

<div class="grid-container">
    <a class="grid-item" href="https://help.juliahub.com/multibody/dev/examples/pendulum/">
        <img src="https://help.juliahub.com/multibody/dev/examples/furuta.gif" alt="Furuta">
    </a>
    <a class="grid-item" href="https://help.juliahub.com/multibody/dev/examples/spring_damper_system/">
        <img src="https://help.juliahub.com/multibody/dev/examples/springdamper.gif" alt="springdamper">
    </a>
    <a class="grid-item" href="https://help.juliahub.com/multibody/dev/examples/wheel/">
        <img src="https://help.juliahub.com/multibody/dev/examples/wheelset.gif" alt="wheels">
    </a>
    <a class="grid-item" href="https://help.juliahub.com/multibody/dev/examples/three_springs/">
        <img src="https://help.juliahub.com/multibody/dev/examples/three_springs.gif" alt="three_springs">
    </a>
    <a class="grid-item" href="https://help.juliahub.com/multibody/dev/examples/space/">
        <img src="https://help.juliahub.com/multibody/dev/examples/space.gif" alt="space">
    </a>
    <a class="grid-item" href="https://help.juliahub.com/multibody/dev/examples/free_motion/#Body-suspended-in-springs">
        <img src="https://help.juliahub.com/multibody/dev/examples/free_body.gif" alt="free_body">
    </a>
    <a class="grid-item" href="https://help.juliahub.com/multibody/dev/examples/ropes_and_cables/">
        <img src="https://help.juliahub.com/multibody/dev/examples/flexible_rope.gif" alt="flexible_rope">
    </a>
    <a class="grid-item" href="https://help.juliahub.com/multibody/dev/examples/ropes_and_cables/">
        <img src="https://help.juliahub.com/multibody/dev/examples/mounted_chain.gif" alt="mounted_chain">
    </a>
    <a class="grid-item" href="https://help.juliahub.com/multibody/dev/examples/quad/">
        <img src="https://help.juliahub.com/multibody/dev/examples/quadrotor.gif" alt="quadrotor">
    </a>
    <a class="grid-item" href="https://help.juliahub.com/multibody/dev/examples/kinematic_loops/">
        <img src="https://help.juliahub.com/multibody/dev/examples/fourbar2.gif" alt="fourbar2">
    </a>
    <a class="grid-item" href="https://help.juliahub.com/multibody/dev/examples/kinematic_loops/">
        <img src="https://help.juliahub.com/multibody/dev/examples/fourbar.gif" alt="fourbar">
    </a>
    <a class="grid-item" href="https://help.juliahub.com/multibody/dev/examples/robot/">
        <img src="https://help.juliahub.com/multibody/dev/examples/robot.gif" alt="robot">
    </a>
</div>

</body>
</html>
```


## Installation
To install this library, first follow the [installation instructions for JuliaSimCompiler](https://juliacomputing.github.io/JuliaSimCompiler.jl/stable/#Installing-and-Using-JuliaSimCompiler). In particular, you need to [add the JuliaHub Registry](https://help.juliahub.com/juliasim/dev/gettingstarted/juliahubregistry/). 

After the registry is added and JuliaSimCompiler is installed, you may install this package using
```julia
import Pkg
Pkg.add("Multibody")
```


## Notable differences from Modelica

- The torque variable in Multibody.jl is typically called `tau` rather than `t` to not conflict with the often used independent variable `t` used to denote time.
- Multibody.jl occasionally requires the user to specify which component should act as the root of the kinematic tree. This only occurs when bodies are connected directly to force components without a joint parallel to the force component.
- In Multibody.jl, the orientation object of a [`Frame`](@ref) is accessed using the function [`ori`](@ref).
- Quaternions in Multibody.jl follow the order ``[s, i, j, k]``, i.e., scalar/real part first.

## 2D and 3D modeling
Multibody.jl offers components for modeling in both 2D and 3D. 2D modeling, often referred to as planar mechanics, is a subset of 3D modeling where the motion is constrained to a plane, the x,y plane. Planar mechanics is sometimes referred to as 3 degrees of freedom (DOF) modeling, referring to the 2 translational DOF and one rotational DOF that the plane offers. Most components in Multibody.jl are aimed at 3D modeling (sometimes referred to as 6 DOF modeling), but components for 2D modeling exist in the submodule `Multibody.PlanarMechanics`.

The components from [`ModelingToolkitStandardLibrary.Mechanical`](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/mechanical/) are 1D, i.e., a single degree of freedom only. These components can be used in both 2D and 3D modeling together with Multibody components that have support for attaching 1D components, such as joints supporting the `axisflange` keyword.

