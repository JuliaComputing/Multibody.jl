# 3D rendering and animations

Multibody.jl has an automatic 3D-rendering feature that draws a mechanism in 3D. This can be used to create animations of the mechanism's motion from a solution trajectory, as well as to create interactive applications where the evolution of time can be controlled by the user.

The functionality requires the user to load any of the Makie frontend packages, e.g., 
```julia
using GLMakie
```
or 
```julia
using CairoMakie
```

After that, the [`render`](@ref) function is the main entry point to create 3D renderings. This function has the following methods:

- `render(model, solution)`: this method creates an animation corresponding to the mechanisms evolution in a simulation trajectory.
- `scene, time = render(model, solution, t::Real)`: this method opens an interactive window with the mechanism in the configuration corresponding to the time `t`. Display `scene` to display the interactive window, and change the time by either dragging the slider in the window, or write to the observable `time[] = new_time`.

## Rendering API

```@docs
render
render!
```