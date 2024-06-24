# 3D rendering and animations

Multibody.jl has an automatic 3D-rendering feature that draws a mechanism in 3D. This can be used to create animations of the mechanism's motion from a solution trajectory, as well as to create interactive applications where the evolution of time can be controlled by the user.

The functionality requires the user to load any of the Makie frontend packages, e.g., 
```julia
using GLMakie # Preferred when running locally
```
or 
```julia
using CairoMakie
```
!!! note "Backend choice"
    GLMakie and WGLMakie produce much nicer-looking animations and is also significantly faster than CairoMakie. CairoMakie may be used to produce the graphics in some web environments if constraints imposed by the web environment do not allow any of the GL alternatives. CairoMakie struggles with the Z-order of drawn objects, sometimes making bodies that should have been visible hidden behind bodies that are further back in the scene.

After that, the [`render`](@ref) function is the main entry point to create 3D renderings. This function has the following methods:

- `render(model, solution)`: this method creates an animation corresponding to the mechanisms evolution in a simulation trajectory.
- `scene, time = render(model, solution, t::Real)`: this method opens an interactive window with the mechanism in the configuration corresponding to the time `t`. Display `scene` to display the interactive window, and change the time by either dragging the slider in the window, or write to the observable `time[] = new_time`.

## Colors
Many components allows the user to select with which color it is rendered. This choice is made by providing a 4-element array with color values in the order (RGBA), where each value is between 0 and 1. The last value is the alpha channel which determines the opacity, i.e., 1 is opaque and 0 is invisible.

## Rendering the world frame
The display of the world frame can be turned off by setting `world.render => false` in the variable map.


## Rendering API

```@docs
render
render!
```