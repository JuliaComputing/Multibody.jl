# URDF import

Multibody.jl supports import of [URDF files](https://wiki.ros.org/urdf) by means of the function [`urdf2multibody`](@ref).
The functionality requires the user to install and load the packages LightXML.jl, Graphs.jl, MetaGraphsNext.jl and JuliaFormatter.jl, e.g.,
```julia
using Pkg
Pkg.add([
    "LightXML",
    "Graphs",
    "MetaGraphsNext",
    "JuliaFormatter"
])
using Multibody, LightXML, Graphs, MetaGraphsNext, JuliaFormatter
```


## Usage
The following example demonstrates how to import a URDF file, the translated model is saved in file `multibody_urdf.jl`. `extras = true` makes the file self contained by including package imports, simulation and plotting.
```julia
filename = joinpath(dirname(pathof(Multibody)), "..", "test/doublependulum.urdf")
out = "multibody_urdf.jl"
urdf2multibody(filename; extras=true, out)

include(out) # Include model, perform simulation and plotting
```

```@docs
Multibody.urdf2multibody
```