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
filename = joinpath(dirname(pathof(Multibody)), "..", "test/urdf/doublependulum.urdf")
out = "multibody_urdf.jl"
urdf2multibody(filename; extras=true, out)

include(joinpath(pwd(), out)) # Include model, perform simulation and plotting
```

## Docstring

```@docs
Multibody.urdf2multibody
```

## Limitations
The URDF import currently has the following limitations:
- Sensors are not imported.
- Transmissions are not imported.
- `friction` is not translated, but `damping` is translated to a 1D `Damping` component attached using an `axisflange`.
- Meshes are not fully supported yet, they will be imported as generic shapes (inertial properties are imported).

## Structure of the translated model
URDF does not store the transformation implied by links in the link itself, instead, the links store visual and inertial geometry information, while the translation between frames is implied by the origin of the following joint(s). Therefore, we do generally not make use of the `r` argument to bodies, and let this be arbitrarily set. The transformation between two joints is instead encoded as a `r` and `R` arguments to each joint, where joints are wrapped in `URDFRevolute` and `URDFPrismatic` components respectively. Internally, these wrapper components are comprised of a transformation, [`FixedTranslation`](@ref) or [`FixedRotation`](@ref), followed by the actual joint. The interface to these special joints are identical to their non-wrapped counterparts, i.e., they have the `frame_a` and `frame_b` connectors as expected. Due to this approach, we always connect to the `frame_a` connector of links/bodies and let `frame_b` be unused.