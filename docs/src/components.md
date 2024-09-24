# Components

The perhaps most fundamental component is a [`Body`](@ref), this component has a single flange, `frame_a`, which is used to connect the body to other components. This component has a mass, a vector `r_cm` from `frame_a` to the center of mass, and a moment of inertia tensor `I` in the center of mass. The body can be thought of as a point mass with a moment of inertia tensor.

A mass with a shape can be modeled using a [`BodyShape`](@ref). The primary difference between a [`Body`](@ref) and a [`BodyShape`](@ref) is that the latter has an additional flange, `frame_b`, which is used to connect the body to other components. The translation between `flange_a` and `flange_b` is determined by the vector `r`. The [`BodyShape`](@ref) is suitable to model, e.g., cylinders, rods, and boxes.

A rod without a mass (just a translation), is modeled using [`FixedTranslation`](@ref).

```@index
```

```@autodocs
Modules = [Multibody, Multibody.PlanarMechanics]
Pages   = ["components.jl", "wheels.jl", "PlanarMechanics/components.jl"]
```
