# Working with orientation and rotation

Orientations and rotations in 3D can be represented in multiple different ways. Components which (may) have a 3D angular state, such as [`Body`](@ref), [`Spherical`](@ref) and [`FreeMotion`](@ref), offer the user to select the orientation representation, either Euler angles or quaternions.

## Euler angles
[Euler angles](https://en.wikipedia.org/wiki/Euler_angles) represent orientation using rotations around three axes, and thus uses three numbers to represent the orientation. The benefit of this representation is that it is minimal (only three numbers used), but the drawback is that any 3-number orientation representation suffers from a kinematic singularity. This representation is beneficial when you know that the singularity will not come into play in your simulation.

Most components that may use Euler angles also allow you to select the sequence of axis around which to perform the rotations, e.g., `sequence = [1,2,3]` performs rotations around ``x`` first, then ``y`` and ``z``.

## Quaternions
A [quaternion](https://en.wikipedia.org/wiki/Quaternion) represents an orientation using 4 numbers. This is a non-minimal representation, but in return it is also singularity free. Multibody.jl uses _non-unit quaternions_[^quat] to represent orientation when `quat = true` is provided to components that support this option. The convention used for quaternions is ``[s, v_1, v_2, v_3]``, sometimes also referred to as ``[w, i, j, k]``, i.e., the real/scalar part comes first, followed by the three imaginary numbers. When quaternions are used, state variables `Q̂` denote non-unit quaternions, while normalized unit quaternions are available as observable variables `Q`. The use of non-unit quaternions allows Multibody to integrate rotations without using dynamic state selection or introducing algebraic equations. 

[^quat]: "Integrating Rotations using Non-Unit Quaternions", Caleb Rucker, https://par.nsf.gov/servlets/purl/10097724

Multibody.jl depends on [Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl) which in turn uses [Quaternions.jl](https://github.com/JuliaGeometry/Quaternions.jl) for quaternion computations. If you manually create quaternions using these packages, you may convert them to a vector to provide, e.g., initial conditions, using `Rotations.params(Q)` (see [Conversion between orientation formats](@ref) below).

### Examples using quaternions
- [Free motions](@ref) (second example on the page)
- [Three springs](@ref)
- [Bodies in space](@ref) (may use, see comment)

## Rotation matrices
Rotation matrices represent orientation using a ``3\times 3`` matrix ``\in SO(3)``. These are used in the equations of multibody components and connectors, but should for the most part be invisible to the user. In particular, they should never appear as state variables after simplification. 


## Conversion between orientation formats
You may convert between different representations of orientation using the appropriate constructors from Rotations.jl, for example:
```@example ORI
using Multibody.Rotations
using Multibody.Rotations: params
using Multibody.Rotations.Quaternions
using LinearAlgebra

R = RotMatrix{3}(I(3))
```

```@example ORI
# Convert R to a quaternion
Q = QuatRotation(R)
```

```@example ORI
# Convert Q to a 4-vector
Qvec = params(Q)
```

```@example ORI
# Convert R to Euler angles in the sequence XYZ
E = RotXYZ(R)
```

```@example ORI
# Convert E to a 3-vector
Evec = params(E)
```

## Conventions for modeling
See [Orientations and directions](@ref)


## Orientation API
See [Orientation utilities](@ref)