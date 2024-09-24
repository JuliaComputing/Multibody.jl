# Joints

A joint restricts the number of degrees of freedom (DOF) of a body. For example, a free floating body has 6 DOF, but if it is attached to a [`Revolute`](@ref) joint, the joint restricts all but one rotational degree of freedom (a revolute joint acts like a hinge). Similarily, a [`Prismatic`](@ref) joint restricts all but one translational degree of freedom (a prismatic joint acts like a slider).

A [`Spherical`](@ref) joints restricts all translational degrees of freedom, but allows all rotational degrees of freedom. It thus transmits no torque. A [`Planar`](@ref) joint moves in a plane, i.e., it restricts one translational DOF and two rotational DOF. A [`Universal`](@ref) joint has two rotational DOF.

Some joints offer the option to add 1-dimensional components to them by providing the keyword `axisflange = true`. This allows us to add, e.g., springs, dampers, sensors, and actuators to the joint.

## Docstrings
```@index
```

```@autodocs
Modules = [Multibody, Multibody.PlanarMechanics]
Pages   = ["joints.jl", "fancy_joints.jl", "PlanarMechanics/joints.jl"]
```
