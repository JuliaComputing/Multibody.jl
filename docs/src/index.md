```@meta
CurrentModule = Multibody
```

# Multibody

Documentation for [Multibody](https://github.com/YingboMa/Multibody.jl).


Welcome to the world of Multibody.jl, a powerful and flexible component of JuliaSim designed to model, analyze, and simulate multibody systems in Julia. As a state-of-the-art tool, Multibody.jl enables users to efficiently study the dynamics of complex mechanical systems in various fields, such as robotics, biomechanics, aerospace, and vehicle dynamics.

Built on top of the Julia language and the JuliaSim suite of tools for modeling, simulation, optimization and control, Multibody.jl harnesses the power of Julia's high-performance computing capabilities, making it a go-to choice for both researchers and engineers who require fast simulations and real-time performance. With an intuitive syntax and a comprehensive set of features, this package seamlessly integrates with other Julia and JuliaSim libraries, enabling users to tackle diverse and sophisticated problems in multibody dynamics.

In this documentation, you will find everything you need to get started with Multibody.jl, from basic component descriptions to detailed examples showcasing the package's capabilities. As you explore this documentation, you'll learn how to create complex models, work with forces and torques, simulate various types of motions, and visualize your results in both 2D and 3D. Whether you are a seasoned researcher or a newcomer to the field, Multibody.jl will empower you to bring your ideas to life and unlock new possibilities in the fascinating world of multibody dynamics.




## Notable difference from Modelica

- The torque variable in Multibody.jl is typically called `tau` rather than `t` to not conflict with the often used independent variable `t` used to denote time.
- Multibody.jl occasionally requires the user to specify which component should act as the root of the kinematic tree. This only occurs when bodies are connected directly to force components without a joint parallel to the force component.
- In Multibody.jl, the orientation object of a [`Frame`](@ref) is accessed using he function [`ori`](@ref).



## Index
```@index
```


## Frames
```@autodocs
Modules = [Multibody]
Pages   = ["frames.jl"]
```

## Joints

A joint restricts the number of degrees of freedom (DOF) of a body. For example, a free floating body has 6 DOF, but if it is attached to a [`Revolute`](@ref) joint, the joint restricts all but one rotational degree of freedom (a revolute joint acts like a hinge). Similarily, a [`Prismatic`](@ref) joint restricts all but one translational degree of freedom (a prismatic joint acts like a slider).

A [`Spherical`](@ref) joints restricts all translational degrees of freedom, but allows all rotational degrees of freedom. It thus transmits no torque.

Some joints offer the option to add 1-dimensional components to them by providing the keyword `useAxisFlange = true`. This allows us to add, e.g., springs, dampers, sensors, and actuators to the joint.

```@autodocs
Modules = [Multibody]
Pages   = ["joints.jl"]
```

## Components

The perhaps most fundamental component is a [`Body`](@ref), this component has a single flange, `frame_a`, which is used to connect the body to other components. This component has a mass, a vector `r_cm` from `frame_a` to the center of mass, and a moment of inertia tensor `I` in the center of mass. The body can be thought of as a point mass with a moment of inertia tensor.

A mass with a shape can be modeled using a [`BodyShape`](@ref). The primary difference between a [`Body`](@ref) and a [`BodyShape`](@ref) is that the latter has an additional flange, `frame_b`, which is used to connect the body to other components. The translation between `flange_a` and `flange_b` is determined by the vector `r`. The [`BodyShape`](@ref) is suitable to model, e.g., cylinders, rods, and boxes.

A rod without a mass (just a translation), is modeled using [`FixedTranslation`](@ref).




```@autodocs
Modules = [Multibody]
Pages   = ["components.jl"]
```

## Forces
```@autodocs
Modules = [Multibody]
Pages   = ["forces.jl"]
```

## Sensors
A sensor is an object that translates quantities in the mechanical domain into causal signals which can interact with causal components from [ModelingToolkitStandardLibrary.Blocks](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/API/blocks/), such as control systems etc.

```@autodocs
Modules = [Multibody]
Pages   = ["sensors.jl"]
```

## Orientation utilities
```@autodocs
Modules = [Multibody]
Pages   = ["orientation.jl"]
```

## Interfaces
```@autodocs
Modules = [Multibody]
Pages   = ["interfaces.jl"]
```