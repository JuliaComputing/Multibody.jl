# Modeling kinematic loops

Kinematic loops can be difficult to simulate since they introduce an over-constrained system. This tutorial demonstrates how to handle a few common cases of kinematic loops.

## A planar kinematic loop

A planar loop is one where the loop is confined to a plane, i.e., all joints in the loop have parallel rotation axes. To simulate a mechanism with such a loop, we break the kinematic loop by replacing one of the revolute joints with a [`RevolutePlanarLoopConstraint`](@ref). The model below contains four bars connected by revolute joints, forming a planar loop. In order to make the animation interesting, we attach dampers to two of the joints such that the mechanism will oscillate for a while before coming to rest.

Perhaps surprisingly, we use 5 joints in total for this mechanism. If we had used four joints only and connected the first frame of the first joint to the world, the mechanism would not be free to rotate around the world frame. We thus have two joints connected to the world frame below.

```@example kinloop
using Multibody
using ModelingToolkit
import ModelingToolkitStandardLibrary.Mechanical.Rotational
using Plots
using OrdinaryDiffEq
using LinearAlgebra
using JuliaSimCompiler

t = Multibody.t
D = Differential(t)
world = Multibody.world

l = 1.5
systems = @named begin
    j1 = Revolute(axisflange=true) # We use an axis flange to attach a damper
    j2 = Revolute(axisflange=true)
    j3 = Revolute()
    j4 = RevolutePlanarLoopConstraint()
    j5 = Revolute()
    b1 = BodyShape(m=1, r = [l, 0, 0], radius=0.03)
    b2 = BodyShape(m=1, r = [0.0, l, 0], radius=0.03)
    b3 = BodyShape(m=1, r = [-l, 0, 0], radius=0.03)
    b4 = BodyShape(m=1, r = [0.0, -l, 0], radius=0.03)
    damper1 = Rotational.Damper(d=0.1)
    damper2 = Rotational.Damper(d=0.1)
end

connections = [
    connect(world.frame_b, j1.frame_a)
    
    connect(j1.frame_b, b1.frame_a)
    connect(b1.frame_b, j2.frame_a)
    connect(j2.frame_b, b2.frame_a)
    connect(b2.frame_b, j3.frame_a)
    connect(j3.frame_b, b3.frame_a)
    connect(b3.frame_b, j4.frame_a)
    
    connect(j4.frame_b, b4.frame_a)
    
    connect(b4.frame_b, j5.frame_a)
    connect(j5.frame_b, world.frame_b) 
    # We need 5 joints since j1.frame_a is rigidly attached to the world, and b4 closing the loop would thus not be able to rotate around j1.

    connect(j1.axis, damper1.flange_a)
    connect(j1.support, damper1.flange_b)

    connect(j2.axis, damper2.flange_a)
    connect(j2.support, damper2.flange_b)
    
]
@named fourbar = ODESystem(connections, t, systems = [world; systems])
fourbar = complete(fourbar)
ssys = structural_simplify(IRSystem(fourbar))
prob = ODEProblem(ssys, [fourbar.j1.phi => 0.1], (0.0, 10.0))
sol = solve(prob, FBDF(autodiff=true))

plot(
    plot(sol, idxs = [j1.phi, j2.phi, j3.phi]),
    plot(sol, idxs = [j1.w, j2.w, j3.w]),
)
```

```@example kinloop
using Test
@test SciMLBase.successful_retcode(sol)
```


### 3D animation
Multibody.jl supports automatic 3D rendering of mechanisms, we use this feature to illustrate the result of the simulation below:

```@example kinloop
import GLMakie
Multibody.render(fourbar, sol, 0:0.05:10; x=4, y=-1, z=4, lookat=[0, -1, 0], filename = "fourbar.gif") # Use "fourbar.mp4" for a video file
nothing # hide
```

![animation](fourbar.gif)


## Using cut joints

The mechanism below is another instance of a 4-bar linkage, this time with 6 revolute joints, 1 prismatic joint and 4 bodies. In order to simulate this mechanism, the user must
1. Use the `iscut=true` keyword argument to one of the `Revolute` joints to indicate that the joint is a cut joint. A cut joint behaves similarly to a regular joint, but it introduces fewer constraints in order to avoid the otherwise over-constrained system resulting from closing the kinematic loop. While almost any joint can be chosen as the cut joint, it might be worthwhile experimenting with this choice in order to get an efficient representation. In this example, cutting `j5` produces an 8-dimensional state realization, while all other joints result in a 17-dimensional state.
2. Increase the `state_priority` of the joint `j1` above the default joint priority 3. This encourages the model compiler to choose the joint coordinate of `j1` as state variable. The joint coordinate of `j1` is the only coordinate that uniquely determines the configuration of the mechanism. The choice of any other joint coordinate would lead to a singular representation in at least one configuration of the mechanism. The joint `j1` is the revolute joint located in the origin, see the animation below where this joint is made larger than the others.


To drive the mechanism, we set the initial velocity of the joint j1 to some non-zero value.

```@example kinloop
systems = @named begin
    j1 = Revolute(n = [1, 0, 0], w0 = 5.235987755982989, state_priority=10.0, radius=0.1f0) # Increase state priority to ensure that this joint coordinate is chosen as state variable
    j2 = Prismatic(n = [1, 0, 0], s0 = -0.2)
    b1 = BodyShape(r = [0, 0.5, 0.1], radius=0.03)
    b2 = BodyShape(r = [0, 0.2, 0], radius=0.03)
    b3 = BodyShape(r = [-1, 0.3, 0.1], radius=0.03)
    rev = Revolute(n = [0, 1, 0])
    rev1 = Revolute()
    j3 = Revolute(n = [1, 0, 0])
    j4 = Revolute(n = [0, 1, 0])
    j5 = Revolute(n = [0, 0, 1], iscut=true)
    b0 = FixedTranslation(r = [1.2, 0, 0], radius=0)
end

connections = [connect(j2.frame_b, b2.frame_a)
               connect(j1.frame_b, b1.frame_a)
               connect(rev.frame_a, b2.frame_b)
               connect(rev.frame_b, rev1.frame_a)
               connect(rev1.frame_b, b3.frame_a)
               connect(world.frame_b, j1.frame_a)
               connect(b1.frame_b, j3.frame_a)
               connect(j3.frame_b, j4.frame_a)
               connect(j4.frame_b, j5.frame_a)
               connect(j5.frame_b, b3.frame_b)
               connect(b0.frame_a, world.frame_b)
               connect(b0.frame_b, j2.frame_a)
               ]
@named fourbar2 = ODESystem(connections, t, systems = [world; systems])
fourbar2 = complete(fourbar2)
ssys = structural_simplify(IRSystem(fourbar2))

prob = ODEProblem(ssys, [], (0.0, 1.4399)) # The end time is chosen to make the animation below appear to loop forever

sol = solve(prob, FBDF(autodiff=true));
@test SciMLBase.successful_retcode(sol)
plot(sol, idxs=[j2.s]) # Plot the joint coordinate of the prismatic joint (green in the animation below)
```

```@example kinloop
import GLMakie
Multibody.render(fourbar2, sol; x=-2, y = 2, z = 3, filename = "fourbar2.gif") # Use "fourbar2.mp4" for a video file
nothing # hide
```

![animation](fourbar2.gif)