# Spherical pendulum
This example models a spherical pendulum. The pivot point is modeled using a [`Spherical`](@ref) joint, the pendulum rod is modeled using a [`FixedTranslation`](@ref) and the mass is modeled using a [`Body`](@ref). In this example, we choose the joint to be the root (joints are often better root objects than bodies).

```@example spring_mass_system
using Multibody
using ModelingToolkit
using Plots
using JuliaSimCompiler
using OrdinaryDiffEq

t = Multibody.t
D = Differential(t)
world = Multibody.world

systems = @named begin
    joint = Spherical(enforceState=true, isroot=true, phi = 1)
    bar = FixedTranslation(r = [0, -1, 0])
    body = Body(; m = 1, isroot = false)
end

connections = [connect(world.frame_b, joint.frame_a)
            connect(joint.frame_b, bar.frame_a)
            connect(bar.frame_b, body.frame_a)]

@named model = ODESystem(connections, t, systems = [world; systems])
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys, [], (0, 10))

sol = solve(prob, Rodas4())
@assert SciMLBase.successful_retcode(sol)

plot(sol, idxs = [body.r_0...])
```
