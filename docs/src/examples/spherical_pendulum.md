# Spherical pendulum
This example models a spherical pendulum. The pivot point is modeled using a [`Spherical`](@ref) joint, the pendulum rod is modeled using a [`FixedTranslation`](@ref) and the mass is modeled using a [`Body`](@ref). In this example, we choose the body to be the root.

```@example spring_mass_system
using Multibody
using ModelingToolkit
using Plots
using SymbolicIR
using OrdinaryDiffEq

t = Multibody.t
D = Differential(t)
world = Multibody.world

@named begin
    joint = Spherical()
    bar = FixedTranslation(r = [0, -1, 0])
    body = Body(; m = 1, isroot = true)
end

connections = [connect(world.frame_b, joint.frame_a)
            connect(joint.frame_b, bar.frame_a)
            connect(bar.frame_b, body.frame_a)]

@named model = ODESystem(connections, t, systems = [world, joint, bar, body])
ssys = structural_simplify(IRSystem(model), alias_eliminate = false)

prob = ODEProblem(ssys,
                [
                    collect((body.phi)) .=> [0.5, 0.5, 0.5];
                    # collect(D.(D.(body.r_0))) .=> 0;
                    collect(D.(body.phi)) .=> 0;
                    collect(D.(body.phid)) .=> 0;
                    # collect(body.frame_a.w₃) .=> 0;
                ], (0, 10))

sol = solve(prob, Rodas4())
@assert SciMLBase.successful_retcode(sol)

plot(sol, idxs = [body.r_0...])
```
