# Free motions
This example demonstrates how a free-floating [`Body`](@ref) can be simulated. The body is attached to the world through a [`FreeMotion`](@ref) joint, i.e., a joint that imposes no constraints. The joint is required to add the appropriate relative states between the world and the body. We choose `enforceStates = true` and `isroot = true` in the [`FreeMotion`](@ref) constructor.

```@example FREE_MOTION
using Multibody
using ModelingToolkit
using Plots
using SymbolicIR
using OrdinaryDiffEq

t = Multibody.t
D = Differential(t)
world = Multibody.world

@named freeMotion = FreeMotion(enforceStates = true, isroot = true)
@named body = Body(m = 1)

eqs = [connect(world.frame_b, freeMotion.frame_a)
       connect(freeMotion.frame_b, body.frame_a)]

@named model = ODESystem(eqs, t,
                         systems = [world;
                                    freeMotion;
                                    body])
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys, [
    D.(freeMotion.r_rel_a) .=> randn();
    D.(D.(freeMotion.r_rel_a)) .=> randn();
    D.(freeMotion.phi) .=> randn();
    D.(D.(freeMotion.phi)) .=> randn();
], (0, 10))

sol = solve(prob, Rodas4())
plot(sol, idxs = body.r_0[2], title="Free falling body")

# Plot analytical solution
tvec = sol.t
plot(-9.81/2 .* tvec .^ 2, lab="Analytical solution")
```

The figure indicates that the body is falling freely, experiencing a constant acceleration of -9.81 m/s² in the ``y`` direction, corresponding to the gravity parameters of the `world`:
```@example FREE_MOTION
show(stdout, MIME"text/plain"(), world)
nothing # hide
```