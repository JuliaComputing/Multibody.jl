# Gearbox

This example models a gearbox in two different ways:
1. Using the 3D [`GearConstraint`](@ref) component from the multibody library.
2. Using the 1D [`IdealGear`](@ref) component from the `Rotational` submodule, together with a [`Mounting1D`](@ref) component.

The [`GearConstraint`](@ref) has two rotational axes which do not have to be parallel. If wou want to select rotational axes, use the keyword arguments `n_a` and `n_b` to [`GearConstraint`](@ref).

```@example gearbox
using Multibody
using ModelingToolkit
using Plots
using JuliaSimCompiler
using OrdinaryDiffEq

t = Multibody.t
D = Differential(t)
world = Multibody.world

systems = @named begin
    gearConstraint = GearConstraint(; ratio = 10)
    cyl1 = Body(; m = 1, r_cm = [0.4, 0, 0])
    cyl2 = Body(; m = 1, r_cm = [0.4, 0, 0])
    torque1 = Torque(resolveInFrame = :frame_b)
    # sine[1:3] = Blocks.Sine(frequency = 1)
    fixed = Fixed() # TODO: implement
    inertia1 = Rotational.Inertia(J = cyl1.I_11)
    idealGear = Rotational.IdealGear(ratio = 10, use_support = true)
    inertia2 = Rotational.Inertia(J = cyl2.I_11)
    torque2 = Rotational.Torque(use_support = true)
    mounting1D = Mounting1D()
end

eqs = [connect(world.frame_b, gearConstraint.bearing)
       connect(cyl1.frame_a, gearConstraint.frame_a)
       connect(gearConstraint.frame_b, cyl2.frame_a)
       connect(torque1.frame_b, cyl1.frame_a)
       connect(torque1.frame_a, world.frame_b)
       # connect(sine.output, torque1.torque)
       torque1.torque.u .~ [2sin(t), 0, 0]
       connect(inertia1.flange_b, idealGear.flange_a)
       connect(idealGear.flange_b, inertia2.flange_a)
       connect(torque2.flange, inertia1.flange_a)
       # connect(sine.output, torque2.tau)
       torque2.tau.u ~ 2sin(t)
       connect(mounting1D.flange_b, idealGear.support)
       connect(mounting1D.flange_b, torque2.support)
       connect(fixed.frame_b, mounting1D.frame_a)]

@named model = ODESystem(eqs, t, systems = [world; systems])
cm = complete(model)
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, [
    D(cm.idealGear.phi_b) => 0
], (0, 10))
sol = solve(prob, Rodas4())
plot(sol)
```