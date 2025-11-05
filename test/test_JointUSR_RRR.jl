using Test
using Multibody
using ModelingToolkit
import ModelingToolkitStandardLibrary.Mechanical.Rotational
# using Plots
using OrdinaryDiffEq
using LinearAlgebra
# using JuliaSimCompiler

world = Multibody.world

t = Multibody.t
D = Differential(t)

# ==============================================================================
## A JointUSR is connected to a prismatic joint, with a Body at their common tip
# ==============================================================================
@component function TestUSR(; name)
    systems = @named begin
        world = World()
        j1 = JointUSR(positive_branch=true, use_arrays=false)
        fixed = FixedTranslation(r=[1,0,0])
        b1 = Body(m=1)
        p1 = Prismatic(state_priority=100, n = [1, 0, 0])
    end

    pars = @parameters begin
    end

    vars = @variables begin
    end

    equations = [
        connect(world.frame_b, j1.frame_a, fixed.frame_a)
        connect(fixed.frame_b, p1.frame_a)
        connect(p1.frame_b, j1.frame_b)
        connect(j1.frame_im, b1.frame_a)
    ]

    return System(equations, t; name, systems)
end

@named model = TestUSR()
model = complete(model)
ssys = structural_simplify(multibody(model))
@test length(unknowns(ssys)) == 2
##

prob = ODEProblem(ssys, [model.p1.v => 0.0], (0.0, 1.4))
sol = solve(prob, FBDF(autodiff=true))
@test sol(1.0, idxs=model.p1.s) ≈ -2.8 rtol=0.01 # test vs. OpenModelica



# ==============================================================================
## A JointRRR is connected to a prismatic joint, with a Body at their common tip
# ==============================================================================
@component function TestRRR(; name)
    systems = @named begin
        world = World()
        j1 = JointRRR(positive_branch=true)
        fixed = FixedTranslation(r=[1,0,0])
        b1 = Body(m=1)
        p1 = Prismatic(state_priority=100, n = [1, 0, 0])
    end

    pars = @parameters begin
    end

    vars = @variables begin
    end

    equations = [
        connect(world.frame_b, j1.frame_a, fixed.frame_a)
        connect(fixed.frame_b, p1.frame_a)
        connect(p1.frame_b, j1.frame_b)
        connect(j1.frame_im, b1.frame_a)
    ]

    return System(equations, t; name, systems)
end

@named model = TestRRR()
model = complete(model)
ssys = structural_simplify(multibody(model))
@test length(unknowns(ssys)) == 2
##

prob = ODEProblem(ssys, [model.p1.v => 0.0], (0.0, 1.4))
sol = solve(prob, FBDF(autodiff=true))
@test sol(1.0, idxs=model.p1.s) ≈ -2.8 rtol=0.01 # test vs. OpenModelica

