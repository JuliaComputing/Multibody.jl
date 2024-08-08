using Test
using Multibody
using ModelingToolkit
import ModelingToolkitStandardLibrary.Mechanical.Rotational
# using Plots
using OrdinaryDiffEq
using LinearAlgebra
using JuliaSimCompiler

t = Multibody.t
D = Differential(t)
world = Multibody.world
W(args...; kwargs...) = Multibody.world

# A JointUSR is connected to a prismatic joint, with a Body at their common tip
@mtkmodel TestUSR begin
    @components begin
        world = W()
        j1 = JointUSR(positive_branch=true, use_arrays=false)
        fixed = FixedTranslation(r=[1,0,0])
        b1 = Body(isroot=false, neg_w=true)
        p1 = Prismatic(state_priority=100)
    end
    @equations begin
        connect(world.frame_b, j1.frame_a, fixed.frame_a)
        connect(fixed.frame_b, p1.frame_a)
        connect(p1.frame_b, j1.frame_b)
        connect(j1.frame_im, b1.frame_a)
    end
end

@named model = TestUSR()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, [], (0.0, 1.0))
sol = solve(prob, FBDF(autodiff=true))


# NOTE: I was working on trying to register the compute_angle function so that there are no symbolic arguments left in the generated code

# foo(x::AbstractArray, p::AbstractArray) = sum(-p .* x .^ 2)
# @register_symbolic foo(x::AbstractArray, p::AbstractArray)
# onetwo = [1, 2]

# @mtkmodel ArrayX begin
#     @variables begin
#         x(t)[1:2] = ones(2)
#     end
#     @parameters begin
#         p[1:2] = onetwo
#     end
#     begin
#         x = collect(x)
#     end
#     @equations begin
#         D(x[1]) ~ foo(x, p)
#         D(x[2]) ~ foo(x, p)
#     end
# end

# @named model = ArrayX()
# model = complete(model)
# ssys = structural_simplify(IRSystem(model))
# prob = ODEProblem(ssys, [], (0.0, 1.0))
# sol = solve(prob, FBDF())


# @variables begin
#     x(t)[1:2] = 1
# end
# @parameters begin
#     p[1:2] = onetwo
# end
# begin
#     x = collect(x)
# end

# foo(x, p)
# # @equations begin
# #     D(x[1]) ~ foo(x, p)[1]
# #     D(x[2]) ~ foo(x, p)[2]
# # end

