using ModelingToolkit
using Multibody
using Test
using LinearAlgebra
using JuliaSimCompiler
using OrdinaryDiffEq
using Random

using Multibody: Rt2T, T2t
t = Multibody.t
D = Differential(t)
W(args...; kwargs...) = Multibody.world

Random.seed!(0)
n = 1
M = 1.0I(6*n)
K = rand(6n, 6n)
K = K'K

D = 0.5M + 0.5K


if false
    n = 2
    M = cat(M,M, dims=(1,2))
    K = cat(K,K, dims=(1,2))
    D = cat(D,D, dims=(1,2))
end


T = [Rt2T(I(3), rand(3)) for i = 1:n]
@named cb = CraigBampton(; M, K, T)


@mtkmodel CBTest begin
    @components begin
        world = W()
        cb = CraigBampton(; M, D, K, T)
        # rev = Spherical()
        # body = Body(m=1)
    end
    @equations begin
        connect(world.frame_b, cb.common_reference_frame)
        # connect(world.frame_b, cb.interface1.frame_b)
        # connect(cb.interface1.frame_b, body.frame_a)

        # connect(world.frame_b, rev.frame_a)
        # connect(rev.frame_b, cb.common_reference_frame)
    end
end

@named model = CBTest()
model = complete(model)
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys, [
    collect(model.cb.x) .=> randn(6n);
    # model.world.g => 100;
], (0.0, 2000.0))
sol = solve(prob, Rodas4())

using Plots
plot(sol)

##

plot(sol, idxs=collect(model.cb.f))


@test sol(prob.tspan[2], idxs=collect(model.cb.x)) ≈ zeros(6n) atol = 1e-6
@test sol(prob.tspan[2], idxs=collect(model.cb.interface1.frame_b.r_0)) ≈ T2t(T[1])