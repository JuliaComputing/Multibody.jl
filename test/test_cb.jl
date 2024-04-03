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

## Test with 1 DOF


# M is generally Diagonal(mI, J)

M = Diagonal([1, 1, 1, 1e-6, 1e-6, 1e-6])
K = Diagonal([1, 1000, 1000, 1000, 1000, 1000])
D = 0.0(M + K)


T = [Rt2T(I(3), zeros(3)) for i = 1:n]
@named cb = CraigBampton(; M, K, T)


@mtkmodel CBTest begin
    @components begin
        world = W()
        cb = CraigBampton(; M, K, T)
        # body = Body(m=1)
    end
    @equations begin
        connect(world.frame_b, cb.common_reference_frame)
        # connect(cb.common_reference_frame, body.frame_a)
    end
end

@named model = CBTest()
model = complete(model)
ssys = structural_simplify(IRSystem(model))



x0 = [1; zeros(5)]
prob = ODEProblem(ssys, [
    collect(model.cb.x) .=> x0;
    # model.world.g => 100;
], (0.0, 20.0))
sol = solve(prob, Rodas4())

using Plots
plot(sol)



#
plot(sol, idxs=collect(model.world.frame_b.f))
plot(sol, idxs=collect(model.cb.f))
plot(sol, idxs=collect(model.cb.interface1.frame_a.f))
plot(sol, idxs=collect(model.cb.common_reference_frame.f))


# After one period (2pi) we are back at the start
@test sol(2pi, idxs=collect(model.cb.x)) ≈ x0 atol = 1e-3
@test sol(2pi, idxs=collect(model.cb.v)) ≈ zeros(6) atol = 1e-3
@test sol(0, idxs=collect(model.cb.interface1.frame_b.r_0)) ≈ T2t(T[1]) + x0[1:3] atol = 1e-3


@test sol(pi/2, idxs=collect(model.cb.x[1:3])) ≈ sol(pi/2, idxs=collect(model.cb.interface1.frame_b.r_0)) atol = 1e-3 # The interface node coords are equal to x in this case when the node transformation is 0










##

Random.seed!(0)
n = 1
M = 1.0I(6*n)
K = rand(6n, 6n)
K = K'K

D = 0.05(M + K)


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
    end
    @equations begin
        connect(world.frame_b, cb.common_reference_frame)
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

@test sol(1, idxs=collect(model.cb.f)) != zeros(6n) 
@test sol(prob.tspan[2], idxs=collect(model.cb.x)) ≈ zeros(6n) atol = 1e-6
@test sol(prob.tspan[2], idxs=collect(model.cb.interface1.frame_b.r_0)) ≈ T2t(T[1])



## test FlexibleCoupling

M = Diagonal([1, 1, 1, 1e-6, 1e-6, 1e-6])
K = Diagonal([1, 1, 1, 1, 1, 1])
# D = 0.0(M + K)

T = Rt2T(I(3), zeros(3))

@mtkmodel TestFlex begin
    @components begin
        world = W()
        fc = Multibody.FlexibleCoupling(; K, T)
        body = Body(m=1)
    end
    @equations begin
        connect(world.frame_b, fc.frame_a)
        connect(fc.frame_b, body.frame_a)
    end
end

@named model = TestFlex()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, [
    # collect(model.fc.x) .=> randn(6);
    # model.world.g => 100;
], (0.0, 1.0))
sol = solve(prob, Rodas4())
plot(sol, idxs=[
    collect(model.fc.f[2])
    collect(model.fc.x[2])
    collect(model.body.r_0[2])
    # collect(model.fc.frame_a.r_0[2])
    collect(model.fc.frame_b.r_0[2])

], layout=1)

# plot(sol, idxs=collect(model.fc.f))

