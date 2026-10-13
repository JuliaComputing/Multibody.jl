using ModelingToolkit
using Multibody
using Test
# using JuliaSimCompiler
using OrdinaryDiffEq
using LinearAlgebra
isdefined(Main, :t) || (t = Multibody.t)
isdefined(Main, :D) || (D = Differential(t))
doplot() = false
world = Multibody.world


# linsys = (; allow_symbolic = true, inline_linear_sccs = true, analytical_linear_scc_limit = 10, reassemble_alg = StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = true, analytical_linear_scc_limit = 10))

linsys = (; reassemble_alg = StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = true, analytical_linear_scc_limit = 1))

@testset "Multibody tests" begin
@testset "world" begin
    @info "Testing world"
    include("test_world.jl")
end

@testset "urdf" begin
    @info "Testing urdf"
    include("test_urdf.jl")
end

@testset "traj" begin
    @info "Testing traj"
    include("test_traj.jl")
end

@testset "robot" begin
    @info "Testing robot"
    include("test_robot.jl")
end

@testset "orientation_getters" begin
    @info "Testing orientation_getters"
    include("test_orientation_getters.jl")
end

@testset "quaternions" begin
    @info "Testing quaternions"
    include("test_quaternions.jl")
end

@testset "worldforces" begin
    @info "Testing worldforces"
    include("test_worldforces.jl")
end

@testset "PlanarMechanics" begin
    @info "Testing PlanarMechanics"
    include("test_PlanarMechanics.jl")
end

@testset "basic_tests" begin
    @info "Testing basic_tests"
    include("basic_tests.jl")
end
end
