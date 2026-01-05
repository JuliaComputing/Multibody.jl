using Test
using ModelingToolkit, Multibody, LightXML, Graphs, MetaGraphsNext, JuliaFormatter
cd(@__DIR__)
# URDF = Base.get_extension(Multibody, :URDF)

urdf2multibody(joinpath(dirname(pathof(Multibody)), "..", "test", "urdf", "doublependulum.urdf"), extras=false, out="doublependulum.jl", modelname="DoublePendulum")

include("doublependulum.jl")



using ModelingToolkit, Multibody, OrdinaryDiffEq#, Plots
import ModelingToolkit: t_nounits as t, D_nounits as D

@named model = DoublePendulum()
model = complete(model)
ssys = multibody(model)
prob = ODEProblem(ssys, [], (0.0, 10.0))
sol = solve(prob, Tsit5())
@test SciMLBase.successful_retcode(sol)
@test !iszero(sol.u[end])
# plot(sol) |> display