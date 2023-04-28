module Multibody

using LinearAlgebra
using ModelingToolkit
using SymbolicIR
import ModelingToolkitStandardLibrary.Mechanical.Rotational
import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica as Translational

export Rotational, Translational

const t = let
    (@variables t)[1]
end
const D = Differential(t)

# Please please please Symbolic arrays, go away
function Base.broadcasted(::typeof(~), lhs::Symbolics.Arr{Num, 1}, rhs)
    collect(lhs) .~ collect(rhs)
end
function Base.broadcasted(::typeof(~), lhs, rhs::Symbolics.Arr{Num, 1})
    collect(lhs) .~ collect(rhs)
end

function Base.broadcasted(::typeof(~), lhs::Symbolics.Arr{Num, 1},
                          rhs::Symbolics.Arr{Num, 1})
    collect(lhs) .~ collect(rhs)
end

function Base.broadcasted(D::Differential, x::Symbolics.Arr{Num, 1})
    collect([D(x) for x in x])
end

"""
    at_variables_t(args)

Emulates the `@variables` macro but does never creates array variables. Also never introuces the variable names into the calling scope.
"""
function at_variables_t(args...; default = nothing)
    xs = Symbolics.variables(args...; T = Symbolics.FnType)
    xs = map(x -> x(t), xs)
    if default !== nothing
        xs = Symbolics.setdefaultval.(xs, default)
    end
    xs
end

using ModelingToolkit.SciMLBase
import SymbolicIR: InitialType

# """
# This method exist only temporarily so that we can set default values for dummy derivatives
# """
function SymbolicIR.SciMLBase.ODEProblem{true}(ss::SymbolicIR.ScheduledSystem, u0, tspan,
                                               p = SciMLBase.NullParameters;
                                               kwargs...)
    fun = ODEFunction{true}(ss)
    defs = ss.state.sys.info.defaults
    if u0 !== nothing && !(u0 isa InitialType)
        error("`u0` must be an array or a dictionary")
    end
    if p !== SciMLBase.NullParameters && !(p isa InitialType)
        error("`p` must be an array or a dictionary")
    end
    if u0 isa InitialType && eltype(u0) <: Pair
        u0 = Dict{SymbolicIR.IRElement, Any}(SymbolicIR.IRElement(SymbolicIR.SymbolicsConversion.extract_ir(k)) => v
                                             for (k, v) in u0)
        defs = merge(defs, u0)
    end
    if p isa InitialType && eltype(p) <: Pair
        p = Dict{SymbolicIR.IRElement, Any}(SymbolicIR.IRElement(SymbolicIR.SymbolicsConversion.SymbolicIR.extract_ir(k)) => v
                                            for (k, v) in p)
        defs = merge(defs, p)
    end
    u0 = Float64[get(defs, SymbolicIR.to_normal(v), 0.0)
                 for v in ModelingToolkit.states(ss)]
    ps = Float64[defs[v] for v in ModelingToolkit.parameters(ss)]
    ODEProblem{true}(fun, u0, tspan, ps; kwargs...)
end

export Orientation, RotationMatrix, ori
include("orientation.jl")

export Frame
include("frames.jl")

export PartialTwoFrames
include("interfaces.jl")

export World, world, Mounting1D, Fixed, FixedTranslation, FixedRotation, Body, BodyShape
include("components.jl")

export Revolute, Prismatic, Spherical, Universal, GearConstraint, RollingWheelJoint,
       RollingWheel, FreeMotion
include("joints.jl")

export Spring, Damper, Torque, Force
include("forces.jl")

export PartialCutForceBaseSensor, BasicCutForce, BasicCutTorque, CutTorque, CutForce
include("sensors.jl")

end
