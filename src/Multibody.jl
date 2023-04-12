module Multibody

using LinearAlgebra
using ModelingToolkit
import ModelingToolkitStandardLibrary.Mechanical.Rotational


const t = let
    (@variables t)[1]
end
const D = Differential(t)

"""
    at_variables_t(args)

Emulates the `@variables` macro but does never create array variables. Also never introuces the variable name into the calling scope.
"""
function at_variables_t(args...; default=nothing)
    xs = Symbolics.variables(args...; T = Symbolics.FnType)
    xs = map(x->x(t), xs)
    if default !== nothing
        xs = Symbolics.setdefaultval.(xs, default)
    end
    xs
end

export Orientation, RotationMatrix
include("orientation.jl")

export Frame
include("frames.jl")

include("interfaces.jl")

export World, world, FixedTranslation, Revolute, Body
include("components.jl")

export Spring
include("forces.jl")

end
