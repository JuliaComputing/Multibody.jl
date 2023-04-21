module Multibody

using LinearAlgebra
using ModelingToolkit
import ModelingToolkitStandardLibrary.Mechanical.Rotational
import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica as Translational

export Rotational, Translational

const t = let
    (@variables t)[1]
end
const D = Differential(t)

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

export Orientation, RotationMatrix
include("orientation.jl")

export Frame
include("frames.jl")

include("interfaces.jl")

export World, world, FixedTranslation, Body
include("components.jl")

export Revolute, Prismatic, Spherical
include("joints.jl")

export Spring, Damper
include("forces.jl")

end
