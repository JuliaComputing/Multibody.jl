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

export Orientation, RotationMatrix, ori
include("orientation.jl")

export Frame
include("frames.jl")

export PartialTwoFrames
include("interfaces.jl")

export World, world, Mounting1D, Fixed, FixedTranslation, Body, BodyShape
include("components.jl")

export Revolute, Prismatic
include("joints.jl")

export Spring, Damper, Torque, Force
include("forces.jl")

export PartialCutForceBaseSensor, BasicCutForce, BasicCutTorque, CutTorque, CutForce
include("sensors.jl")

end
