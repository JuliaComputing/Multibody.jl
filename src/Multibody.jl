module Multibody

using LinearAlgebra
using ModelingToolkit
const t = let
    (@variables t)[1]
end
const D = Differential(t)

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
