module Multibody

using ModelingToolkit
const t = let
    (@variables t)[1]
end
const D = Differential(t)

export Frame, Orientation
include("frames.jl")

include("interfaces.jl")

export World, world, FixedTranslation, Revolute, Body
include("components.jl")

export Spring
include("forces.jl")

end
