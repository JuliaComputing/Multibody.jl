module Multibody

using ModelingToolkit
const t = let
    (@variables t)[1]
end
const D = Differential(t)

include("frames.jl")

end
