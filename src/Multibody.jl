module Multibody

using LinearAlgebra
using ModelingToolkit
using JuliaSimCompiler
import ModelingToolkitStandardLibrary.Mechanical.Rotational
import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica as Translational

export Rotational, Translational

export render, render!

"""
    scene, time = render(model, sol, t::Real; framerate = 30)
    path        = render(model, sol, timevec = range(sol.t[1], sol.t[end], step = 1 / framerate); framerate = 30)

Create a 3D animation of a multibody system

# Arguments:
- `model`: The unsimplified multibody model
- `sol`: The `ODESolution` produced by simulating the system using `solve`
- `t`: If a single number `t` is provided, the mechanism at this time is rendered and a scene is returned together with the time as an `Observable`. Modify `time[] = new_time` to change the rendering.
- `timevec`: If a vector of times is provided, an animation is created and the path to the file on disk is returned.
- `framerate`: Number of frames per second.
"""
function render end

"""
    did_render::Bool = render!(scene, ::typeof(ComponentConstructor), sys, sol, t)

Each component that can be rendered must have a `render!` method. This method is called by `render` for each component in the system.

This method is responsible for drawing the component onto the scene the way it's supposed to look at time `t` in the solution `sol`.
`t` is an Observable. It's recommended to follow the pattern
```julia
thing = @lift begin
    acces relevant coordinates from sol at time t
    create a geometric object that can be rendered
end
mesh!(scene, thing; style...)
```

# Returns
A boolean indicating whether or not the component performed any rendering. Typically, all custom methods of this function should return `true`, while the default fallback method is the only one returning false.
"""
function render! end

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
    collect([D(x) for x in collect(x)])
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

# using ModelingToolkit.SciMLBase
# import SymbolicIR: InitialType
# """
# This method exist only temporarily so that we can set default values for dummy derivatives
# """
# function SymbolicIR.SciMLBase.ODEProblem{true}(ss::SymbolicIR.ScheduledSystem, u0, tspan,
#                                                p = SciMLBase.NullParameters;
#                                                kwargs...)
#     fun = ODEFunction{true}(ss)
#     defs = ss.state.sys.info.defaults
#     if u0 !== nothing && !(u0 isa InitialType)
#         error("`u0` must be an array or a dictionary")
#     end
#     if p !== SciMLBase.NullParameters && !(p isa InitialType)
#         error("`p` must be an array or a dictionary")
#     end
#     if u0 isa InitialType && eltype(u0) <: Pair
#         u0 = Dict{SymbolicIR.IRElement, Any}(SymbolicIR.IRElement(SymbolicIR.SymbolicsConversion.extract_ir(k)) => v
#                                              for (k, v) in u0)
#         defs = merge(defs, u0)
#     end
#     if p isa InitialType && eltype(p) <: Pair
#         p = Dict{SymbolicIR.IRElement, Any}(SymbolicIR.IRElement(SymbolicIR.SymbolicsConversion.SymbolicIR.extract_ir(k)) => v
#                                             for (k, v) in p)
#         defs = merge(defs, p)
#     end
#     u0 = Float64[get(defs, SymbolicIR.to_normal(v), 0.0)
#                  for v in ModelingToolkit.states(ss)]
#     ps = Float64[defs[v] for v in ModelingToolkit.parameters(ss)]
#     ODEProblem{true}(fun, u0, tspan, ps; kwargs...)
# end

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

export Spring, Damper, SpringDamperParallel, Torque, Force
include("forces.jl")

export PartialCutForceBaseSensor, BasicCutForce, BasicCutTorque, CutTorque, CutForce
include("sensors.jl")

include("robot/path_planning.jl")
include("robot/robot_components.jl")
include("robot/FullRobot.jl")



end
