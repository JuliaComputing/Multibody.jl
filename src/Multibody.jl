# Find variables that are both array form and scalarized / collected
# foreach(println, sort(unknowns(IRSystem(model)), by=string))
module Multibody

using LinearAlgebra
using ModelingToolkit
using JuliaSimCompiler
import ModelingToolkitStandardLibrary.Mechanical.Rotational
import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica as Translational
using StaticArrays
export Rotational, Translational

export render, render!

"""
Find parameters that occur both scalarized and not scalarized
"""
function find_arry_problems(model)
    # foreach(println, sort(unknowns(IRSystem(model)), by=string))
    xs = string.(unknowns(IRSystem(model)))
    for x in xs
        endswith(x, ']') && continue # Only look at non-array vars
        l = ncodeunits(x)
        inds = findall(y->startswith(y, x*"["), xs)
        isempty(inds) && continue
        println(x)
        println(xs[inds])
    end
end

"""
    benchmark_f(prob)
    
Benchmark the invocation of the function `prob.f` on the initial condition `prob.u0`
"""
function benchmark_f(prob)
    x = prob.u0
    dx = similar(x)
    p = prob.p
    if !isdefined(Main, :btime)
        @eval Main using BenchmarkTools
    end
    Main.@btime $(prob.f)($dx, $x, $p, 0.0)
end

"""
    scene, time = render(model, sol, t::Real; framerate = 30, traces = [])
    path        = render(model, sol, timevec = range(sol.t[1], sol.t[end], step = 1 / framerate); framerate = 30, timescale=1, display=false, loop=1)

Create a 3D animation of a multibody system

# Arguments:
- `model`: The _unsimplified_ multibody model, i.e., this is the model _before_ any call to `structural_simplify`.
- `sol`: The `ODESolution` produced by simulating the system using `solve`
- `t`: If a single number `t` is provided, the mechanism at this time is rendered and a scene is returned together with the time as an `Observable`. Modify `time[] = new_time` to change the rendering.
- `timevec`: If a vector of times is provided, an animation is created and the path to the file on disk is returned.
- `framerate`: Number of frames per second.
- `timescale`: Scaling of the time vector. This argument can be made to speed up animations (`timescale < 1`) or slow them down (`timescale > 1`). A value of `timescale = 2` will be 2x slower than real time.
- `loop`: The animation will be looped this many times. Please note: looping the animation using this argument is only recommended when `display = true` for camera manipulation purposes. When the camera is not manipulated, looping the animation by other means is recommended to avoid an increase in the file size.
- `filename` controls the name and the file type of the resulting animation
- `traces`: An optional array of frames to show the trace of.
- `show_axis = false`: Whether or not to show the plot axes, including background grid.

# Camera control
The following keyword arguments are available to control the camera pose:
- `x = 2`
- `y = 0.5`
- `z = 2`
- `lookat = [0,0,0]`: a three-vector of coordinates indicating the point at which the camera looks.
- `up = [0,1,0]`: A vector indicating the direction that is up.
- `display`: if `true`, the figure will be displayed during the recording process and time will advance in real-time. This allows the user to manipulate the camera options using the mouse during the recording.

See also [`loop_render`](@ref)
"""
function render end

"""
    loop_render(model, sol; framerate = 30, timescale = 1, max_loop = 5, kwargs...)

Similar to the method of [`render`](@ref) that produces an animation, but instead opens an interactive window where the time is automatically advanced in real time. This allows the user to manually manipulate the camera using the mouse is a live animation.
"""
function loop_render end

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
    (@independent_variables t)[1]
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
function at_variables_t(args...; default = nothing, state_priority = nothing)
    xs = Symbolics.variables(args...; T = Symbolics.FnType)
    xs = map(x -> x(t), xs)
    if default !== nothing
        xs = Symbolics.setdefaultval.(xs, default)
    end
    if state_priority !== nothing
        xs = Symbolics.setmetadata.(xs, ModelingToolkit.VariableStatePriority, state_priority)
    end
    xs
end

encode(s) = Float64.(codeunits(s)) # Used to store strings as vectors of floats in parameters. useful for providing paths to shapefiles for 3D rendering
decode(s) = String(UInt8.(s))

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

export Orientation, RotationMatrix, ori, get_rot, get_trans, get_frame
include("orientation.jl")

export Frame
include("frames.jl")

export PartialTwoFrames
include("interfaces.jl")

export World, world, Mounting1D, Fixed, FixedTranslation, FixedRotation, Body, BodyShape, BodyCylinder, BodyBox, Rope
include("components.jl")

export Revolute, Prismatic, Planar, Spherical, Universal,
GearConstraint, FreeMotion, RevolutePlanarLoopConstraint, Cylindrical
include("joints.jl")

export SphericalSpherical, UniversalSpherical, JointUSR, JointRRR
include("fancy_joints.jl")

export RollingWheelJoint, RollingWheel, RollingWheelSet, RollingWheelSetJoint, RollingConstraintVerticalWheel
include("wheels.jl")

export Spring, Damper, SpringDamperParallel, Torque, Force, WorldForce, WorldTorque
include("forces.jl")

export PartialCutForceBaseSensor, BasicCutForce, BasicCutTorque, CutTorque, CutForce, Power
include("sensors.jl")

export point_to_point, traj5, KinematicPTP, Kinematic5
include("robot/path_planning.jl")
include("robot/robot_components.jl")
include("robot/FullRobot.jl")



end
