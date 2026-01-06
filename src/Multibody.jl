# Find variables that are both array form and scalarized / collected
# foreach(println, sort(unknowns(IRSystem(model)), by=string))

# Problem: ERROR: Could not evaluate value of parameter rod₊body₊I. Missing values for variables in expression rod₊I.
# solution: explicitly pass I in pars to System constructor
module Multibody
# Find variables that are both array form and scalarized / collected
# foreach(println, sort(unknowns(IRSystem(model)), by=string))
using LinearAlgebra
using ModelingToolkit
import ModelingToolkitStandardLibrary.Mechanical.Rotational
import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica as Translational
import ModelingToolkitStandardLibrary.Blocks
using SparseArrays
using StaticArrays
export Rotational, Translational

export render, render!
export subs_constants


## JuliaSimCompiler transition helpers
export IRSystem
IRSystem(x) = x


"""
A helper function that calls `collect` on all parameters in a vector of parameters in a smart way
"""
function collect_all(pars)
    pc = map(pars) do p
        if p isa AbstractArray || !(p isa SymbolicUtils.BasicSymbolic{<:Real})
            collect(p)
        else
            p
        end
    end
    reduce(vcat, pc)
end

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
    @eval Main @btime $(prob.f)($dx, $x, $p, 0.0)
end

"get_systemtype(sys): Get the constructor of a component for dispatch purposes. This only supports components that have the `gui_metadata` property set. If no metadata is available, nothing is returned."
function get_systemtype(sys)
    meta = getfield(sys, :gui_metadata)
    meta === nothing && return nothing
    eval(meta.type)
end

"""
    subs_constants(model, c=[0, 1]; ssys = multibody(model), defs = defaults(model))

A value-dependent compile-time optimization. Replace parameters in the model that have a default value contained in `c` with their value.

Many parameters in multibody models have a sparse structure, e.g., rotation axes like `[1, 0, 0]`, or diagonal mass matrices etc. Specializing the generated code to this structure may increase performance, sometimes up to 2x. This function replaces all parameters that have a value contained in `c`, which defaults to `[0, 1]`, with their value.

This performance optimization is primarily beneficial when the runtime of the simulation exceeds the compilation time of the model. For non-multibody models, sparse parameter structure is less common and this optimization is this not likely to be beneficial.

# Drawbacks
There are two main drawbacks to performing this optimization
- Parameters that have been replaced **cannot be changed** after the optimization has been performed, without recompiling the model using `multibody`.
- The value of the repalced parameters are no longer accessible in the solution object. This typically means that **3D animations cannot be rendered** for models with replaced parameters.

# Example usage
```
@named robot = Robot6DOF()
robot = complete(robot)
ssys = multibody(robot)
ssys = Multibody.subs_constants(model, [0, 1]; ssys)
```

If this optimization is to be performed repeatedly for several simulations of the same model, the indices of the substituted parameters can be stored and reused, call the lower-level function `Multibody.find_defaults_with_val` with the same signature as this function to obtain these indices, and then call `JuliaSimCompiler.freeze_parameters(ssys, inds)` with the indices to freeze the parameters.
"""
function subs_constants(model, c=[0, 1]; ssys = multibody(model), kwargs...)
    inds = find_defaults_with_val(model, c; ssys, kwargs...)
    # ssys = JuliaSimCompiler.freeze_parameters(ssys, inds)
    @error "JuliaSimCompiler.freeze_parameters is no longer available. This optimization is currently disabled."
    return ssys
end

function find_defaults_with_val(model, c=[0, 1]; defs = defaults(model), ssys = multibody(model))
    kvpairs = map(collect(pairs(defs))) do (key, val)
        if val isa AbstractArray
            string.(collect(key)) .=> collect(val)
        else
            string(key) => val
        end
    end
    kvpairs = reduce(vcat, kvpairs)

    p = parameters(ssys)
    stringp = string.(p)

    inds = map(c) do ci
        sdefs = map(kvpairs) do (key, val)
            if (val isa Bool) && !(ci isa Bool)
                # We do not want to treat bools as an integer to prevent subsituting false when c contains 0
                return key => false
            end
            key => isequal(val, ci)
        end |> Dict
        map(eachindex(p)) do i
            get(sdefs, stringp[i], false)
        end |> findall
    end
    sort(reduce(vcat, inds))
end


"""
    guesses_for_all_parameters(ssys, guesses = Dict{Any, Any}())

Generate a dictionary of NaN guesses for all parameters that do not already have a guess specified in `guesses` or a default guess in the system `ssys`. This is useful as a debugging tool when parameter initialization or optimization is not giving the expected results.
"""
function guesses_for_all_parameters(ssys, guesses = Dict{Any, Any}())
    guesses = Dict(deepcopy(guesses))
    for p in ModelingToolkit.parameters(ssys)
        haskey(guesses, p) && continue
        ModelingToolkit.hasguess(p) && continue
        ModelingToolkit.symtype(p) <: Number || continue
        guesses[p] = NaN
    end
    return guesses
end



export multibody

"""
    multibody(model)

Perform validity checks on the model, such as the precense of exactly one world component in the top level of the model, and call `mtkcompile` with simplification options suitable for multibody systems.
"""
function multibody(model, level=0; reassemble_alg = StructuralTransformations.DefaultReassembleAlgorithm(; inline_linear_sccs = true, analytical_linear_scc_limit = 3), kwargs...)
    found_world = false
    found_planar = false
    for subsys in getfield(model, :systems)
        system_type = get_systemtype(subsys)
        subsys_ns = getproperty(model, getfield(subsys, :name))
        isworld = system_type == World
        isplanar = system_type !== nothing && parentmodule(system_type) == PlanarMechanics
        found_world = found_world || isworld
        found_planar = found_planar || isplanar
        multibody(subsys_ns, level + 1)
    end
    if level == 0 && !found_world && !found_planar
        @warn("No world found in the top level of the model, this may lead to missing equations")
    elseif level != 0 && found_world
        @warn("World found in a non-top level component ($(nameof(model))) of the model, this may lead to extra equations. Consider using the component `Fixed` instead of `World` in component models.")
    end
    if level == 0
        return mtkcompile(model; reassemble_alg, kwargs...)
    else
        return nothing
    end
end

"""
    scene       = render(model, prob)
    scene, time = render(model, sol, t::Real; framerate = 30, traces = [])
    path        = render(model, sol, timevec = range(sol.t[1], sol.t[end], step = 1 / framerate); framerate = 30, timescale=1, display=false, loop=1)

Create a 3D animation of a multibody system

# Arguments:
- `model`: The _unsimplified_ multibody model, i.e., this is the model _before_ any call to `structural_simplify`.
- `prob`: If an `ODEProblem` is passed, a static rendering of the system at the initial condition is generated.
- `sol`: If an `ODESolution` produced by simulating the system using `solve` is passed, an animation or dynamic rendering of the system is generated.
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

"""
    urdf2multibody(filename::AbstractString; extras=false, out=nothing, worldconnection = :rigid)

Translate a URDF file into a Multibody model. Only available if LightXML.jl, Graphs.jl, MetaGraphs.jl and JuliaFormatter.jl are manually installed and loaded by the user.

Example usage:
```
using Multibody, ModelingToolkit, LightXML, Graphs, MetaGraphsNext, JuliaFormatter
urdf2multibody(joinpath(dirname(pathof(Multibody)), "..", "test/doublependulum.urdf"), extras=true, out="/tmp/urdf_import.jl")
```

## Keyword arguments
- `extras=false`: If `true`, the generated code will include package imports, a simulation of the model and a rendering of the model.
- `out=nothing`: If provided, the generated code will be written to this file, otherwise the string will only be returned.
- `worldconnection=:rigid`: If `:rigid`, the world frame will be connected to the root link with a rigid connection. If a joint constructor is provided, this component will be instantiated and the root link is connected to the world through this, e.g., `worldconnection = FreeMotion`, `()->Prismatic(n=[0, 1, 0])` etc.
`render_fixed = false`: Whether or not to render "fixed" joints. These joints aren't actually joints (no degrees of freedom), they are translated to FixedTranslation or FixedRotation components.
"""
function urdf2multibody end
export urdf2multibody, URDFRevolute, URDFPrismatic, NullJoint

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
    xs = Symbolics.variables(args...; T = Symbolics.FnType{Tuple, Real, Nothing})
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

export World, world, Mounting1D, Fixed, Position, Pose, FixedTranslation, FixedRotation, Body, BodyShape, BodyCylinder, BodyBox, Rope
include("components.jl")

export Revolute, Prismatic, Planar, Spherical, Universal,
GearConstraint, FreeMotion, RevolutePlanarLoopConstraint, Cylindrical
include("joints.jl")

export SphericalSpherical, UniversalSpherical, JointUSR, JointRRR, PrismaticConstraint
include("fancy_joints.jl")

export RollingWheelJoint, RollingWheel, SlipWheelJoint, SlippingWheel, RollingWheelSet, RollingWheelSetJoint, RollingConstraintVerticalWheel
include("wheels.jl")

export Spring, Damper, SpringDamperParallel, Torque, Force, WorldForce, WorldTorque
include("forces.jl")

export PartialCutForceBaseSensor, BasicCutForce, BasicCutTorque, CutTorque, CutForce, Power
include("sensors.jl")

export point_to_point, traj5, KinematicPTP, Kinematic5
include("robot/path_planning.jl")
include("robot/robot_components.jl")
include("robot/FullRobot.jl")


export PlanarMechanics
include("PlanarMechanics/PlanarMechanics.jl")


# These are extended in the render module
function get_rot_fun end
function get_fun end
function get_frame_fun end
function get_color end
function get_shapefile end
function get_shape end

export SphereVisualizer, CylinderVisualizer, BoxVisualizer
include("visualizers.jl")

end
