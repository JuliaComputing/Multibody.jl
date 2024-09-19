"""
    Revolute(; name, phi = 0.0, tau = 0.0, axisflange = false)
A revolute joint

# Parameters:
- `axisflange=false`: If `true`, a force flange is enabled, otherwise implicitly grounded"
- `phi`: [rad] Initial angular position for the flange
- `tau`: [Nm] Initial Cut torque in the flange

# Variables:
- `phi(t)`: [rad] angular position
- `w(t)`: [rad/s] angular velocity
- `α(t)`: [rad/s²] angular acceleration
- `tau(t)`: [Nm] torque

# Connectors
- `frame_a` [Frame](@ref)
- `frame_b` [Frame](@ref)
- `fixed` [Fixed](@ref) if `axisflange == false`
- `flange_a` [Flange](@ref) if `axisflange == true`
- `support` [Support](@ref) if `axisflange == true`
"""
@component function Revolute(;
        name,
        axisflange = false, render = true, radius = 0.1, color = [1.0, 0.0, 0.0, 1.0], phi=nothing, w=nothing,
        iscut = false,
        state_priority = 10)
    @named partial_frames = PartialTwoFrames()
    @unpack frame_a, frame_b = partial_frames
    systems = [frame_a, frame_b]

    vars = @variables begin
        (phi(t) = phi), [state_priority=state_priority]
        (w(t) = w), [state_priority=state_priority]
        α(t), [state_priority=state_priority]
        tau(t)
    end

    pars = @parameters begin
        render = render, [description = "Render the joint in animations"]
        radius = radius, [description = "Radius of the body in animations"]
        color[1:4] = color, [description = "Color of the body in animations"]
    end

    eqs = [
        w ~ D(phi),
        α ~ D(w),
        # rigidly connect positions
        frame_a.x ~ frame_b.x,
        frame_a.y ~ frame_b.y,
        # balance forces
        frame_a.fx + frame_b.fx ~ 0,
        frame_a.fy + frame_b.fy ~ 0,
        # balance torques
        frame_a.tau + frame_b.tau ~ 0,
        frame_a.tau ~ tau
    ]
    if !iscut
        push!(eqs, frame_a.phi + phi ~ frame_b.phi)
    end

    if axisflange
        @named fixed = Rotational.Fixed()
        @named flange_a = Rotational.Flange(; phi, tau)
        @named support = Rotational.Support()
        push!(systems, fixed)
        push!(systems, flange_a)
        push!(systems, support)
        append!(eqs, [
            connect(fixed.flange, support)
            flange_a.phi ~ phi
            flange_a.tau ~ tau
            ])
    else
        # actutation torque
        push!(eqs, tau ~ 0)
    end


    return compose(ODESystem(eqs, t, vars, pars; name = name),
        systems...)
end

"""
    Prismatic(; name, f, s = 0, axisflange = false)
A prismatic joint

# Parameters
- `r`: [m, m] x,y-direction of the rod wrt. body system at phi=0
- `axisflange=false`: If `true`, a force flange is enabled, otherwise implicitly grounded"
- `render`: Render the joint in animations
- `radius`: Radius of the body in animations
- `color`: Color of the body in animations

# Variables
  - `s(t)`: [m] Elongation of the joint
  - `v(t)`: [m/s] Velocity of elongation
  - `a(t)`: [m/s²] Acceleration of elongation
  - `f(t)`: [N] Force in direction of elongation

# Connectors
  - `frame_a` [Frame](@ref)
  - `frame_b` [Frame](@ref)
  - `fixed` [Fixed](@ref) if `axisflange == false`
  - `flange_a` [Flange](@ref) if `axisflange == true`
  - `support` [Support](@ref) if `axisflange == true`
"""
@component function Prismatic(;
        name,
        r = [0,0],
        s = nothing,
        v = nothing,
        axisflange = false,
        render = true,
        radius = 0.1,
        color = [0,0.8,1,1],
        state_priority = 2,
        iscut = false,
    )
    @named partial_frames = PartialTwoFrames()
    @unpack frame_a, frame_b = partial_frames

    systems = @named begin
        fixed = TranslationalModelica.Support()
    end

    if axisflange
        more_systems = @named begin
            flange_a = TranslationalModelica.Flange()
            support = TranslationalModelica.Flange()
        end
        systems = [systems; more_systems]
    end

    pars = @parameters begin
        (r[1:2] = r), [description="Direction of the rod wrt. body system at phi=0"]
        render = render, [description="Render the joint in animations"]
        radius = radius, [description="Radius of the body in animations"]
        color[1:4] = color, [description="Color of the body in animations"]
    end

    vars = @variables begin
        (s(t) = s), [state_priority = state_priority, description="Joint coordinate"]
        (v(t) = v), [state_priority = state_priority]
        a(t)
        f(t)
        e0(t)[1:2]
        (r0(t)[1:2]=r), [description="Translation vector of the prismatic rod resolved w.r.t. inertial frame"]
    end

    e = Multibody._normalize(r)
    R = ori_2d(frame_a)
    
    eqs = [
        e0 .~ R * e
        r0 .~ e0 * s
        v ~ D(s)
        a ~ D(v)
        # rigidly connect positions
        frame_a.x + r0[1] ~ frame_b.x
        frame_a.y + r0[2] ~ frame_b.y
        frame_a.fx + frame_b.fx ~ 0
        frame_a.fy + frame_b.fy ~ 0
        frame_a.tau + frame_b.tau + r0[1] * frame_b.fy - r0[2] * frame_b.fx ~ 0
        e0[1] * frame_a.fx + e0[2] * frame_a.fy ~ f
    ]
    if !iscut
        push!(eqs, frame_a.phi ~ frame_b.phi)
    end

    if axisflange
        append!(eqs, [
            connect(fixed, support)
            flange_a.s ~ s
            flange_a.f ~ f
            ])
    else
        # actutation torque
        push!(eqs, f ~ 0)
    end
    return extend(ODESystem(eqs, t, vars, pars; name, systems), partial_frames)
end
