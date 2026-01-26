purple = Multibody.purple
"""
    Fixed(; name, r = [0.0, 0.0], phi = 0.0)

Frame fixed in the planar world frame at a given position and orientation

# Parameters:
- `r`: [m, m] Fixed absolute x,y-position, resolved in world frame
- `phi`: [rad] Fixed angle

# Connectors:
- `frame_b`: 2-dim. Coordinate system
"""
@component function Fixed(; name, r = [0.0, 0.0], phi = 0.0)
    pars = @parameters begin
        r[1:2] = r, [description = "Fixed absolute xy-position, resolved in planarWorld frame"]
        phi = phi, [description = "Fixed angle"]
    end

    systems = @named begin
        frame_b = Frame()
    end

    vars = @variables begin
    end

    equations = Equation[
        [frame_b.x
        frame_b.y] ~ r
        frame_b.phi ~ phi
    ]

    return System(equations, t; name, systems)
end

"""
    Body(; name, m=1, I=0.1, r=0, gy=-9.80665, radius=0.1, render=true, color=Multibody.purple)

Body component with mass and inertia

# Parameters:
- `m`: [kg] mass of the body
- `I`: [kg.m²] inertia of the body with respect to the origin of `frame` along the z-axis of `frame`
- `r`: [m, m] Translational position x,y-position
- `gy`: [m/s²] gravity field acting on the mass in the y-direction, positive value acts in the positive direction defaults to -9.80665
- `radius`: [m] Radius of the body in animations
- `render`: [Bool] Render the body in animations
- `color`: [Array{Float64,1}] Color of the body in animations

# Variables:
  - `r`: [m, m] x,y position
  - `v`: [m/s, m/s] x,y velocity
  - `a`: [m/s², m/s²] x,y acceleration
  - `phi`: [rad] rotation angle (counterclockwise)
  - `w`: [rad/s] angular velocity
  - `α`: [rad/s²] angular acceleration

# Connectors:
  - `frame`: 2-dim. Coordinate system
"""
@component function Body(; name, m, I, r = nothing, v=nothing, phi = nothing, w=nothing, gy = -9.80665, radius=0.1, render=true, color=Multibody.purple, state_priority=2)
    systems = @named begin
        frame_a = Frame()
    end
    pars = @parameters begin
        m = m, [description = "Mass of the body"]
        I = I, [description = "Inertia of the body with respect to the origin of frame_a along the z-axis of frame_a"]
        gy = gy, [description = "Gravity field acting on the mass in the y-direction, positive value acts in the positive direction"]
        radius = radius, [description = "Radius of the body in animations"]
        render = render, [description = "Render the body in animations"]
        color[1:4] = color, [description = "Color of the body in animations"]
        z = 0, [description = "Fixed z-position"]
    end

    vars = @variables begin
        f(t)[1:2], [description = "Force"]
        (r(t)[1:2] = r), [state_priority=state_priority, description = "x,y position"]
        (v(t)[1:2] = v), [state_priority=state_priority, description = "x,y velocity"]
        a(t)[1:2], [description = "x,y acceleration"]
        (phi(t) = phi), [state_priority=state_priority, description = "Rotation angle"]
        (w(t) = w), [state_priority=state_priority, description = "Angular velocity"]
        α(t), [description = "Angular acceleration"]
    end

    eqs = [
        # velocity is the time derivative of position
        r ~ [frame_a.x, frame_a.y]
        v ~ D(r)
        phi ~ frame_a.phi
        w ~ D(phi)
        # acceleration is the time derivative of velocity
        a ~ D(v)
        α ~ D(w)
        # newton's law
        f ~ [frame_a.fx, frame_a.fy]
        f + [0, m*gy] ~ m*a#ifelse(gy !== nothing, fy / m + gy, fy / m),
        I * α ~ frame_a.tau
    ]

    return System(eqs, t, vars, pars; name, systems)
end

"""
    BodyShape(; name, r = [1,0], r_cm = 0.5*r, gy = -9.80665)

The `BodyShape` component is similar to a [`Body`](@ref), but it has two frames and a vector `r` that describes the translation between them, while the body has a single frame only.

# Parameters
- `r`: (Structural) Vector from `frame_a` to `frame_b` resolved in `frame_a`
- `r_cm`: (Structural) Vector from `frame_a` to the center of mass resolved in `frame_a`

# Subsystems
- `translation`: [FixedTranslation](@ref) Fixed translation between `frame_a` and `frame_b`
- `translation_cm`: [FixedTranslation](@ref) Fixed translation between `frame_a` and the center of mass
- `body`: [Body](@ref) Body component placed at center of mass. This component holds the inertial properties

# Connectors
- `frame_a`
- `frame_b`
"""
@component function BodyShape(; name, r0=nothing, r = [1, 0], r_cm = 0.5*r, state_priority = 2,
                                gy = -9.80665, m = 1, I = 0.1, radius = 0.1,
                                render = true, color = purple, v=nothing, phi=nothing, w=nothing)
    pars = @parameters begin
        r0[1:2] = r0, [description = "Vector from frame_a to frame_b resolved in frame_a"]
        r[1:2] = r, [description = "Vector from frame_a to frame_b resolved in frame_a"]
        r_cm[1:2] = r_cm, [description = "Vector from frame_a to center of mass, resolved in frame_a"]
        gy = gy, [description = "Gravity field acting on the mass in the y-direction"]
        m = m, [description = "Mass of the body"]
        I = I, [description = "Inertia of the body with respect to the center of mass"]
        radius = radius, [description = "Radius of the body in animations"]
        render = render, [description = "Render the body in animations"]
        color[1:4] = color, [description = "Color of the body in animations"]
    end

    systems = @named begin
        translation = FixedTranslation(; r, render=false)
        translation_cm = FixedTranslation(; r=r_cm, render=false)
        body = Body(; r=r0, I, m, gy, state_priority, v, phi, w)
        frame_a = Frame()
        frame_b = Frame()
    end

    vars = @variables begin
    end

    equations = Equation[
        connect(frame_a, translation.frame_a, translation_cm.frame_a)
        connect(frame_b, translation.frame_b)
        connect(translation_cm.frame_b, body.frame_a)
    ]

    return System(equations, t, vars, pars; name, systems)
end

"""
    FixedTranslation(; name, r::AbstractArray, l)

A fixed translation between two components (rigid rod)

# Parameters:
- `rx`: [m] Fixed x-length of the rod resolved w.r.t to body frame_a at phi = 0
- `ry`: [m] Fixed y-length of the rod resolved w.r.t to body frame_a at phi = 0
- `radius`: [m] Radius of the rod in animations
- `render`: [Bool] Render the rod in animations

# Connectors:

- `frame_a` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
- `frame_b` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque

"""
@component function FixedTranslation(; name, r = [1.0, 0], radius = 0.1, render = true)
    pars = @parameters begin
        r[1:2] = r, [description = "Fixed x,y-length of the rod resolved w.r.t to body frame_a at phi = 0"]
        radius = radius, [description = "Radius of the rod in animations"]
        render = render, [description = "Render the rod in animations"]
    end

    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
    end

    vars = @variables begin
        phi(t), [state_priority=1, description = "Angle"]
        w(t), [state_priority=1, description = "Angular velocity"]
    end

    # Calculations from begin blocks
    # r = collect(r)
    R = ori_2d(phi)
    r0 = R * r

    equations = Equation[
        phi ~ frame_a.phi
        w ~ D(phi)
        # rigidly connect positions
        [frame_a.x
         frame_a.y] + r0 ~ [frame_b.x, frame_b.y]
        frame_a.phi ~ frame_b.phi
        # balancing force including lever principle
        frame_a.fx + frame_b.fx ~ 0
        frame_a.fy + frame_b.fy ~ 0
        frame_a.tau + frame_b.tau + dot(r0, [frame_b.fy, -frame_b.fx]) ~ 0
    ]

    return System(equations, t, vars, pars; name, systems)
end

@component function FixedRotation(; name, alpha = 0, render = true)
    pars = @parameters begin
        alpha = alpha, [description = "Fixed rotation angle between frame_a and frame_b"]
        render = render, [description = "Render the rotation in animations"]
    end

    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
    end

    vars = @variables begin
    end

    equations = Equation[
        # Positions are the same
        frame_a.x ~ frame_b.x
        frame_a.y ~ frame_b.y

        # Fixed rotation offset
        frame_a.phi + alpha ~ frame_b.phi

        # Force balance
        frame_a.fx + frame_b.fx ~ 0
        frame_a.fy + frame_b.fy ~ 0
        frame_a.tau + frame_b.tau ~ 0
    ]

    return System(equations, t; name, systems)
end

"""
    Spring(; name, c_x = 1, c_y = 1, c_phi = 1e5, s_relx0 = 0, s_rely0 = 0, phi_rel0 = 0, s_small = 1.e-10)

Linear 2D translational spring

# Parameters:
- `c_x`: [N/m] Spring constant in x dir
- `c_y`: [N/m] Spring constant in y dir
- `c_phi`: [N.m/rad] Spring constant in phi dir
- `s_relx0`: [m] Unstretched spring length
- `s_rely0`: [m] Unstretched spring length
- `phi_rel0`: [rad] Unstretched spring angle
- `s_small`: [m] Prevent zero-division if distance between frame_a and frame_b is zero
- `num_windings`: [Int] Number of windings of the coil when rendered
- `color = [0,0,1,1]` Color of the spring in animations
- `render = true` Render the spring in animations
- `radius = 0.1` Radius of spring when rendered
- `N = 200` Number of points in mesh when rendered

# Connectors:
- `frame_a` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
- `frame_b` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
"""
@component function Spring(; name, c_x = 1, c_y = 1, c_phi = 1.0e5, s_relx0 = 0,
                            s_rely0 = 0, phi_rel0 = 0, s_small = 1.0e-10,
                            num_windings = 6, color = [0, 0.0, 1, 1],
                            render = true, radius = 0.1, N = 200)
    pars = @parameters begin
        c_x = c_x, [description = "Spring constant in x dir"]
        c_y = c_y, [description = "Spring constant in y dir"]
        c_phi = c_phi, [description = "Spring constant"]
        s_relx0 = s_relx0, [description = "Unstretched spring length"]
        s_rely0 = s_rely0, [description = "Unstretched spring length"]
        phi_rel0 = phi_rel0, [description = "Unstretched spring angle"]
        s_small = s_small, [description = "Prevent zero-division if distance between frame_a and frame_b is zero"]
        num_windings = num_windings, [description = "Number of windings of the coil when rendered"]
        color[1:4] = color, [description = "Color of the spring in animations"]
        render = render, [description = "Render the spring in animations"]
        radius = radius, [description = "Radius of spring when rendered"]
        N = N, [description = "Number of points in mesh when rendered"]
    end

    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
    end

    vars = @variables begin
        s_relx(t), [description = "Relative x position"]
        s_rely(t), [description = "Relative y position"]
        phi_rel(t), [description = "Relative angle"]
        f_x(t), [description = "Force in x direction"]
        f_y(t), [description = "Force in y direction"]
    end

    # Calculations from begin block
    r_rel_0 = [s_relx, s_rely, 0]
    l = sqrt(r_rel_0' * r_rel_0)
    e_rel_0 = r_rel_0 / max(l, s_small)

    equations = Equation[
        phi_rel ~ frame_b.phi - frame_a.phi
        frame_a.tau ~ 0
        frame_b.tau ~ 0
        s_relx ~ frame_b.x - frame_a.x
        s_rely ~ frame_b.y - frame_a.y
        f_x ~ c_x * (s_relx - s_relx0)
        f_y ~ c_y * (s_rely - s_rely0)
        frame_a.fx ~ -f_x
        frame_b.fx ~ f_x
        frame_a.fy ~ -f_y
        frame_b.fy ~ f_y
    ]

    return System(equations, t; name, systems)
end

"""
    Damper(; name, d = 1, s_small = 1.e-10)

Linear (velocity dependent) damper

# Parameters:
- `d`: [N.s/m] Damping constant
- `s_small`: [m] Prevent zero-division if distance between frame_a and frame_b is zero


# Connectors:
- `frame_a` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
- `frame_b` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
"""
@component function Damper(; name, d = 1, s_small = 1.0e-10)
    pars = @parameters begin
        d = d, [description = "Damping constant"]
        s_small = s_small, [description = "Prevent zero-division if distance between frame_a and frame_b is zero"]
    end

    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
    end

    vars = @variables begin
        r0x(t), [description = "Relative position x component"]
        r0y(t), [description = "Relative position y component"]
        d0x(t), [description = "Direction x component"]
        d0y(t), [description = "Direction y component"]
        vx(t), [description = "Relative velocity x component"]
        vy(t), [description = "Relative velocity y component"]
        v(t), [description = "Relative velocity magnitude"]
        f(t), [description = "Damping force"]
    end

    # Calculations from begin block
    r0 = [r0x, r0y]
    l = sqrt(r0' * r0)

    equations = Equation[
        frame_a.x + r0x ~ frame_b.x
        frame_a.y + r0y ~ frame_b.y
        D(frame_a.x) + vx ~ D(frame_b.x)
        D(frame_a.y) + vy ~ D(frame_b.y)
        v ~ [vx, vy]' * [d0x, d0y]
        f ~ -d * v
        d0x ~ ifelse(l < s_small, r0[1], r0[1] / l)
        d0y ~ ifelse(l < s_small, r0[2], r0[2] / l)
        frame_a.fx ~ d0x * f
        frame_a.fy ~ d0y * f
        frame_a.tau ~ 0
        frame_a.fx + frame_b.fx ~ 0
        frame_a.fy + frame_b.fy ~ 0
        frame_a.tau + frame_b.tau ~ 0

        # lossPower ~ -f * v
    ]

    return System(equations, t; name, systems)
end

"""
    SpringDamper(; name, c_x = 1, c_y = 1, c_phi = 1e5, d_x = 1, d_y = 1, d_phi = 1, s_relx0 = 0, s_rely0 = 0, phi_rel0 = 0, s_small = 1.e-10)

Linear 2D translational spring damper model

# Parameters:
- `c_x`: [N/m] Spring constant in x dir
- `c_y`: [N/m] Spring constant in y dir
- `c_phi`: [N.m/rad] Spring constant in phi dir
- `d_x`: [N.s/m] Damping constant in x dir
- `d_y`: [N.s/m] Damping constant in y dir
- `d_phi`: [N.m.s/rad] Damping constant in phi dir
- `s_relx0`: [m] Unstretched spring length
- `s_rely0`: [m] Unstretched spring length
- `phi_rel0`: [rad] Unstretched spring angle
- `s_small`: [m] Prevent zero-division if distance between frame_a and frame_b is zero
- `num_windings`: [Int] Number of windings of the coil when rendered
- `color = [0,0,1,1]` Color of the spring in animations
- `render = true` Render the spring in animations
- `radius = 0.1` Radius of spring when rendered
- `N = 200` Number of points in mesh when rendered

# Connectors:
- `frame_a` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
- `frame_b` [Frame](@ref) Coordinate system fixed to the component with one cut-force and cut-torque
"""
@component function SpringDamper(; name, c_x = 1, c_y = 1, c_phi = 1.0e5, d_x = 1, d_y = 1, d_phi = 1,
                                  s_relx0 = 0, s_rely0 = 0, phi_rel0 = 0, s_small = 1.0e-10,
                                  num_windings = 6, color = [0, 0.0, 1, 1],
                                  render = true, radius = 0.1, N = 200)
    pars = @parameters begin
        c_x = c_x, [description = "Spring constant in x dir"]
        c_y = c_y, [description = "Spring constant in y dir"]
        c_phi = c_phi, [description = "Spring constant in phi dir"]
        d_x = d_x, [description = "Damping constant in x dir"]
        d_y = d_y, [description = "Damping constant in y dir"]
        d_phi = d_phi, [description = "Damping constant in phi dir"]
        s_relx0 = s_relx0, [description = "Unstretched spring length"]
        s_rely0 = s_rely0, [description = "Unstretched spring length"]
        phi_rel0 = phi_rel0, [description = "Unstretched spring angle"]
        s_small = s_small, [description = "Prevent zero-division if distance between frame_a and frame_b is zero"]
        num_windings = num_windings, [description = "Number of windings of the coil when rendered"]
        color[1:4] = color, [description = "Color of the spring in animations"]
        render = render, [description = "Render the spring in animations"]
        radius = radius, [description = "Radius of spring when rendered"]
        N = N, [description = "Number of points in mesh when rendered"]
    end

    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
    end

    vars = @variables begin
        v_relx(t), [description = "Relative velocity in x direction"]
        v_rely(t), [description = "Relative velocity in y direction"]
        w_rel(t), [description = "Relative angular velocity"]
        s_relx(t), [description = "Relative x position"]
        s_rely(t), [description = "Relative y position"]
        phi_rel(t), [description = "Relative angle"]
        f_x(t), [description = "Force in x direction"]
        f_y(t), [description = "Force in y direction"]
        tau(t), [description = "Torque"]
    end

    # Calculations from begin block
    r_rel_0 = [s_relx, s_rely, 0]
    l = sqrt(r_rel_0' * r_rel_0)
    e_rel_0 = r_rel_0 / max(l, s_small)

    equations = Equation[
        s_relx ~ frame_b.x - frame_a.x
        s_rely ~ frame_b.y - frame_a.y
        phi_rel ~ frame_b.phi - frame_a.phi
        v_relx ~ D(s_relx)
        v_rely ~ D(s_rely)
        w_rel ~ D(phi_rel)

        tau ~ c_phi * (phi_rel - phi_rel0) + d_phi * w_rel
        frame_a.tau ~ -tau
        frame_b.tau ~ tau
        f_x ~ c_x * (s_relx - s_relx0) + d_x * v_relx
        f_y ~ c_y * (s_rely - s_rely0) + d_y * v_rely
        frame_a.fx ~ -f_x
        frame_b.fx ~ f_x
        frame_a.fy ~ -f_y
        frame_b.fy ~ f_y

        # lossPower ~ d_x * v_relx * v_relx + d_y * v_rely * v_rely
    ]

    return System(equations, t; name, systems)
end


"""
    SimpleWheel(; name, radius = 0.3, color = [1, 0, 0, 1], μ = 1e9)

Simple wheel model with viscous lateral friction and a driving torque

# Connectors:
- `frame_a` (Frame) Coordinate system fixed to the component with one cut-force and cut-torque
- `thrust` (RealInput) Input for the longitudinal force applied to the wheel

# Parameters:
- `μ`: [Ns/m] Viscous friction coefficient
- `radius`: [m] Radius of the wheel
- `color`: Color of the wheel in animations

# Variables:
- `θ`: [rad] Wheel angle
- `Vx`: [m/s] Longitudinal velocity (resolved in local frame)
- `Vy`: [m/s] Lateral velocity (resolved in local frame)
- `Fy`: [N] Lateral friction force
- `Fx`: [N] Applied longitudinal wheel force
"""
@component function SimpleWheel(; name, friction_model = :viscous, radius = 0.3,
                                 color = [1, 0, 0, 1], μ = 1e9)
    pars = @parameters begin
        radius = radius, [description = "Radius of the wheel"]
        color[1:4] = color, [description = "Color of the wheel in animations"]
        μ = μ, [description = "Viscous friction coefficient"]
        # Fy0 = 1e4, [description = "Lateral friction force at zero longitudinal force"]
        # μx = Fy0, [description = "Maximum longitudinal friction force"]
    end

    systems = @named begin
        frame_a = Frame()
        thrust = Blocks.RealInput()
    end

    vars = @variables begin
        θ(t), [guess=0, description="Wheel angle"]
        Vx(t), [guess=0, description="Longitudinal velocity (resolved in local frame)"]
        Vy(t), [guess=0, description="Lateral velocity (resolved in local frame)"]
        Fy(t), [guess=0, description="Lateral friction force"]
        Fx(t), [guess=0, description="Applied longitudinal wheel force"]
    end

    # Calculations from begin block
    R_W_F = ori_2d(frame_a) # rotation matrix, local to global
    veqs = R_W_F'*[D(frame_a.x), D(frame_a.y)] #~ [Vx, Vy]
    feqs = R_W_F'*[frame_a.fx, frame_a.fy] #~ [Fx, Fy]


    equations = Equation[
        θ ~ frame_a.phi

        # road friction
        Fx ~ thrust.u
        # if friction_model == :viscous
            Fy ~ μ*Vy
        # elseif friction_model == :ellipse
            # Fy ~ Fy0*sqrt(1 - (Fx/μ)^2)*Vy
        # end
        veqs[1] ~ Vx
        veqs[2] ~ Vy
        feqs[1] ~ Fx
        feqs[2] ~ Fy

        # R'*[D(frame.x), D(frame.y)] ~ [Vx, Vy]
        # R'*[frame.fx, frame.fy] ~ [Fx, Fy]

        frame_a.tau ~ 0 # Assume that wheel does not transmit any torque
    ]

    return System(equations, t; name, systems)
end

"""
    limit_S_triple(x_max, x_sat, y_max, y_sat, x)

Returns a point-symmetric Triple S-Function

A point symmetric interpolation between points `(0, 0), (x_max, y_max) and (x_sat, y_sat)`, provided `x_max < x_sat`. The approximation is done in such a way that the function's 1st derivative is zero at points `(x_max, y_max)` and `(x_sat, y_sat)`. Thus, the function's 1st derivative is continuous for all `x`. The higher derivatives are discontinuous at these points.

```
x_max = 0.2
x_sat = 0.5
y_max = 1.4
y_sat = 1.2

plot(x->Multibody.PlanarMechanics.limit_S_triple(x_max, x_sat, y_max, y_sat, x), -1, 1)
vline!([x_max x_sat], label=["x_max" "x_sat"])

            ┌────────────────────────────────────────┐ 
    1.48385 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⢀⡔⠢⠤⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
            │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⡜⠀⠀⠀⠈⠉⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠂⠀│ 
            │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⢠⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
            │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⡜⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
            │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣇⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
            │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
            │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
            │⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⢤⡧⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤│ 
            │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
            │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
            │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⠁⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
            │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡜⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
            │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠃⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
            │⠀⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⣀⡀⠀⠀⠀⡜⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
   -1.48377 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠒⠦⠼⠁⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
            └────────────────────────────────────────┘ 
            ⠀-1.06⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀1.06⠀ 
```
"""
function limit_S_triple(x_max, x_sat, y_max, y_sat, x)
    if x > x_max
        return limit_S_form(x_max, x_sat, y_max, y_sat, x)
    elseif x < -x_max
        return limit_S_form(-x_max, -x_sat, -y_max, -y_sat, x)
    else
        return limit_S_form(-x_max, x_max, -y_max, y_max, x)
    end
end


"""
    limit_S_form(x_min, x_max, y_min, y_max, x)

Returns a S-shaped transition

A smooth transition between points `(x_min, y_min)` and `(x_max, y_max)`. The transition is done in such a way that the function's 1st derivative is continuous for all `x`. The higher derivatives are discontinuous at input points.

```
x_min = -0.4
x_max = 0.6
y_max = 1.4
y_min = 1.2

julia> plot(x->Multibody.PlanarMechanics.limit_S_form(x_min, x_max, y_min, y_max, x), -1, 1, legend=false)
         ┌────────────────────────────────────────┐ 
   1.406 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⡠⠔⠒⠒⠒⠒⠒⠒⠒⠂⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⢠⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⡴⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⡰⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⡜⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⡰⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⢰⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⢠⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣧⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⡏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠃⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠇⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠃⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡔⠁⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
   1.194 │⠀⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠖⠉⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         └────────────────────────────────────────┘ 
         ⠀-1.06⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀1.06⠀ 
```
"""
function limit_S_form(x_min, x_max, y_min, y_max, x)
    x2 = x - x_max/2 - x_min/2
    x3 = x2*2/(x_max-x_min)
    y1 = if x3 > 1
        1
    elseif x3 < -1
        -1
    else
        -0.5*x3^3 + 1.5*x3
    end
    y2 = y1*(y_max-y_min)/2
    y3 = y2 + y_max/2 + y_min/2
    return y3
end


@register_symbolic limit_S_triple(x_max::Real, x_sat::Real, y_max::Real, y_sat::Real, x::Real)
@register_symbolic limit_S_form(x_min::Real, x_max::Real, y_min::Real, y_max::Real, x::Real)

"""
    SlipBasedWheelJoint(;
        name,
        r = [1, 0],
        N,
        vAdhesion_min,
        vSlide_min,
        sAdhesion,
        sSlide,
        mu_A,
        mu_S,
        render = true,
        color = [0.1, 0.1, 0.1, 1],
        z = 0,
        diameter = 0.1,
        width = diameter * 0.6,
        radius = 0.1,
        w_roll = nothing,
    )

Slip-based wheel joint

The ideal wheel joint models the behavior of a wheel rolling on a x,y-plane whose contact patch has slip-dependent friction characteristics. This is an approximation for wheels with a rim and a rubber tire.

The force depends with friction characteristics on the slip. The slip is split into two components:

- lateral slip: the lateral velocity divided by the rolling velocity.
- longitudinal slip: the longitudinal slip velocity divided by the rolling velocity.

For low rolling velocity this definition become ill-conditioned. Hence a dry-friction model is used for low rolling velocities. For **zero rolling velocity**, the intitialization might fail if automatic differentiation is used. Either start with a non-zero (but tiny) rolling velocity or pass `autodiff=false` to the solver.

The radius of the wheel can be specified by the parameter `radius`. The driving direction (for `phi = 0`) can be specified by the parameter `r`. The normal load is set by `N`.

The wheel contains a 2D connector `frame_a` for the steering on the plane. The rolling motion of the wheel can be actuated by the 1D connector `flange_a`.

In addition there is an input `dynamicLoad` for a dynamic component of the normal load.

# Connectors:
- `frame_a` (Frame) Coordinate system fixed to the component with one cut-force and cut-torque
- `flange_a` (Rotational.Flange) Flange for the rolling motion
- `dynamicLoad` (Blocks.RealInput) Input for the dynamic component of the normal load (must be connected)

# Terminology: 
- _Adhesion_ refers to the peak of the traction curve, where the slip is such that the maximum amount of traction is generated.
- _Sliding velocity_ refers to the velocity at which the traction curve saturates and stays constant with increased slip velocity.
"""
@component function SlipBasedWheelJoint(;
    name,
    r = [1, 0],
    N,
    vAdhesion_min,
    vSlide_min,
    sAdhesion,
    sSlide,
    mu_A,
    mu_S,
    render = true,
    color = [0.1, 0.1, 0.1, 1],
    z = 0,
    diameter = 0.1,
    width = diameter * 0.6,
    radius = 0.1,
    phi_roll = nothing,
    w_roll = nothing,
)
    systems = @named begin
        frame_a = Frame()
        flange_a = Rotational.Flange()
        dynamicLoad = Blocks.RealInput()
    end
    pars = @parameters begin
        r[1:2] = r, [description = "Driving direction of the wheel at angle phi = 0"]
        N = N, [description = "Base normal load"]
        vAdhesion_min = vAdhesion_min, [description = "Minimum adhesion velocity"]
        vSlide_min = vSlide_min, [description = "Minimum sliding velocity"]
        sAdhesion = sAdhesion, [description = "Adhesion slippage"]
        sSlide = sSlide, [description = "Sliding slippage"]
        mu_A = mu_A, [description = "Friction coefficient at adhesion"]
        mu_S = mu_S, [description = "Friction coefficient at sliding"]
        render = render, [description = "Render the wheel in animations"]
        color[1:4] = color, [description = "Color of the wheel in animations"]
        z = 0, [description = "Position z of the body"]
        diameter = diameter, [description = "Diameter of the rims"]
        width = width, [description = "Width of the wheel"]
        radius = radius, [description = "Radius of the wheel"]
    end
    r = collect(r)
    l = Multibody._norm(r)
    e = Multibody._normalize(r) # Unit vector in direction of r

    vars = @variables begin
        e0(t)[1:2], [description="Unit vector in direction of r resolved w.r.t. inertial frame"]
        (phi_roll(t) = phi_roll), [guess=0, description="wheel angle"] # wheel angle
        (w_roll(t)=w_roll), [guess=0, description="Roll velocity of wheel"]
        v(t)[1:2], [description="velocity"]
        v_lat(t), [guess=0, description="Velocity in lateral direction"]
        v_long(t), [guess=0, description="Velocity in longitudinal direction"]
        v_slip_long(t), [guess=0, description="Slip velocity in longitudinal direction"]
        v_slip_lat(t), [guess=0, description="Slip velocity in lateral direction"]
        v_slip(t), [description="Slip velocity, norm of component slip velocities"]
        f(t), [description="Total traction force"]
        f_lat(t), [description="Lateral force"]
        f_long(t), [description="Longitudinal force"]
        fN(t), [description="Base normal load"]
        vAdhesion(t), [description="Adhesion velocity"]
        vSlide(t), [description="Sliding velocity"]
    end
    e0 = collect(e0)
    v = collect(v)
    vars = [
        e0; phi_roll; w_roll; v; v_lat; v_long; v_slip_long; v_slip_lat; v_slip; f; f_lat; f_long; fN; vAdhesion; vSlide
    ]

    R = ori_2d(frame_a)

    equations = Equation[
        e0 .~ R * e
        v .~ D.([frame_a.x, frame_a.y])

        phi_roll ~ flange_a.phi
        w_roll ~ D(phi_roll)
        v_long ~ v' * e0
        v_lat ~ -v[1] * e0[2] + v[2] * e0[1]
        v_slip_lat ~ v_lat - 0
        v_slip_long ~ v_long - radius * w_roll
        v_slip ~ sqrt(v_slip_long^2 + v_slip_lat^2) + 0.0001
        -f_long * radius ~ flange_a.tau
        frame_a.tau ~ 0
        vAdhesion ~ max(vAdhesion_min, sAdhesion * abs(radius * w_roll))
        vSlide ~ max(vSlide_min, sSlide * abs(radius * w_roll))
        fN ~ max(0, N + dynamicLoad.u)
        f ~ fN * limit_S_triple(vAdhesion, vSlide, mu_A, mu_S, v_slip) # limit_S_triple(x_max, x_sat, y_max, y_sat, x)
        f_long ~ f * v_slip_long / v_slip
        f_lat ~ f * v_slip_lat / v_slip
        f_long ~ [frame_a.fx, frame_a.fy]'e0
        f_lat ~ [frame_a.fy, -frame_a.fx]'e0
    ]

    return System(equations, t, vars, pars; name, systems)

end


"""
    OneDOFWheelJoint(;
        name,
        vAdhesion_min,
        vSlide_min,
        sAdhesion,
        sSlide,
        mu_A,
        mu_S,
        render = true,
        color = [0.1, 0.1, 0.1, 1],
        z = 0,
        diameter = 0.1,
        width = diameter * 0.6,
        radius = 0.1,
        phi_roll = nothing,
        w_roll = nothing,
    )

A simplified wheel joint that constrains motion to only the x-axis (forward/backward) with slip-dependent friction.

Unlike `SlipBasedWheelJoint` which allows motion in any direction based on the frame orientation, this joint:
- Fixes the driving direction to the global x-axis `[1, 0]`
- Constrains the y-position to equal the wheel radius (ground contact)
- Eliminates lateral velocity and forces

This makes it suitable for simplified 1-DOF models like a planar Segway where the wheel should only move forward/backward.

The slip-dependent friction model is the same as `SlipBasedWheelJoint`:
- For low velocities, a dry-friction model is used
- The normal load is determined automatically from constraint dynamics (mass and gravity of attached body)
- Friction transitions smoothly between adhesion and sliding regimes

# Parameters:
- `radius`: Wheel radius (also determines y-position constraint)
- `mu_A`: Friction coefficient at adhesion (peak traction)
- `mu_S`: Friction coefficient at sliding (saturated traction)
- `sAdhesion`: Adhesion slip ratio threshold
- `sSlide`: Sliding slip ratio threshold
- `vAdhesion_min`: Minimum adhesion velocity (for low-speed stability)
- `vSlide_min`: Minimum sliding velocity (for low-speed stability)

# Connectors:
- `frame_a` (Frame) Coordinate system fixed to the wheel. Attach a body here for mass/inertia.

# Example
```julia
wheelJoint = OneDOFWheelJoint(
    radius = 0.025,
    mu_A = 1,
    mu_S = 0.7,
    sAdhesion = 0.04,
    sSlide = 0.12,
    vAdhesion_min = 0.05,
    vSlide_min = 0.15,
)
```
"""
@component function OneDOFWheelJoint(;
    name,
    vAdhesion_min,
    vSlide_min,
    sAdhesion,
    sSlide,
    mu_A,
    mu_S,
    render = true,
    color = [0.1, 0.1, 0.1, 1],
    x = nothing,
    v = nothing,
    z = 0,
    diameter = 0.1,
    width = diameter * 0.6,
    radius = 0.1,
    phi_roll = nothing,
    w_roll = nothing,
)
    systems = @named begin
        frame_a = Frame()
    end
    pars = @parameters begin
        vAdhesion_min = vAdhesion_min, [description = "Minimum adhesion velocity"]
        vSlide_min = vSlide_min, [description = "Minimum sliding velocity"]
        sAdhesion = sAdhesion, [description = "Adhesion slippage"]
        sSlide = sSlide, [description = "Sliding slippage"]
        mu_A = mu_A, [description = "Friction coefficient at adhesion"]
        mu_S = mu_S, [description = "Friction coefficient at sliding"]
        render = render, [description = "Render the wheel in animations"]
        color[1:4] = color, [description = "Color of the wheel in animations"]
        z = z, [description = "Position z of the body"]
        diameter = diameter, [description = "Diameter of the rims"]
        width = width, [description = "Width of the wheel"]
        radius = radius, [description = "Radius of the wheel"]
    end

    vars = @variables begin
        (phi_roll(t) = phi_roll), [guess=0, description="Wheel rolling angle"]
        (w_roll(t) = w_roll), [guess=0, description="Wheel rolling velocity"]
        x(t)=x, [description = "Position in x direction"]
        v(t)=v, [guess=0, description="Velocity in longitudinal (x) direction"]
        v_slip_long(t), [guess=0, description="Slip velocity in longitudinal direction"]
        v_slip(t), [description="Slip velocity magnitude"]
        f(t), [description="Total traction force magnitude"]
        f_long(t), [description="Longitudinal friction force"]
        f_n(t), [guess=10.0, description="Normal constraint force (determined by dynamics)"]
        vAdhesion(t), [description="Adhesion velocity threshold"]
        vSlide(t), [description="Sliding velocity threshold"]
    end

    equations = Equation[
        # Velocity in x-direction (fixed driving direction along global x-axis)
        x ~ frame_a.x
        v ~ D(x)

        # Wheel angle coupling
        phi_roll ~ -frame_a.phi
        w_roll ~ D(phi_roll)

        # Longitudinal slip (difference between ground velocity and wheel surface velocity)
        v_slip_long ~ v - radius * w_roll
        v_slip ~ abs(v_slip_long) + 0.0001

        # Slip-dependent friction (like 3D SlipWheelJoint)
        # f_n is the normal force, determined by constraint dynamics when a body is attached
        vAdhesion ~ max(vAdhesion_min, sAdhesion * abs(radius * w_roll))
        vSlide ~ max(vSlide_min, sSlide * abs(radius * w_roll))
        f ~ f_n * limit_S_triple(vAdhesion, vSlide, mu_A, mu_S, v_slip)
        f_long ~ - f * v_slip_long / v_slip

        # Frame forces from contact
        frame_a.fx ~ f_long
        frame_a.fy ~ f_n
        frame_a.tau ~ radius * f_long

        # Position constraint
        frame_a.y ~ radius      # Wheel center at ground level + radius
    ]

    return System(equations, t, vars, pars; name, systems)
end


"""
    RollingWheelJoint(;
        name,
        radius = 0.1,
        render = true,
        color = [0.1, 0.1, 0.1, 1],
        z = 0,
        diameter = 0.1,
        width = diameter * 0.6,
        phi_roll = nothing,
        w_roll = nothing,
    )

An ideal rolling wheel joint constrained to move only along the x-axis (forward/backward) without slip.

Unlike `OneDOFWheelJoint` which uses slip-dependent friction, this joint enforces a perfect rolling constraint:
`D(x) = radius * w_roll` (ground velocity equals wheel surface velocity).

This is suitable for simplified models where:
- Tire slip can be neglected
- The kinematic relationship between wheel rotation and translation is exact
- No friction parameters need to be tuned

# Parameters:
- `radius`: Wheel radius (determines y-position constraint and rolling kinematics)

# Connectors:
- `frame_a` (Frame) Coordinate system fixed to the wheel. Attach a revolute joint and a body inertia here.

# Example
```julia
wheelJoint = RollingWheelJoint(radius = 0.025)
```
"""
@component function RollingWheelJoint(;
    name,
    radius = 0.1,
    render = true,
    color = [0.1, 0.1, 0.1, 1],
    z = 0,
    x = nothing,
    x0 = 0,
    v = nothing,
    diameter = 0.1,
    width = diameter * 0.6,
    phi_roll = nothing,
    w_roll = nothing,
)
    systems = @named begin
        frame_a = Frame()
        # flange_a = Rotational.Flange()
    end
    pars = @parameters begin
        render = render, [description = "Render the wheel in animations"]
        color[1:4] = color, [description = "Color of the wheel in animations"]
        z = z, [description = "Position z of the body"]
        x0 = x0, [description = "x position at zero roll angle"]
        diameter = diameter, [description = "Diameter of the rims"]
        width = width, [description = "Width of the wheel"]
        radius = radius, [description = "Radius of the wheel"]
    end

    vars = @variables begin
        (x(t) = x), [description = "Position x of the body"]
        (v(t) = v), [description = "Velocity x of the body"]
        # (a(t))
        (phi_roll(t) = phi_roll), [guess=0, description="Wheel rolling angle"]
        (w_roll(t) = w_roll), [guess=0, description="Wheel rolling velocity"]
    end

    equations = Equation[
        # Wheel angle coupling to flange
        phi_roll ~ frame_a.phi
        w_roll ~ D(phi_roll)

        x ~ frame_a.x
        v ~ D(x)
        # a ~ D(v)
        # Ideal rolling constraint: ground velocity = wheel surface velocity
        x ~ radius * phi_roll + x0
        # v ~ radius * w_roll

        # Force/torque relationship (from constraint)
        frame_a.fx * radius ~ frame_a.tau

        # 1-DOF constraints
        frame_a.y ~ radius      # Wheel center at ground level + radius
        # frame_a.tau + flange_a.tau ~ 0
        # frame_a.phi ~ flange_a.phi
        # frame_a.fy  ~ 0          # Vertical force handled externally
    ]

    return System(equations, t, vars, pars; name, systems)
end


"""
    IdealPlanetary(; name, ratio = 2)

The IdealPlanetary gear box is an ideal gear without inertia, elasticity, damping or backlash consisting of an inner sun wheel, an outer ring wheel and a planet wheel located between sun and ring wheel. The bearing of the planet wheel shaft is fixed in the planet carrier. The component can be connected to other elements at the sun, ring and/or carrier flanges. It is not possible to connect to the planet wheel. If inertia shall not be neglected, the sun, ring and carrier inertias can be easily added by attaching inertias (= model Inertia) to the corresponding connectors. The inertias of the planet wheels are always neglected.

The ideal planetary gearbox is uniquely defined by the ratio of the number of ring teeth ``z_r`` with respect to the number of sun teeth ``z_s``. For example, if there are 100 ring teeth and 50 sun teeth then ratio = ``z_r/z_s = 2``. The number of planet teeth ``z_p`` has to fulfill the following relationship:
```math
z_p = (z_r - z_s) / 2
```
Therefore, in the above example ``z_p = 25`` is required.

According to the overall convention, the positive direction of all vectors, especially the absolute angular velocities and cut-torques in the flanges, are along the axis vector displayed in the icon.

# Parameters:
- `ratio`: Number of ring teeth/sun teeth

# Connectors:
- `sun` (Rotational.Flange) Sun wheel
- `carrier` (Rotational.Flange) Planet carrier
- `ring` (Rotational.Flange) Ring wheel
"""
@component function IdealPlanetary(; name, ratio = 2)
    pars = @parameters begin
        ratio = ratio, [description = "Number of ring_teeth/sun_teeth"]
    end

    systems = @named begin
        sun = Rotational.Flange()
        carrier = Rotational.Flange()
        ring = Rotational.Flange()
    end

    vars = @variables begin
    end

    equations = Equation[
        (1 + ratio)*carrier.phi ~ sun.phi + ratio*ring.phi
        ring.tau ~ ratio*sun.tau
        carrier.tau ~ -(1 + ratio)*sun.tau
    ]

    return System(equations, t; name, systems)
end

"""
    DifferentialGear(; name)

A 1D-rotational component that is a variant of a planetary gear and can be used to distribute the torque equally among the wheels on one axis.

# Connectors:
- `flange_b` (Rotational.Flange) Flange for the input torque
- `flange_left` (Rotational.Flange) Flange for the left output torque
- `flange_right` (Rotational.Flange) Flange for the right output torque
"""
@component function DifferentialGear(; name)
    pars = @parameters begin
    end

    systems = @named begin
        ideal_planetary = IdealPlanetary(ratio=-2)
        flange_b = Rotational.Flange()
        flange_left = Rotational.Flange()
        flange_right = Rotational.Flange()
    end

    vars = @variables begin
    end

    equations = Equation[
        connect(flange_b, ideal_planetary.ring)
        connect(ideal_planetary.carrier, flange_right)
        connect(ideal_planetary.sun, flange_left)
    ]

    return System(equations, t; name, systems)
end
