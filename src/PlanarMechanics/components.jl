purple = Multibody.purple
"""
    Fixed(; name, r = (0.0, 0.0), phi = 0.0)

Frame fixed in the planar world frame at a given position and orientation

# Parameters:

  - `x`: [m] Fixed absolute x-position, resolved in planarWorld frame
  - `y`: [m] Fixed absolute y-position, resolved in planarWorld frame
  - `phi`: [rad] Fixed angle

# Connectors:

  - `frame: 2-dim. Coordinate system

"""
@mtkmodel Fixed begin
    @parameters begin
        x = 0, [description = "Fixed absolute x-position, resolved in planarWorld frame"]
        y = 0, [description = "Fixed absolute y-position, resolved in planarWorld frame"]
        phi = 0, [description = "Fixed angle"]
    end

    @components begin
        frame = Frame()
    end

    @equations begin
        frame.x ~ x
        frame.y ~ y
        frame.phi ~ phi
    end
end

"""
    Body(; name, m=1, I=0.1, r=0, gy=-9.807, radius=0.1, render=true, color=Multibody.purple)

Body component with mass and inertia

# Parameters:
- `m`: [kg] mass of the body
- `I`: [kg.m²] inertia of the body with respect to the origin of `frame` along the z-axis of `frame`
- `r`: [m, m] Translational position x,y-position
- `gy`: [m/s²] gravity field acting on the mass in the y-direction, positive value acts in the positive direction defaults to -9.807
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
@component function Body(; name, m, I, r = zeros(2), phi = 0, gy = -9.807, radius=0.1, render=true, color=Multibody.purple)
    @named frame_a = Frame()
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
        (r(t)[1:2] = r), [description = "x,y position"]
        v(t)[1:2], [description = "x,y velocity"]
        a(t)[1:2], [description = "x,y acceleration"]
        (phi(t) = phi), [description = "Rotation angle"]
        w(t), [description = "Angular velocity"]
        α(t), [description = "Angular acceleration"]
    end

    eqs = [
        # velocity is the time derivative of position
        r .~ [frame_a.x, frame_a.y]
        v .~ D.(r)
        phi ~ frame_a.phi
        w .~ D.(phi)
        # acceleration is the time derivative of velocity
        a .~ D.(v)
        α .~ D.(w)
        # newton's law
        f .~ [frame_a.fx, frame_a.fy]
        f + [0, m*gy] .~ m*a#ifelse(gy !== nothing, fy / m + gy, fy / m),
        I * α ~ frame_a.tau
    ]

    return compose(ODESystem(eqs, t, vars, pars; name),
        frame_a)
end

"""
    BodyShape(; name, r = [1,0], r_cm = 0.5*r, gy = -9.807)

The `BodyShape` component is similar to a [`Body`](@ref), but it has two frames and a vector `r` that describes the translation between them, while the body has a single frame only.

# Parameters
- `r`: (Structural) Vector from `frame_a` to `frame_b` resolved in `frame_a`
- `r_cm`: (Structural) Vector from `frame_a` to the center of mass resolved in `frame_a`
"""
@mtkmodel BodyShape begin
    @structural_parameters begin
        r = [1,0]
        r_cm = 0.5*r
        gy = -9.807
    end
    @parameters begin
        # r[1:2] = [1,0], [description = "Fixed x,y-length of the rod resolved w.r.t to body frame_a at phi = 0"]
        # r_cm[1:2] = 0.5*r, [description = "Vector from frame_a to center of mass, resolved in frame_a"]
        m = 1, [description = "mass of the body"]
        I = 0.1, [description = "inertia of the body with respect to the center of mass"]
        radius = 0.1, [description = "Radius of the body in animations"]
        render = true, [description = "Render the body in animations"]
        (color[1:4] = purple), [description = "Color of the body in animations"]
    end
    @components begin
        translation = FixedTranslation(; r)
        translation_cm = FixedTranslation(; r=r_cm)
        body = Body(; r=r_cm, I, m, gy)
        frame_a = Frame()
        frame_b = Frame()
    end
    @equations begin
        connect(frame_a, translation.frame_a, translation_cm.frame_a)
        connect(frame_b, translation.frame_b)
        connect(translation_cm.frame_b, body.frame_a)
    end
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
@mtkmodel FixedTranslation begin
    @extend frame_a, frame_b = partial_frames = PartialTwoFrames()

    @parameters begin
        r[1:2] = [1.0, 0],
        [
            description = "Fixed x,y-length of the rod resolved w.r.t to body frame_a at phi = 0"
        ]
        radius = 0.1, [description = "Radius of the rod in animations"]
        render = true, [description = "Render the rod in animations"]
    end
    begin
        r = collect(r)
    end

    begin
        R = [cos(frame_a.phi) -sin(frame_a.phi);
             sin(frame_a.phi) cos(frame_a.phi)]
        r0 = R * r
    end

    @equations begin
        # rigidly connect positions
        frame_a.x + r0[1] ~ frame_b.x
        frame_a.y + r0[2] ~ frame_b.y
        frame_a.phi ~ frame_b.phi
        # balancing force including lever principle
        frame_a.fx + frame_b.fx ~ 0
        frame_a.fy + frame_b.fy ~ 0
        frame_a.tau + frame_b.tau + r0' * [frame_b.fy, -frame_b.fx] ~ 0
    end
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
@mtkmodel Spring begin
    @extend frame_a, frame_b = partial_frames = PartialTwoFrames()

    @parameters begin
        c_x = 1, [description = "Spring constant in x dir"]
        c_y = 1, [description = "Spring constant in y dir"]
        c_phi = 1.0e5, [description = "Spring constant"]
        s_relx0 = 0, [description = "Unstretched spring length"]
        s_rely0 = 0, [description = "Unstretched spring length"]
        phi_rel0 = 0, [description = "Unstretched spring angle"]
        s_small = 1.e-10,
        [
            description = "Prevent zero-division if distance between frame_a and frame_b is zero"
        ]
        num_windings = 6, [description = "Number of windings of the coil when rendered"]
        color[1:4] = [0, 0.0, 1, 1], [description = "Color of the spring in animations"]
        render = true, [description = "Render the spring in animations"]
        radius = 0.1, [description = "Radius of spring when rendered"]
        N = 200, [description = "Number of points in mesh when rendered"]
    end

    @variables begin
        s_relx(t) = 0
        s_rely(t) = 0
        phi_rel(t) = 0
        f_x(t)
        f_y(t)
    end

    begin
        r_rel_0 = [s_relx, s_rely, 0]
        l = sqrt(r_rel_0' * r_rel_0)
        e_rel_0 = r_rel_0 / max(l, s_small)
    end

    @equations begin
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
    end
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
@mtkmodel Damper begin
    @extend frame_a, frame_b = partial_frames = PartialTwoFrames()

    @parameters begin
        d = 1, [description = "damping constant"]
        s_small = 1.e-10,
        [
            description = "Prevent zero-division if distance between frame_a and frame_b is zero"
        ]
    end

    @variables begin
        r0x(t) = 0
        r0y(t) = 0
        d0x(t) = 0
        d0y(t) = 0
        vx(t) = 0
        vy(t) = 0
        v(t)
        f(t)
    end

    begin
        r0 = [r0x, r0y]
        l = sqrt(r0' * r0)
    end

    @equations begin
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
    end
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
@mtkmodel SpringDamper begin
    @extend frame_a, frame_b = partial_frames = PartialTwoFrames()

    @parameters begin
        c_x = 1, [description = "Spring constant in x dir"]
        c_y = 1, [description = "Spring constant in y dir"]
        c_phi = 1.0e5, [description = "Spring constant in phi dir"]
        d_x = 1, [description = "Damping constant in x dir"]
        d_y = 1, [description = "Damping constant in y dir"]
        d_phi = 1, [description = "Damping constant in phi dir"]
        s_relx0 = 0, [description = "Unstretched spring length"]
        s_rely0 = 0, [description = "Unstretched spring length"]
        phi_rel0 = 0, [description = "Unstretched spring angle"]
        s_small = 1.e-10,
        [
            description = "Prevent zero-division if distance between frame_a and frame_b is zero"
        ]
        num_windings = 6, [description = "Number of windings of the coil when rendered"]
        color[1:4] = [0, 0.0, 1, 1], [description = "Color of the spring in animations"]
        render = true, [description = "Render the spring in animations"]
        radius = 0.1, [description = "Radius of spring when rendered"]
        N = 200, [description = "Number of points in mesh when rendered"]
    end

    @variables begin
        v_relx(t)
        v_rely(t)
        ω_rel(t) = 0
        s_relx(t)
        s_rely(t)
        phi_rel(t) = 0
        f_x(t)
        f_y(t)
        tau(t)
    end

    begin
        r_rel_0 = [s_relx, s_rely, 0]
        l = sqrt(r_rel_0' * r_rel_0)
        e_rel_0 = r_rel_0 / max(l, s_small)
    end

    @equations begin
        s_relx ~ frame_b.x - frame_a.x
        s_rely ~ frame_b.y - frame_a.y
        phi_rel ~ frame_b.phi - frame_a.phi
        v_relx ~ D(s_relx)
        v_rely ~ D(s_rely)
        ω_rel ~ D(phi_rel)

        tau ~ c_phi * (phi_rel - phi_rel0) + d_phi * ω_rel
        frame_a.tau ~ -tau
        frame_b.tau ~ tau
        f_x ~ c_x * (s_relx - s_relx0) + d_x * v_relx
        f_y ~ c_y * (s_rely - s_rely0) + d_y * v_rely
        frame_a.fx ~ -f_x
        frame_b.fx ~ f_x
        frame_a.fy ~ -f_y
        frame_b.fy ~ f_y

        # lossPower ~ d_x * v_relx * v_relx + d_y * v_rely * v_rely
    end
end
