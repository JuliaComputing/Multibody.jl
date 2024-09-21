"""
    RollingWheelJoint(; name, radius, angles, x0, y0, z0)

Joint (no mass, no inertia) that describes an ideal rolling wheel (rolling on the plane y=0). See [`RollingWheel`](@ref) for a realistic wheel model with inertia.

A joint for a wheel rolling on the x-z plane of the world frame.
The rolling contact is considered being ideal, i.e. there is no
slip between the wheel and the ground. This is simply
gained by two non-holonomic constraint equations on velocity level
defined for both longitudinal and lateral direction of the wheel.
There is also a holonomic constraint equation on position level
granting a permanent contact of the wheel to the ground, i.e.
the wheel can not take off.

The origin of the frame `frame_a` is placed in the intersection
of the wheel spin axis with the wheel middle plane and rotates
with the wheel itself. The z-axis of `frame_a` is identical with
the wheel spin axis, i.e. the wheel rotates about z-axis of `frame_a`.
A wheel body collecting the mass and inertia should be connected to
this frame.

# Arguments and parameters:

- `radius`: Radius of the wheel
- `angles`: Angles to rotate world-frame into frame_a around y-, z-, x-axis
- `surface`: By default, the wheel is rolling on a flat xz plane. A function `surface = (x, z)->y` may be provided to define a road surface. The function should return the height of the road at `(x, z)`.

# Variables:
- `x`: x-position of the wheel axis
- `y`: y-position of the wheel axis
- `z`: z-position of the wheel axis
- `angles`: Angles to rotate world-frame into `frame_a` around y-, z-, x-axis
- `der_angles`: Derivatives of angles
- `r_road_0`: Position vector from world frame to contact point on road, resolved in world frame
- `f_wheel_0`: Force vector on wheel, resolved in world frame
- `f_n`: Contact force acting on wheel in normal direction
- `f_lat`: Contact force acting on wheel in lateral direction
- `f_long`: Contact force acting on wheel in longitudinal direction
- `err`: Absolute value of `(r_road_0 - frame_a.r_0) - radius` (must be zero; used for checking)
- `e_axis_0`: Unit vector along wheel axis, resolved in world frame
- `delta_0`: Distance vector from wheel center to contact point
- `e_n_0`: Unit vector in normal direction of road at contact point, resolved in world frame
- `e_lat_0`: Unit vector in lateral direction of road at contact point, resolved in world frame
- `e_long_0`: Unit vector in longitudinal direction of road at contact point, resolved in world frame
- `s`: Road surface parameter 1
- `w`: Road surface parameter 2
- `e_s_0`: Road heading at `(s,w)`, resolved in world frame (unit vector)
- `v_0`: Velocity of wheel center, resolved in world frame
- `w_0`: Angular velocity of wheel, resolved in world frame
- `vContact_0`: Velocity of contact point, resolved in world frame

# Connector frames
- `frame_a`: Frame for the wheel joint
"""
@component function RollingWheelJoint(; name, radius, angles = zeros(3), der_angles=zeros(3), x0=0, y0 = radius, z0=0, sequence = [2, 3, 1], iscut=false, surface = nothing)
    @parameters begin radius = radius, [description = "Radius of the wheel"] end
    @variables begin
        (x(t) = x0), [state_priority = 15, description = "x-position of the wheel axis"]
        (y(t) = y0), [state_priority = 0, description = "y-position of the wheel axis"]
        (z(t) = z0), [state_priority = 15, description = "z-position of the wheel axis"]
        (angles(t)[1:3] = angles),
        [state_priority = 5, description = "Angles to rotate world-frame into frame_a around z-, y-, x-axis"]
        (der_angles(t)[1:3] = der_angles), [state_priority = 5, description = "Derivatives of angles"]
        (r_road_0(t)[1:3] = zeros(3)),
        [
            description = "Position vector from world frame to contact point on road, resolved in world frame",
        ]
        (f_wheel_0(t)[1:3] = zeros(3)),
        [description = "Force vector on wheel, resolved in world frame"]
        (f_n(t) = 0), [description = "Contact force acting on wheel in normal direction"]
        (f_lat(t) = 0), [
            description = "Contact force acting on wheel in lateral direction",
        ]
        (f_long(t) = 0),
        [description = "Contact force acting on wheel in longitudinal direction"]
        # (err(t) = 0),
        # [
        #     description = "|r_road_0 - frame_a.r_0| - radius (must be zero; used for checking)",
        # ]
        (e_axis_0(t)[1:3] = zeros(3)),
        [description = "Unit vector along wheel axis, resolved in world frame"]
        (delta_0(t)[1:3] = [0,-radius, 0]),
        [description = "Distance vector from wheel center to contact point"]
        (e_n_0(t)[1:3] = zeros(3)),
        [
            description = "Unit vector in normal direction of road at contact point, resolved in world frame",
        ]
        (e_lat_0(t)[1:3] = zeros(3)),
        [
            description = "Unit vector in lateral direction of road at contact point, resolved in world frame",
        ]
        (e_long_0(t)[1:3] = zeros(3)),
        [
            description = "Unit vector in longitudinal direction of road at contact point, resolved in world frame",
        ]

        (s(t) = 0), [description = "Road surface parameter 1"]
        (w(t) = 0), [description = "Road surface parameter 2"]
        (e_s_0(t)[1:3] = zeros(3)),
        [description = "Road heading at (s,w), resolved in world frame (unit vector)"]

        (v_0(t)[1:3] = zeros(3)),
        [description = "Velocity of wheel center, resolved in world frame"]
        (w_0(t)[1:3] = zeros(3)),
        [description = "Angular velocity of wheel, resolved in world frame"]
        (vContact_0(t)[1:3] = zeros(3)),
        [description = "Velocity of contact point, resolved in world frame"]

        (aux(t)[1:3] = zeros(3)), [description = "Auxiliary variable"]
    end

    angles,der_angles,r_road_0,f_wheel_0,e_axis_0,delta_0,e_n_0,e_lat_0,e_long_0,e_s_0,v_0,w_0,vContact_0,aux = collect.((angles,der_angles,r_road_0,f_wheel_0,e_axis_0,delta_0,e_n_0,e_lat_0,e_long_0,e_s_0,v_0,w_0,vContact_0,aux))

    @named frame_a = Frame(varw=true)
    Ra = ori(frame_a, true)

    Rarot = axes_rotations(sequence, angles, -der_angles) # The - is the neg_w change

    equations = if surface === nothing
        [ # Road description
            r_road_0 .~ [s, 0, w]
            e_n_0 .~ [0, 1, 0]
            e_s_0 .~ [1, 0, 0]
        ]
     else
        sy = surface(s, w)
        e_w_0 = _normalize([0, expand_derivatives(Differential(w)(sy)), 1])
        # @show sy, expand_derivatives(Differential(s)(sy)), expand_derivatives(Differential(w)(sy))
        [
            r_road_0 .~ [s, sy, w]
            e_s_0 .~ _normalize([1, expand_derivatives(Differential(s)(sy)), 0])
            e_n_0 .~ _normalize(cross(e_w_0, e_s_0))
        ]
     end

    equations = [
                equations;
                connect_orientation(Ra, Rarot; iscut)   # Ra ~ Rarot
                Ra.w ~ Rarot.w

                 # frame_a.R is computed from generalized coordinates
                 collect(frame_a.r_0) .~ [x, y, z]
                 der_angles .~ D.(angles)


                 # Coordinate system at contact point (e_long_0, e_lat_0, e_n_0)
                 e_axis_0 .~ resolve1(Ra, [0, 0, 1])
                 aux .~ (cross(e_n_0, e_axis_0))
                 e_long_0 .~ (aux ./ _norm(aux))
                 e_lat_0 .~ (cross(e_long_0, e_n_0))

                 # Determine point on road where the wheel is in contact with the road
                 delta_0 .~ r_road_0 - frame_a.r_0
                 0 ~ delta_0'e_axis_0
                 0 ~ delta_0'e_long_0

                 # One holonomic positional constraint equation (no penetration in to the ground)
                 0 ~ radius - delta_0'cross(e_long_0, e_axis_0)

                 # only for testing
                #  err ~ norm(delta_0) - radius

                 # Slip velocities
                 v_0 .~ D.(frame_a.r_0)
                 w_0 .~ angular_velocity1(Ra)
                 vContact_0 .~ v_0 + cross(w_0, delta_0)

                 # Two non-holonomic constraint equations on velocity level (ideal rolling, no slippage)
                 0 ~ vContact_0'e_long_0
                 0 ~ vContact_0'e_lat_0

                 # Contact force
                 f_wheel_0 .~ f_n * e_n_0 + f_lat * e_lat_0 + f_long * e_long_0

                 # Force and torque balance at the wheel center
                 zeros(3) .~ collect(frame_a.f) + resolve2(Ra, f_wheel_0)
                 zeros(3) .~ collect(frame_a.tau) +
                             resolve2(Ra, cross(delta_0, f_wheel_0))]
    compose(ODESystem(equations, t; name), frame_a)
end

"""
    RollingWheel(; name, radius, m, I_axis, I_long, width=0.035, x0, y0, kwargs...)

Ideal rolling wheel on flat surface y=0 (5 positional, 3 velocity degrees of freedom)

A wheel rolling on the x-z plane of the world frame including wheel mass.
The rolling contact is considered being ideal, i.e. there is no
slip between the wheel and the ground.
The wheel can not take off but it can incline toward the ground.
The frame `frame_a` is placed in the wheel center point and rotates
with the wheel itself. A [`Revolute`](@ref) joint rotationg around `n = [0, 1, 0]` is required to attach the wheel to a wheel axis.

# Arguments and parameters:
- `name`: Name of the rolling wheel component
- `radius`: Radius of the wheel
- `m`: Mass of the wheel
- `I_axis`: Moment of inertia of the wheel along its axis
- `I_long`: Moment of inertia of the wheel perpendicular to its axis
- `width`: Width of the wheel (default: 0.035)
- `x0`: Initial x-position of the wheel axis
- `z0`: Initial z-position of the wheel axis
- `kwargs...`: Additional keyword arguments passed to the `RollingWheelJoint` function

# Variables:
- `x`: x-position of the wheel axis
- `z`: z-position of the wheel axis
- `angles`: Angles to rotate world-frame into `frame_a` around y-, z-, x-axis
- `der_angles`: Derivatives of angles  (y: like rotational velocity of a spinning coin, z: wheel forward spin speed, x: wheel falling over speed)

# Named components:
- `frame_a`: Frame for the wheel component
- `wheeljoint`: Rolling wheel joint representing the wheel's contact with the road surface

!!! tip "Rendering tip"
    Due to the symmetry of the wheel, it can be hard to discern how the wheel is rotating in animations. Try enabling rendering of the frame of the wheel by setting 
    ```
    wheel.frame_a.render => true;
    wheel.frame_a.length => 1.1radius;
    wheel.frame_a.radius => 0.02radius;
    ```
"""
@component function RollingWheel(; name, radius, m, I_axis, I_long, width = 0.035, x0=0, z0=0,
                      angles = zeros(3), der_angles = zeros(3), kwargs...)

    @named wheeljoint = RollingWheelJoint(; radius, angles, x0, z0, der_angles, kwargs...)
    @named begin
        frame_a = Frame()
        body = Body(r_cm = [0, 0, 0],
                    state_priority = 0,
                    m = m,
                    I_11 = I_long,
                    I_22 = I_long,
                    I_33 = I_axis,
                    I_21 = 0,
                    I_31 = 0,
                    I_32 = 0)
    end
    pars = @parameters begin
        radius = radius, [description = "Radius of the wheel"]
        m = m, [description = "Mass of the wheel"]
        I_axis = I_axis, [description = "Moment of inertia of the wheel along its axis"]
        I_long = I_long,
                 [description = "Moment of inertia of the wheel perpendicular to its axis"]
        width = width, [description = "Width of the wheel"]
    end
    sts = @variables begin
        (x(t) = x0), [state_priority = 20, description = "x-position of the wheel axis"]
        (z(t) = z0), [state_priority = 20, description = "z-position of the wheel axis"]
        (angles(t)[1:3] = angles),
        [state_priority = 30, description = "Angles to rotate world-frame into frame_a around y-, z-, x-axis"]
        (der_angles(t)[1:3] = der_angles), [state_priority = 30, description = "Derivatives of angles"]
    end
    # sts = reduce(vcat, collect.(sts))

    equations = Equation[wheeljoint.x ~ x
                         wheeljoint.z ~ z
                         collect(wheeljoint.angles) .~ collect(angles)
                         collect(wheeljoint.der_angles) .~ collect(der_angles)
                         connect(body.frame_a, frame_a)
                         connect(wheeljoint.frame_a, frame_a)]
    compose(ODESystem(equations, t; name), frame_a, wheeljoint, body)
end



"""
    SlipWheelJoint(; name, radius, angles = zeros(3), der_angles = zeros(3), x0 = 0, y0 = radius, z0 = 0, sequence, iscut = false, surface = nothing, vAdhesion_min = 0.1, vSlide_min = 0.1, sAdhesion = 0.1, sSlide = 0.1, mu_A = 0.8, mu_S = 0.6, phi_roll = 0, w_roll = 0)

Joint for a wheel with slip rolling on a surface.

# Parameters
- `radius`: Radius of the wheel
- `vAdhesion_min`: Minimum adhesion velocity
- `vSlide_min`: Minimum sliding velocity
- `sAdhesion`: Adhesion slippage
- `sSlide`: Sliding slippage
- `mu_A`: Friction coefficient at adhesion
- `mu_S`: Friction coefficient at sliding
"""
@component function SlipWheelJoint(; name, radius, angles = zeros(3), der_angles=zeros(3), x0=0, y0 = radius, z0=0, sequence = [2, 3, 1], iscut=false, surface = nothing, vAdhesion_min = 0.1, vSlide_min = 0.1, sAdhesion = 0.1, sSlide = 0.1, mu_A = 0.8, mu_S = 0.6, phi_roll = 0, w_roll = 0)
    @parameters begin
        radius = radius, [description = "Radius of the wheel"]
        vAdhesion_min = vAdhesion_min, [description = "Minimum adhesion velocity"]
        vSlide_min = vSlide_min, [description = "Minimum sliding velocity"]
        sAdhesion = sAdhesion, [description = "Adhesion slippage"]
        sSlide = sSlide, [description = "Sliding slippage"]
        mu_A = mu_A, [description = "Friction coefficient at adhesion"]
        mu_S = mu_S, [description = "Friction coefficient at sliding"]
    end
    @variables begin
        (x(t) = x0), [state_priority = 15, description = "x-position of the wheel axis"]
        (y(t) = y0), [state_priority = 0, description = "y-position of the wheel axis"]
        (z(t) = z0), [state_priority = 15, description = "z-position of the wheel axis"]
        (angles(t)[1:3] = angles),
        [state_priority = 5, description = "Angles to rotate world-frame into frame_a around z-, y-, x-axis"]
        (der_angles(t)[1:3] = der_angles), [state_priority = 5, description = "Derivatives of angles"]
        (phi_roll(t) = phi_roll), [guess=0, description="wheel angle"] # wheel angle
        (w_roll(t)=w_roll), [guess=0, description="Roll velocity of wheel"]
        (r_road_0(t)[1:3] = zeros(3)),
        [
            description = "Position vector from world frame to contact point on road, resolved in world frame",
        ]
        (f_wheel_0(t)[1:3] = zeros(3)),
        [description = "Force vector on wheel, resolved in world frame"]
        (f_n(t) = 0), [description = "Contact force acting on wheel in normal direction"]
        (f_lat(t) = 0), [
            description = "Contact force acting on wheel in lateral direction",
        ]
        (f_long(t) = 0),
        [description = "Contact force acting on wheel in longitudinal direction"]
        # (err(t) = 0),
        # [
        #     description = "|r_road_0 - frame_a.r_0| - radius (must be zero; used for checking)",
        # ]
        (e_axis_0(t)[1:3] = zeros(3)),
        [description = "Unit vector along wheel axis, resolved in world frame"]
        (delta_0(t)[1:3] = [0,-radius, 0]),
        [description = "Distance vector from wheel center to contact point"]
        (e_n_0(t)[1:3] = zeros(3)),
        [
            description = "Unit vector in normal direction of road at contact point, resolved in world frame",
        ]
        (e_lat_0(t)[1:3] = zeros(3)),
        [
            description = "Unit vector in lateral direction of road at contact point, resolved in world frame",
        ]
        (e_long_0(t)[1:3] = zeros(3)),
        [
            description = "Unit vector in longitudinal direction of road at contact point, resolved in world frame",
        ]

        (s(t) = 0), [description = "Road surface parameter 1"]
        (w(t) = 0), [description = "Road surface parameter 2"]
        (e_s_0(t)[1:3] = zeros(3)),
        [description = "Road heading at (s,w), resolved in world frame (unit vector)"]

        (v_0(t)[1:3] = zeros(3)),
        [description = "Velocity of wheel center, resolved in world frame"]
        (w_0(t)[1:3] = zeros(3)),
        [description = "Angular velocity of wheel, resolved in world frame"]
        (vContact_0(t)[1:3] = zeros(3)),
        [description = "Velocity of contact point, resolved in world frame"]

        (aux(t)[1:3] = zeros(3)), [description = "Auxiliary variable"]

        # New variables ========================================================
        v_lat(t), [guess=0, description="Velocity in lateral direction"]
        v_long(t), [guess=0, description="Velocity in longitudinal direction"]
        v_slip_long(t), [guess=0, description="Slip velocity in longitudinal direction"]
        v_slip_lat(t), [guess=0, description="Slip velocity in lateral direction"]
        v_slip(t), [description="Slip velocity, norm of component slip velocities"]
        f(t), [description="Total traction force"]
        vAdhesion(t), [description="Adhesion velocity"]
        vSlide(t), [description="Sliding velocity"]
    end

    angles,der_angles,r_road_0,f_wheel_0,e_axis_0,delta_0,e_n_0,e_lat_0,e_long_0,e_s_0,v_0,w_0,vContact_0,aux = collect.((angles,der_angles,r_road_0,f_wheel_0,e_axis_0,delta_0,e_n_0,e_lat_0,e_long_0,e_s_0,v_0,w_0,vContact_0,aux))

    @named frame_a = Frame(varw=true)
    Ra = ori(frame_a, true)

    Rarot = axes_rotations(sequence, angles, -der_angles) # The - is the neg_w change

    equations = if surface === nothing
        [ # Road description
            r_road_0 .~ [s, 0, w]
            e_n_0 .~ [0, 1, 0]
            e_s_0 .~ [1, 0, 0]
        ]
    else
        sy = surface(s, w)
        e_w_0 = _normalize([0, expand_derivatives(Differential(w)(sy)), 1])
        # @show sy, expand_derivatives(Differential(s)(sy)), expand_derivatives(Differential(w)(sy))
        [
            r_road_0 .~ [s, sy, w]
            e_s_0 .~ _normalize([1, expand_derivatives(Differential(s)(sy)), 0])
            e_n_0 .~ _normalize(cross(e_w_0, e_s_0))
        ]
    end


    equations = [
                equations;
                connect_orientation(Ra, Rarot; iscut)   # Ra ~ Rarot
                Ra.w ~ Rarot.w

                phi_roll ~ angles[2]
                w_roll ~ D(phi_roll)

                # frame_a.R is computed from generalized coordinates
                collect(frame_a.r_0) .~ [x, y, z]
                der_angles .~ D.(angles)


                # Coordinate system at contact point (e_long_0, e_lat_0, e_n_0)
                e_axis_0 .~ resolve1(Ra, [0, 0, 1])
                aux .~ (cross(e_n_0, e_axis_0))
                e_long_0 .~ (aux ./ _norm(aux))
                e_lat_0 .~ (cross(e_long_0, e_n_0))

                # Determine point on road where the wheel is in contact with the road
                delta_0 .~ r_road_0 - frame_a.r_0
                0 ~ delta_0'e_axis_0
                0 ~ delta_0'e_long_0

                # One holonomic positional constraint equation (no penetration in to the ground)
                0 ~ radius - delta_0'cross(e_long_0, e_axis_0)

                # Slip velocities
                v_0 .~ D.(frame_a.r_0)
                w_0 .~ angular_velocity1(Ra)
                vContact_0 .~ v_0 + cross(w_0, delta_0)

                # Contact dynamics =============================================

                v_slip_lat ~ vContact_0' * e_lat_0
                v_slip_long ~ vContact_0' * e_long_0
                # v_slip_lat ~ v_lat - 0
                # v_slip_long ~ v_long - radius * w_roll

                v_slip ~ sqrt(v_slip_long^2 + v_slip_lat^2) + 0.00001
                # -f_long * radius ~ flange_a.tau # No longer needed?
                # frame_a.tau ~ 0
                vAdhesion ~ max(vAdhesion_min, sAdhesion * abs(radius * w_roll))
                vSlide ~ max(vSlide_min, sSlide * abs(radius * w_roll))

                f ~ f_n * PlanarMechanics.limit_S_triple(vAdhesion, vSlide, mu_A, mu_S, v_slip) # limit_S_triple(x_max, x_sat, y_max, y_sat, x)
                f_long ~ -f * v_slip_long / v_slip
                f_lat ~ -f * v_slip_lat / v_slip

                # Contact force
                f_wheel_0 .~ f_n * e_n_0 + f_lat * e_lat_0 + f_long * e_long_0

                # Force and torque balance at the wheel center
                zeros(3) .~ collect(frame_a.f) + resolve2(Ra, f_wheel_0)
                zeros(3) .~ collect(frame_a.tau) +
                            resolve2(Ra, cross(delta_0, f_wheel_0))]
    compose(ODESystem(equations, t; name), frame_a)
end


"""
    SlippingWheel(; name, radius, m, I_axis, I_long, width = 0.035, x0=0, z0=0,
                      angles = zeros(3), der_angles = zeros(3), kwargs...)

Wheel with slip rolling on a surface.

# Parameters
- `radius`: Radius of the wheel
- `m`: Mass of the wheel
- `I_axis`: Moment of inertia of the wheel along its axis
- `I_long`: Moment of inertia of the wheel perpendicular to its axis
- `width`: Width of the wheel (for rendering)
- `x0`: Initial x-position of the wheel axis
- `z0`: Initial z-position of the wheel axis

# Variables
- `x`: x-position of the wheel axis
- `z`: z-position of the wheel axis
- `angles`: Angles to rotate world-frame into `frame_a` around y-, z-, x-axis
- `der_angles`: Derivatives of angles

# Connectors
- `frame_a`: Frame for the wheel component

# Examples
See [Docs: Wheels](https://help.juliahub.com/multibody/dev/examples/wheel/)
"""
@component function SlippingWheel(; name, radius, m, I_axis, I_long, width = 0.035, x0=0, z0=0,
                      angles = zeros(3), der_angles = zeros(3), kwargs...)
    @named wheeljoint = SlipWheelJoint(; radius, angles, x0, z0, der_angles, kwargs...)
    @named begin
        frame_a = Frame()
        body = Body(r_cm = [0, 0, 0],
                    state_priority = 0,
                    m = m,
                    I_11 = I_long,
                    I_22 = I_long,
                    I_33 = I_axis,
                    I_21 = 0,
                    I_31 = 0,
                    I_32 = 0)
    end
    pars = @parameters begin
        radius = radius, [description = "Radius of the wheel"]
        m = m, [description = "Mass of the wheel"]
        I_axis = I_axis, [description = "Moment of inertia of the wheel along its axis"]
        I_long = I_long,
                 [description = "Moment of inertia of the wheel perpendicular to its axis"]
        width = width, [description = "Width of the wheel"]
    end
    sts = @variables begin
        (x(t) = x0), [state_priority = 20, description = "x-position of the wheel axis"]
        (z(t) = z0), [state_priority = 20, description = "z-position of the wheel axis"]
        (angles(t)[1:3] = angles),
        [state_priority = 30, description = "Angles to rotate world-frame into frame_a around y-, z-, x-axis"]
        (der_angles(t)[1:3] = der_angles), [state_priority = 30, description = "Derivatives of angles"]
    end
    # sts = reduce(vcat, collect.(sts))

    equations = Equation[wheeljoint.x ~ x
                         wheeljoint.z ~ z
                         collect(wheeljoint.angles) .~ collect(angles)
                         collect(wheeljoint.der_angles) .~ collect(der_angles)
                         connect(body.frame_a, frame_a)
                         connect(wheeljoint.frame_a, frame_a)]
    compose(ODESystem(equations, t; name), frame_a, wheeljoint, body)
end

"""
    RollingConstraintVerticalWheel(;
        name,
        radius = 0.3,
        lateral_sliding_constraint = true,
    )

Rolling constraint for wheel that is always perpendicular to x-z plane

Joint for a wheel rolling on the x-z plane of the world frame intended for an idealized wheelset. To meet this objective, the wheel always runs upright and enables no slip in the longitudinal direction of the wheel/ground contact. On the contrary, the wheel can optionally slip in the lateral direction which is reasonable for the wheelset where just one of the wheels should be laterally constrained. The frame `frame_a` is placed in the intersection of the wheel spin axis with the wheel middle plane and rotates with the wheel itself. A wheel body collecting the mass and inertia

# Arguments and parameters:
- `name`: Name of the rolling wheel joint component
- `radius`: Wheel radius
- `lateral_sliding_constraint`: true, if lateral sliding constraint taken into account, = false if lateral force = 0 (needed to avoid overconstraining if two ideal rolling wheels are connect on one axis)

# Connectors:
- `frame_a`: Frame for the wheel joint
"""
function RollingConstraintVerticalWheel(;
    name,
    radius = 0.3,
    lateral_sliding_constraint = true,
)
    @named frame_a = Frame()
    @parameters begin
        radius = radius, [description = "Wheel radius"]
        # lateral_sliding_constraint = true,
        # [description = "true, if lateral sliding constraint taken into account, = false if lateral force = 0 (needed to avoid overconstraining if two ideal rolling wheels are connect on one axis)"]
    end
    @variables begin
        (f_wheel_0(t)[1:3] = zeros(3)),
        [description = "Contact force acting on wheel, resolved in world frame"]
        (f_lat(t) = 0), [description = "Contact force acting on wheel in lateral direction"]
        (f_long(t) = 0),
        [description = "Contact force acting on wheel in longitudinal direction"]
        (e_axis_0(t)[1:3] = zeros(3)),
        [description = "Unit vector along wheel axis, resolved in world frame"]
        # (rContact_0(t)[1:3] = [0, -radius, 0]),
        # [description = "Distance vector from wheel center to contact point, resolved in world frame"]
        # (e_n_0(t)[1:3] = [0, 1, 0]),
        # [description = "Unit vector in normal direction of road at contact point, resolved in world frame"]
        (e_lat_0(t)[1:3] = zeros(3)),
        [description = "Unit vector in lateral direction of wheel at contact point, resolved in world frame"]
        (e_long_0(t)[1:3] = zeros(3)),
        [description = "Unit vector in longitudinal direction of wheel at contact point, resolved in world frame"]
        (v_0(t)[1:3] = zeros(3)),
        [description = "Velocity of wheel center, resolved in world frame"]
        (w_0(t)[1:3] = zeros(3)),
        [description = "Angular velocity of wheel, resolved in world frame"]
        (vContact_0(t)[1:3] = zeros(3)),
        [description = "Velocity of wheel contact point, resolved in world frame"]
        (aux(t)[1:3] = zeros(3)),
        [description = "Auxiliary variable"]
    end
    e_n_0 = [0, 1, 0]
    rContact_0 = [0, -radius, 0]

    # e_n_0 = [0, 0, 1]
    # rContact_0 = [0, 0, -radius]

    f_wheel_0, e_lat_0, e_long_0, vContact_0, aux, v_0 = collect.((f_wheel_0, e_lat_0, e_long_0, vContact_0, aux, v_0))

    equations = Equation[
        # Coordinate system at contact point (e_long_0, e_lat_0, e_n_0)
        e_axis_0 .~ resolve1(ori(frame_a), [0, 0, 1])
        aux .~ cross(e_n_0, e_axis_0)
        e_long_0 .~ aux ./ _norm(aux)
        e_lat_0 .~ cross(e_long_0, e_n_0)

        # Slip velocities
        v_0 .~ D.(frame_a.r_0)
        w_0 .~ angular_velocity1(ori(frame_a))
        vContact_0 .~ v_0 + cross(w_0, rContact_0)

        # Two non-holonomic constraint equations on velocity level (ideal rolling, no slippage)
        0 ~ vContact_0'e_long_0
        if lateral_sliding_constraint
            [
                0 ~ vContact_0'e_lat_0
                f_wheel_0 .~ collect(f_lat * e_lat_0) + collect(f_long * e_long_0)
            ]
        else
            [
                0 ~ f_lat
                f_wheel_0 .~ f_long .* e_long_0
            ]
        end
        zeros(3) .~ collect(frame_a.f) + resolve2(ori(frame_a), f_wheel_0)
        zeros(3) .~ collect(frame_a.tau) + resolve2(ori(frame_a), cross(rContact_0, f_wheel_0))
    ]
    ODESystem(equations, t; name, systems = [frame_a])
end


"""
    RollingWheelSetJoint(;
        name,
        radius = 0.3,
        track = 1.0,
        state_priority = 1,
        x0 = 0,
        z0 = 0,
        phi0 = 0,
        theta1_0 = 0,
        theta2_0 = 0,
        der_theta1_0 = 0,
        der_theta2_0 = 0,
        render = true,
        iscut = false,
    )

Joint (no mass, no inertia) that describes an ideal rolling wheel set (two ideal rolling wheels connected together by an axis)

An assembly joint for a wheelset rolling on the x-z plane of the world frame. The frames `frame1` and `frame2` are connected to rotating wheels; the `frame_middle` moves in a plane parallel to the x-z plane of the world and should be connected to the vehicle body.

To work properly, the gravity acceleration vector g of the world must point in the negative y-axis (default)

# Arguments and parameters:
- `iscut`: if more than one wheel set is connected to the same rigid body, `iscut` must be set to true for all but one set. This avoids overconstraining the system by replacing the planar joint giving the set coordinates by an unconstrained FreeMotion joint.
- `radius`: Radius of one wheel
- `track`: Distance between the two wheels (= axle track)

# Connectors:
- `frame_middle`: Frame fixed in middle of axis connecting both wheels (z-axis: along wheel axis, y-axis: upwards)
- `frame1`: Frame fixed in center point of left wheel, rotating with the wheel (z-axis: along wheel axis, y-axis: upwards when wheel angle is zero)
- `frame2`: Frame fixed in center point of right wheel, rotating with the wheel (z-axis: along wheel axis, y-axis: upwards when wheel angle is zero)
- `axis1`: 1-dim. Rotational flange that drives the joint
- `axis2`: 1-dim. Rotational flange that drives the joint
- `support`: Support of 1-dim axes
"""
function RollingWheelSetJoint(;
    name,
    radius = 0.3,
    track = 1.0,
    state_priority = 100,
    x0 = 0,
    z0 = 0,
    phi0 = 0,
    theta1_0 = 0,
    theta2_0 = 0,
    der_theta1_0 = 0,
    der_theta2_0 = 0,
    render = true,
    color = [0, 1, 1, 0.1],
    width_wheel = 0.15*radius,
    iscut = false,
)
    pars = @parameters begin
        # radius = radius, [description = "Radius of one wheel"]
        # track = track, [description = "Distance between the two wheels (= axle track)"]
        width_wheel = width_wheel, [description = "Width of one wheel"]
        # color[1:4] = color, [description = "Color of the wheel set in animations"]
    end
    
    systems = @named begin
        frame_middle = Frame()
        frame1 = Frame()
        frame2 = Frame()
        fixed = Fixed(; r = [0, radius, 0], render)
        rod1 = FixedTranslation(; r = [0, 0, track / 2], render, color)
        rod2 = FixedTranslation(; r = [0, 0, -track / 2], render, color)
        revolute1 = Revolute(; n = [0, 0, 1], axisflange = true, render, state_priority=11, radius, length=width_wheel, color)
        revolute2 = Revolute(; n = [0, 0, 1], axisflange = true, render, state_priority=11, radius, length=width_wheel, color)
        rolling1 = RollingConstraintVerticalWheel(; radius, lateral_sliding_constraint = true)
        rolling2 = RollingConstraintVerticalWheel(; radius, lateral_sliding_constraint = false)
        axis1 = Rotational.Flange()
        axis2 = Rotational.Flange()
        mounting1D = Mounting1D()
        support = Rotational.Flange()
    end
    if iscut
        @named freemotion = FreeMotion(iscut=false, isroot=true, state=false)
        push!(systems, freemotion)
    else
        more = @named begin
            revolute = Revolute(; render, n = [0, 1, 0], color)
            prismatic1 = Prismatic(; n = [1, 0, 0], render, state_priority=10, color=[0,1,1,0.1], radius=0.01)
            prismatic2 = Prismatic(; n = [0, 0, 1], render, state_priority=10, color=[0,1,1,0.1], radius=0.01)
        end
        append!(systems, more)
    end

    sts = @variables begin
        (x(t) = x0), [description = "x coordinate for center between wheels", state_priority = state_priority]
        (z(t) = z0), [description = "z coordinate for center between wheels", state_priority = state_priority]
        (phi(t) = phi0), [description = "Orientation angle of wheel axis along y-axis", state_priority = state_priority]
        (theta1(t) = theta1_0), [description = "Angle of wheel 1", state_priority = state_priority]
        (theta2(t) = theta2_0), [description = "Angle of wheel 2", state_priority = state_priority]
        (der_theta1(t) = der_theta1_0), [description = "Derivative of theta 1", state_priority = state_priority]
        (der_theta2(t) = der_theta2_0), [description = "Derivative of theta 2", state_priority = state_priority]
    end
    equations = Equation[
        if iscut
            [
                connect(freemotion.frame_a, fixed.frame_b)
                connect(freemotion.frame_b, frame_middle)
                freemotion.r_rel_a[1] ~ x
                freemotion.r_rel_a[3] ~ z
                # freemotion.phi[2] ~ phi
                phi ~ 0
            ]
        else    
            [prismatic1.s ~ x
            prismatic2.s ~ z
            revolute.phi ~ phi
            connect(revolute.frame_b, frame_middle)
            connect(prismatic1.frame_a, fixed.frame_b)
            connect(prismatic1.frame_b, prismatic2.frame_a)
            connect(prismatic2.frame_b, revolute.frame_a)]
        end

        revolute1.phi ~ theta1
        revolute2.phi ~ theta2
        der_theta1 ~ D(theta1)
        der_theta2 ~ D(theta2)

        connect(rod1.frame_a, frame_middle)
        connect(rod2.frame_a, frame_middle)
        connect(rod1.frame_b, revolute1.frame_a)
        connect(revolute1.frame_b, frame1)
        connect(revolute2.frame_a, rod2.frame_b)
        connect(revolute2.frame_b, frame2)
       
        connect(rolling1.frame_a, revolute1.frame_b)
        connect(rolling2.frame_a, revolute2.frame_b)
        connect(revolute1.axis, axis1)
        connect(revolute2.axis, axis2)
        connect(frame_middle, mounting1D.frame_a)
        connect(mounting1D.flange_b, support)
    ]
    sys = ODESystem(equations, t; name=:nothing, systems)
    add_params(sys, [width_wheel]; name)
end

"""
    RollingWheelSet(;
        name,
        radius = 0.3,
        m_wheel = 1.0,
        I_axis = 1.0,
        I_long = 1.0,
        track = 1.0,
        state_priority = 1,
        x0 = 0,
        z0 = 0,
        phi0 = 0,
        theta1_0 = 0,
        theta2_0 = 0,
        der_theta1_0 = 0,
        der_theta2_0 = 0,
        width_wheel = 0.01,
        color = [0.3, 0.3, 0.3, 1],
        render = true,
        iscut = false,
    )

Ideal rolling wheel set consisting of two ideal rolling wheels connected together by an axis

Two wheels are connected by an axis and can rotate around this axis. The wheels are rolling on the x-z plane of the world frame. The coordinate system attached to the center of the wheel axis (`frame_middle`) is constrained so that it is always parallel to the x-z plane. If all generalized coordinates are zero, `frame_middle` is parallel to the world frame.

# Arguments and parameters:
- `iscut`: if more than one wheel set is connected to the same rigid body, `iscut` must be set to true for all but one set. This avoids overconstraining the system by replacing the planar joint giving the set coordinates by an unconstrained FreeMotion joint.
- `radius`: Radius of one wheel
- `m_wheel`: Mass of one wheel
- `I_axis`: Moment of inertia of one wheel around its rotation axis
- `I_long`: Moment of inertia of one wheel perpendicular to its rotation axis
- `track`: Distance between the two wheels (= axle track)

# Connectors
- `frame_middle`: Frame fixed in middle of axis connecting both wheels (z-axis: along wheel axis, y-axis: upwards)
- `frame1`: Frame fixed in center point of left wheel, rotating with the wheel (z-axis: along wheel axis, y-axis: upwards when wheel angle is zero)
- `frame2`: Frame fixed in center point of right wheel, rotating with the wheel (z-axis: along wheel axis, y-axis: upwards when wheel angle is zero)
- `axis1`: 1-dim. Rotational flange that drives the left wheel
- `axis2`: 1-dim. Rotational flange that drives the right wheel
- `support`: Support of 1D axes

To connect driving torques or friction to rotation of the wheels, connect between `axis1` and `support`, and `axis2` and `support` respectively. To connect the wheel set to, e.g., a body, connect the `frame_middle` to the body frame.
"""
function RollingWheelSet(;
    name,
    radius = 0.3,
    m_wheel = 1.0,
    I_axis = 1.0,
    I_long = 1.0,
    track = 1.0,
    state_priority = 1,
    x0 = 0,
    z0 = 0,
    phi0 = 0,
    theta1_0 = 0,
    theta2_0 = 0,
    der_theta1_0 = 0,
    der_theta2_0 = 0,
    width_wheel = 0.01,
    hollow_fraction = 0.8,
    color = [0.3, 0.3, 0.3, 1],
    render = true,
    kwargs...
)
    pars = @parameters begin
        # radius = radius, [description = "Radius of one wheel"]
        width_wheel = width_wheel, [description = "Width of one wheel"]
        m_wheel = m_wheel, [description = "Mass of one wheel"]
        I_axis = I_axis, [description = "Inertia along one wheel axis"]
        I_long = I_long, [description = "Inertia perpendicular to one wheel axis"]
        # track = track, [description = "Distance between the two wheels (= axle track)"]
        # hollow_fraction = hollow_fraction,
        # [description = "For ring-like wheel visualization: wheel radius / inner hole radius; i.e. 1.0: completely hollow, 0.0: full disc"]
        # color[1:4] = color, [description = "Color of wheels"]
    end
    systems = @named begin
        frame_middle = Frame()
        frame1 = Frame()
        frame2 = Frame()
        body2 = Body(r_cm = [0, 0, 0],
                    m = m_wheel,
                    I_11 = I_long,
                    I_22 = I_long,
                    I_33 = I_axis,
                    sparse_I = true,
                    render = false)
        body1 = Body(r_cm = [0, 0, 0],
                    m = m_wheel,
                    I_11 = I_long,
                    I_22 = I_long,
                    I_33 = I_axis,
                    sparse_I = true,
                    render = false)
        axis1 = Rotational.Flange()
        axis2 = Rotational.Flange()
        wheelSetJoint = RollingWheelSetJoint(; radius, track, state_priority, x0, z0, phi0, theta1_0, theta2_0, der_theta1_0, der_theta2_0, render, width_wheel, color, kwargs...)
        support = Rotational.Flange()
    end

    sts = @variables begin
        (x(t) = x0), [description = "x coordinate of center between wheels", state_priority = state_priority]
        (z(t) = z0), [description = "z coordinate of center between wheels", state_priority = state_priority]
        (phi(t) = phi0), [description = "Orientation angle of wheel axis along z-axis", state_priority = state_priority]
        (theta1(t) = theta1_0), [description = "Angle of wheel 1", state_priority = state_priority]
        (theta2(t) = theta2_0), [description = "Angle of wheel 2", state_priority = state_priority]
        (der_theta1(t) = der_theta1_0), [description = "Derivative of theta 1", state_priority = state_priority]
        (der_theta2(t) = der_theta2_0), [description = "Derivative of theta 2", state_priority = state_priority]
    end

    equations = [
        wheelSetJoint.x ~ x
        wheelSetJoint.z ~ z
        wheelSetJoint.phi ~ phi
        wheelSetJoint.theta1 ~ theta1
        wheelSetJoint.theta2 ~ theta2
        der_theta1 ~ D(theta1)
        der_theta2 ~ D(theta2)

        connect(body2.frame_a, frame2)
        connect(body1.frame_a, frame1)
        connect(wheelSetJoint.frame2, frame2)
        connect(wheelSetJoint.frame1, frame1)
        connect(wheelSetJoint.axis1, axis1)
        connect(wheelSetJoint.axis2, axis2)
        connect(wheelSetJoint.support, support)
        connect(wheelSetJoint.frame_middle, frame_middle)
    ]

    sys = ODESystem(equations, t; name=:nothing, systems)
    add_params(sys, [width_wheel]; name)

end

