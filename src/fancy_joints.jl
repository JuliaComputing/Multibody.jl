import ModelingToolkitStandardLibrary.Blocks

"""
    SphericalSpherical(; name, state = false, isroot = true, iscut=false, w_rel_a_fixed = false, r_0 = [0,0,0], color = [1, 1, 0, 1], m = 0, radius = 0.1, kinematic_constraint=true)

Joint that has a spherical joint on each of its two ends. The rod connecting the two spherical joints is approximated by a point mass that is located in the middle of the rod. When the mass is set to zero (default), special code for a massless body is generated. 

This joint introduces one constraint defining that the distance between the origin of `frame_a` and the origin of `frame_b` is constant (= rodLength). It is highly recommended to use this joint in loops whenever possible, because this enhances the efficiency considerably due to smaller systems of non-linear algebraic equations.

It is not possible to connect other components, such as a body with mass properties or a special visual shape object to the rod connecting the two spherical joints. If this is needed, use instead joint [`UniversalSpherical`](@ref) that has the additional frame `frame_ia` for this.

# Connectors:
- `frame_a`: Frame for the first spherical joint
- `frame_b`: Frame for the second spherical joint

# Rendering parameters:
- `radius`: Radius of the joint in animations
- `color`: Color of the joint in animations (RGBA)
"""
@component function SphericalSpherical(; name, state = false, isroot = true, iscut=false, w_rel_a_fixed = false,
                    r_0 = state ? [0,0,0] : nothing,
                   color = [1, 1, 0, 1],
                   m = 0,
                   radius = 0.1,
                   kinematic_constraint=true
                   )

    @named begin
        ptf = PartialTwoFrames()
    end
    pars = @parameters begin
        radius = radius, [description = "radius of the joint in animations"]
        color[1:4] = color, [description = "color of the joint in animations (RGBA)"]
    end
    @unpack frame_a, frame_b = ptf
    # @parameters begin # Currently not using parameters due to these appearing in if statements
    #     sequence[1:3] = sequence
    # end

    @variables f_rod(t), [description="Constraint force in direction of the rod (positive on frame_a, when directed from frame_a to frame_b)";]
    @variables rRod_0(t)[1:3]=r_0, [description="Position vector from frame_a to frame_b resolved in world frame";]
    @variables rRod_a(t)[1:3], [description="Position vector from frame_a to frame_b resolved in frame_a";]
    @variables eRod_a(t)[1:3], [description="Unit vector in direction from frame_a to frame_b, resolved in frame_a";]
    @variables r_cm_0(t)[1:3], [description="Dummy if m==0, or position vector from world frame to mid-point of rod, resolved in world frame";]
    @variables v_cm_0(t)[1:3], [description="First derivative of r_cm_0";]
    @variables f_cm_a(t)[1:3], [description="Dummy if m==0, or inertial force acting at mid-point of rod due to mass point acceleration, resolved in frame_a";]
    @variables f_cm_e(t)[1:3], [description="Dummy if m==0, or projection of f_cm_a onto eRod_a, resolved in frame_a";]
    @variables f_b_a1(t)[1:3], [description="Force acting at frame_b, but without force in rod, resolved in frame_a";]

    rodlength = _norm(r_0)
    constraint_residue = rRod_0'rRod_0 - rodlength^2


    eqs = [
        # Determine relative position vector between the two frames
  if kinematic_constraint
    rRod_0 ~ transpose(ori(frame_b).R)*(ori(frame_b)*frame_b.r_0) - transpose(ori(frame_a).R)*(ori(frame_a)*frame_a.r_0)
  else
    rRod_0 ~ frame_b.r_0 - frame_a.r_0
  end

  #rRod_0 = frame_b.r_0 - frame_a.r_0;
  rRod_a ~ resolve2(ori(frame_a), rRod_0)
  eRod_a ~ rRod_a/rodlength

  # Constraint equation
  constraint_residue ~ 0

  # Cut-torques at frame_a and frame_b
  frame_a.tau ~ zeros(3)
  frame_b.tau ~ zeros(3)

  #= Force and torque balance of rod
     - Kinematics for center of mass CM of mass point
       r_cm_0 = frame_a.r_0 + rRod_0/2;
       v_cm_0 = der(r_cm_0);
       a_cm_a = resolve2(ori(frame_a), der(v_cm_0) - world.gravity_acceleration(r_cm_0));
     - Inertial and gravity force in direction (f_cm_e) and orthogonal (f_cm_n) to rod
       f_cm_a = m*a_cm_a
       f_cm_e = f_cm_a*eRod_a;           # in direction of rod
       f_cm_n = rodlength(f_cm_a - f_cm_e);  # orthogonal to rod
     - Force balance in direction of rod
       f_cm_e = fa_rod_e + fb_rod_e;
     - Force balance orthogonal to rod
       f_cm_n = fa_rod_n + fb_rod_n;
     - Torque balance with respect to frame_a
       0 = (-f_cm_n)*rodlength/2 + fb_rod_n*rodlength
     The result is:
     fb_rod_n = f_cm_n/2;
     fa_rod_n = fb_rod_n;
     fb_rod_e = f_cm_e - fa_rod_e;
     fa_rod_e is the unknown computed from loop
  =#

    # f_b_a1 is needed in aggregation joints to solve kinematic loops analytically
  if m > 0
    [r_cm_0 ~ frame_a.r_0 + rRod_0/2;
    v_cm_0 ~ D(r_cm_0);
    f_cm_a ~ m*resolve2(ori(frame_a), D(v_cm_0) - gravity_acceleration(r_cm_0))
    f_cm_e ~ (f_cm_a'eRod_a)*eRod_a
    frame_a.f ~ (f_cm_a - f_cm_e)./2 + f_rod*eRod_a
    f_b_a1 ~ (f_cm_a + f_cm_e)./2
    frame_b.f ~ resolve_relative(f_b_a1 - f_rod*eRod_a, ori(frame_a),
      ori(frame_b));
    ]
  else
    [r_cm_0 ~ zeros(3);
    v_cm_0 ~ zeros(3);
    f_cm_a ~ zeros(3);
    f_cm_e ~ zeros(3);
    f_b_a1 ~ zeros(3);
    frame_a.f ~ f_rod*eRod_a;
    frame_b.f ~ -resolve_relative(frame_a.f, ori(frame_a), ori(frame_b));
    ]
  end
    ]


    sys = extend(System(eqs, t; name=:nothing), ptf)
    add_params(sys, pars; name)
end

"""
    SphericalConstraint(; name, color = [1, 1, 0, 1], radius = 0.1, x_locked = true, y_locked = true, z_locked = true)

Spherical cut joint and translational directions may be constrained or released

This model does not use explicit variables e.g. state variables in order to describe the relative motion of `frame_b` with to respect to `frame_a`, but defines kinematic constraints between the `frame_a` and `frame_b`. The forces and torques at both frames are then evaluated in such a way that the constraints are satisfied. Sometimes this type of formulation is also called an implicit joint in literature.

As a consequence of the formulation the relative kinematics between `frame_a` and `frame_b` cannot be initialized.

In complex multibody systems with closed loops this may help to simplify the system of non-linear equations. Please compare state realization chosen by `structural_simplify` using the classical joint formulation and the alternative formulation used here in order to check whether this fact applies to the particular system under consideration. In systems without closed loops the use of this implicit joint is not recommended.

# Arguments
- `x_locked`: Set to false if the translational motion in x-direction shall be free
- `y_locked`: Set to false if the translational motion in y-direction shall be free
- `z_locked`: Set to false if the translational motion in z-direction shall be free

# Rendering parameters
- `color`: Color of the joint in animations (RGBA)
- `radius`: Radius of the joint in animations
"""
@component function SphericalConstraint(; name,
                   color = [1, 1, 0, 1],
                   radius = 0.1,
                   x_locked = true, y_locked = true, z_locked = true,
)

    @named begin
        ptf = PartialTwoFrames()
        Rrel = NumRotationMatrix()
        Rrel_inv = NumRotationMatrix()
    end
    pars = @parameters begin
        radius = radius, [description = "radius of the joint in animations"]
        color[1:4] = color, [description = "color of the joint in animations (RGBA)"]
    end
    @unpack frame_a, frame_b = ptf
    @variables begin (r_rel_a(t)[1:3] = zeros(3)),
                     [
                         description = "Position vector from origin of frame_a to origin of frame_b, resolved in frame_a",
                     ] end

    Rrel = relative_rotation(frame_a, frame_b)

    eqs = [
        if x_locked
            r_rel_a[1] ~ 0
        else
            frame_a.f[1] ~ 0
        end

        if y_locked
            r_rel_a[2] ~ 0
        else
            frame_a.f[2] ~ 0
        end

        if z_locked
            r_rel_a[3] ~ 0
        else
            frame_a.f[3] ~ 0
        end
        r_rel_a ~ resolve2(ori(frame_a), frame_b.r_0 - frame_a.r_0);
        zeros(3) ~ frame_b.tau;
        frame_b.f ~ -resolve2(Rrel, frame_a.f);
        zeros(3) ~ frame_a.tau + resolve1(Rrel, frame_b.tau) - cross(r_rel_a, frame_a.f);
    ]

    sys = extend(System(eqs, t; name=:nothing), ptf)
    add_params(sys, pars; name)
end

"""
    PrismaticConstraint(; name, color, radius = 0.05, x_locked = true, y_locked = true, z_locked = true, render = true)

This model does not use explicit variables e.g. state variables in order to describe the relative motion of `frame_b` with respect to `frame_a`, but defines kinematic constraints between the `frame_a` and `frame_b`. The forces and torques at both frames are then evaluated in such a way that the constraints are satisfied. Sometimes this type of formulation is called an implicit joint in literature.

As a consequence of the formulation, the relative kinematics between `frame_a` and `frame_b` cannot be initialized.

In complex multibody systems with closed loops this may help to simplify the system of non-linear equations. Compare the simplification result using the classical joint formulation and this alternative formulation to check which one is more efficient for the particular system under consideration.

In systems without closed loops the use of this implicit joint does not make sense or may even be disadvantageous.

# Parameters
- `color`: Color of the joint in animations (RGBA)
- `radius`: Radius of the joint in animations
- `x_locked`: Set to false if the translational motion in x-direction shall be free
- `y_locked`: Set to false if the translational motion in y-direction shall be free
- `z_locked`: Set to false if the translational motion in z-direction shall be free
- `render`: Whether or not the joint is rendered in animations
"""
@component function PrismaticConstraint(; name, color = [1, 1, 0, 1], radius = 0.05, x_locked = true, y_locked = true, z_locked = true, render = true)
    @named begin
        ptf = PartialTwoFrames()
    end
    pars = @parameters begin
        radius = radius, [description = "radius of the joint in animations"]
        color[1:4] = color, [description = "color of the joint in animations (RGBA)"]
        render = render, [description = "Set to false if the joint shall not be rendered"]
    end
    @unpack frame_a, frame_b = ptf
    @variables (r_rel_a(t)[1:3] = zeros(3)), [description = "Position vector from origin of frame_a to origin of frame_b, resolved in frame_a"]

    Rrel = relative_rotation(frame_a, frame_b)

    eqs = [
        if x_locked
            r_rel_a[1] ~ 0
        else
            frame_a.f[1] ~ 0
        end

        if y_locked
            r_rel_a[2] ~ 0
        else
            frame_a.f[2] ~ 0
        end

        if z_locked
            r_rel_a[3] ~ 0
        else
            frame_a.f[3] ~ 0
        end
        r_rel_a ~ resolve2(ori(frame_a), frame_b.r_0 - frame_a.r_0)
        zeros(3) ~ frame_a.tau + resolve1(Rrel, frame_b.tau) + cross(r_rel_a, resolve1(Rrel, frame_b.f))
        zeros(3) ~ resolve1(Rrel, frame_b.f) + frame_a.f
        # orientation_constraint(ori(frame_a), ori(frame_b)) ~ 0
        residue(ori(frame_a), ori(frame_b)) ~ 0
    ]

    Main.eqs = eqs

    sys = extend(System(eqs, t, collect(r_rel_a), pars; name), ptf)
end


"""
    UniversalSpherical(; name, n1_a, rRod_ia, sphere_diameter = 0.1, sphere_color, rod_width = 0.1, rod_height = 0.1, rod_color, cylinder_length = 0.1, cylinder_diameter = 0.1, cylinder_color, kinematic_constraint = true)

Universal - spherical joint aggregation (1 constraint, no potential states)

This component consists of a universal joint at `frame_a` and a spherical joint at `frame_b` that are connected together with a rigid rod.

This joint aggregation has no mass and no inertia and introduces the constraint that the distance between the origin of `frame_a` and the origin of `frame_b` is constant (= `length(rRod_ia)`). The universal joint is defined in the following way:
- The rotation axis of revolute joint 1 is along parameter vector `n1_a` which is fixed in `frame_a`.
- The rotation axis of revolute joint 2 is perpendicular to axis 1 and to the line connecting the universal and the spherical joint.

Note, there is a singularity when axis 1 and the connecting rod are parallel to each other. Therefore, if possible `n1_a` should be selected in such a way that it is perpendicular to `rRod_ia` in the initial configuration (i.e., the distance to the singularity is as large as possible).

An additional `frame_ia` is present. It is fixed in the connecting rod at the origin of `frame_a`. The placement of `frame_ia` on the rod is implicitly defined by the universal joint (`frame_a` and `frame_ia` coincide when the angles of the two revolute joints of the universal joint are zero) and by parameter vector `rRod_ia`, the position vector from the origin of `frame_a` to the origin of `frame_b`, resolved in `frame_ia`.

This joint aggregation can be used in cases where in reality a rod with spherical joints at end are present. Such a system has an additional degree of freedom to rotate the rod along its axis. In practice this rotation is usually of no interest and is mathematically removed by replacing one of the spherical joints by a universal joint. Still, in most cases the [`SphericalSpherical`](@ref) joint aggregation can be used instead of the UniversalSpherical joint since the rod is animated and its mass properties are approximated by a point mass in the middle of the rod. The [`SphericalSpherical`](@ref) joint has the advantage that it does not have a singular configuration.

# Arguments
- `n1_a` Axis 1 of universal joint resolved in frame_a (axis 2 is orthogonal to axis 1 and to rod)
- `rRod_ia` Vector from origin of frame_a to origin of frame_b, resolved in `frame_ia` (if computeRodLength=true, rRod_ia is only an axis vector along the connecting rod)
- `kinematic_constraint = true` Set to false if no constraint shall be defined, due to analytically solving a kinematic loop
- `constraint_residue` If set to `:external`, an equation in the parent system is expected to define this variable, e.g., `rod.constraint_residue ~ rod.f_rod - f_rod` where `rod` is the name of the UniversalSpherical joint. If unspecified, the length constraint `rRod_0'rRod_0 - rodLength'rodLength` is used

# Connectors
- `frame_a`: Frame for the universal joint
- `frame_b`: Frame for the spherical joint
- `frame_ia`: Frame fixed in the rod at the origin of `frame_a`

# Rendering parameters
- `sphere_diameter`: Diameter of spheres representing the universal and the spherical joint
- `sphere_color`: Color of spheres representing the universal and the spherical joint (RGBA)
- `rod_width`: Width of rod shape in direction of axis 2 of universal joint
- `rod_height`: Height of rod shape in direction that is orthogonal to rod and to axis 2
- `rod_color`: Color of rod shape connecting the universal and the spherical joints (RGBA)
- `cylinder_length`: Length of cylinders representing the two universal joint axes
- `cylinder_diameter`: Diameter of cylinders representing the two universal joint axes
- `cylinder_color`: Color of cylinders representing the two universal joint axes (RGBA)
"""
@component function UniversalSpherical(; name,
                    n1_a = [0, 0, 1],
                    rRod_ia = [1, 0, 0],
                    sphere_diameter = 0.1,
                    sphere_color = [1, 0.2, 1, 0.9],
                    rod_width = 0.1,
                    rod_height = 0.1,
                    rod_color = [0, 0.1, 1, 0.9],
                    constraint_residue = nothing,
                    cylinder_length = 0.1,
                    cylinder_diameter = 0.1,
                    cylinder_color = [1, 0.2, 0, 1], 
                    kinematic_constraint = true,
                    render = true,
)

    residue = constraint_residue
    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
        frame_ia = Frame()
        # R_rel_ia_frame = Frame(varw=true)
    end
    pars = @parameters begin
        sphere_diameter = sphere_diameter, [description = "Diameter of spheres representing the universal and the spherical joint"]
        sphere_color[1:4] = sphere_color, [description = "Color of spheres representing the universal and the spherical joint (RGBA)"]
        rod_width = rod_width, [description = "Width of rod shape in direction of axis 2 of universal joint"]
        rod_height = rod_height, [description = "Height of rod shape in direction that is orthogonal to rod and to axis 2"]
        rod_color[1:4] = rod_color, [description = "Color of rod shape connecting the universal and the spherical joints (RGBA)"]
        cylinder_length = cylinder_length, [description = "Length of cylinders representing the two universal joint axes"]
        cylinder_diameter = cylinder_diameter, [description = "Diameter of cylinders representing the two universal joint axes"]
        cylinder_color[1:4] = cylinder_color, [description = "Color of cylinders representing the two universal joint axes (RGBA)"]
        n1_a[1:3] = n1_a, [description = "Axis 1 of universal joint resolved in frame_a (axis 2 is orthogonal to axis 1 and to rod)"]
        rRod_ia[1:3] = rRod_ia, [description = "Vector from origin of frame_a to origin of frame_b, resolved in frame_ia (if computeRodLength=true, rRod_ia is only an axis vector along the connecting rod)"]
        rodLength = _norm(rRod_ia), [description="Length of rod (distance between origin of frame_a and origin of frame_b)"]
        render = render, [description="Whether or not to render the joint in animations"]
    end

    # pars = collect_all(pars)

    vars = @variables begin
        f_rod(t), [description="Constraint force in direction of the rod (positive on frame_a, when directed from frame_a to frame_b)"]
        # eRod_ia(t)[1:3], [description="Unit vector from origin of frame_a to origin of frame_b, resolved in frame_ia"]
        # e2_ia(t)[1:3], [description="Unit vector in direction of axis 2 of universal joint, resolved in frame_ia (orthogonal to n1_a and eRod_ia; note: frame_ia is parallel to frame_a when the universal joint angles are zero)"]
        # e3_ia(t)[1:3], [description="Unit vector perpendicular to eRod_ia and e2_ia, resolved in frame_ia"]
        f_b_a1(t)[1:3], [description="frame_b.f without f_rod part, resolved in frame_a (needed for analytic loop handling)"]
        eRod_a(t)[1:3], [description="Unit vector in direction of rRod_a, resolved in frame_a (needed for analytic loop handling)"]
        rRod_0(t)[1:3], [guess = rRod_ia, description="Position vector from frame_a to frame_b resolved in world frame"]
        rRod_a(t)[1:3], [guess = rRod_ia, description="Position vector from frame_a to frame_b resolved in frame_a"]
        (constraint_residue(t)), [guess = 0, description="Constraint equation of joint in residue form: Either length constraint (= default) or equation to compute rod force (for analytic solution of loops in combination with Internal.RevoluteWithLengthConstraint/PrismaticWithLengthConstraint)"]
        f_b_a(t)[1:3], [description="frame_b.f resolved in frame_a"]
        f_ia_a(t)[1:3], [description="frame_ia.f resolved in frame_a"]
        t_ia_a(t)[1:3], [description="frame_ia.t resolved in frame_a"]
        n2_a(t)[1:3], [description="Vector in direction of axis 2 of the universal joint (e2_ia), resolved in frame_a"]
        (length2_n2_a(t)), [guess = 1, description="Square of length of vector n2_a"]
        length_n2_a(t), [description="Length of vector n2_a"]
        e2_a(t)[1:3], [description="Unit vector in direction of axis 2 of the universal joint (e2_ia), resolved in frame_a"]
        e3_a(t)[1:3], [description="Unit vector perpendicular to eRod_ia and e2_a, resolved in frame_a"]
        der_rRod_a_L(t)[1:3], [description="= der(rRod_a)/rodLength"]
        w_rel_ia1(t)[1:3], [description="Angular velocity of frame_ia1 with respect to frame_a"]
        # R_rel_ia1(t), [description="Rotation object from frame_a to frame_ia1 (frame that is fixed in frame_ia such that x-axis is along the rod axis)"]
        # R_rel_ia2(t), [description="Fixed rotation object from frame_ia1 to frame_ia"]
        # R_rel_ia(t), [description="Rotation from frame_a to frame_ia"]
    end

    eRod_ia = normalize(rRod_ia)
    e2_ia = cross(n1_a, eRod_ia)
    e3_ia = cross(eRod_ia, e2_ia)


    R_rel_ia1 = RotationMatrix(transpose([eRod_a e2_a e3_a]), w_rel_ia1)
    R_rel_ia2 = RotationMatrix([eRod_ia e2_ia e3_ia], zeros(3))
    # R_rel_ia = ori(R_rel_ia_frame, true)
    # R_rel_ia = ori(R_rel_ia_frame)
    R_rel_ia = absolute_rotation(R_rel_ia1, R_rel_ia2)

    Ra = ori(frame_a)
    
    Ria = absolute_rotation(Ra, R_rel_ia)

    eqs = [

        if kinematic_constraint
            rRod_0 ~ ori(frame_b).R'*(ori(frame_b).R*frame_b.r_0) - Ra.R'*(Ra.R*frame_a.r_0)
        else
            rRod_0 ~ frame_b.r_0 - frame_a.r_0
        end

        rRod_a ~ resolve2(Ra, rRod_0)

        constraint_residue ~ 0

        eRod_a ~ rRod_a/rodLength
        n2_a ~ cross(n1_a, eRod_a)
        length2_n2_a ~ dot(n2_a,n2_a)

        # assert(length2_n2_a > 1e-10, "A MultiBody.Joints.UniversalSpherical joint (consisting of a universal joint and a spherical joint connected together by a rigid rod) is in the singular configuration of the universal joint. This means that axis 1 of the universal joint defined via parameter \"n1_a\" is parallel to vector \"rRod_ia\" that is directed from the origin of frame_a to the origin of frame_b. You may try to use another \"n1_a\" vector. If this fails, use instead MultiBody.Joints.SphericalSpherical, if this is possible, because this joint aggregation does not have a singular configuration.")

        length_n2_a ~ sqrt(length2_n2_a)
        e2_a ~ n2_a/length_n2_a
        e3_a ~ cross(eRod_a, e2_a)

        der_rRod_a_L ~ (resolve2(Ra, D(rRod_0)) - cross(Ra.w, rRod_a))/rodLength
        w_rel_ia1 ~ [e3_a'cross(n1_a, der_rRod_a_L)/length_n2_a, -dot(e3_a,der_rRod_a_L), dot(e2_a,der_rRod_a_L)]

        # R_rel_ia ~ Rrelia# absolute_rotation(R_rel_ia1, R_rel_ia2)
        # R_rel_ia.w ~ Rrelia.w

        frame_ia.r_0 ~ frame_a.r_0
        ori(frame_ia) ~ Ria # absolute_rotation(frame_a, R_rel_ia)
        # ori(frame_ia).w ~ Ria.w

        f_ia_a ~ resolve1(R_rel_ia, frame_ia.f)
        t_ia_a ~ resolve1(R_rel_ia, frame_ia.tau)

        f_b_a1 ~ -e2_a*(dot(n1_a,t_ia_a)/(rodLength*dot(n1_a,e3_a))) + collect(e3_a)*(dot(e2_a,t_ia_a)/rodLength)
        f_b_a ~ -f_rod*eRod_a + f_b_a1
        frame_b.f ~ resolve_relative(f_b_a, Ra, ori(frame_b))
        frame_b.tau ~ zeros(3)
        zeros(3) ~ frame_a.f + f_b_a + f_ia_a
        zeros(3) ~ frame_a.tau + t_ia_a + cross(rRod_a, f_b_a)
    ]

    if residue === nothing
        push!(eqs, constraint_residue ~ dot(rRod_0,rRod_0) - dot(rodLength,rodLength))
    else
        # See implementation of JointUSR for how to use the constraint_residue=:external
        residue === :external || error("Unknown value for constraint_residue, expected nothing or :external")
    end

    sys = System(eqs, t; name=:nothing, systems)
    add_params(sys, pars; name)
end


@component function RevoluteWithLengthConstraint(; name, n = Float64[0, 0, 1], axisflange = false,
    positive_branch = true, radius = 0.05, length = radius, color = [0.5019608f0,0.0f0,0.5019608f0,1.0f0], state_priority = 1.0, phi_offset=0, phi_guess=0, length_constraint=1, use_arrays = false)
    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
        axis = Rotational.Flange()
        bearing = Rotational.Flange()
        position_a = RealInput(nin=3) # Position vector from frame_a to frame_a side of length constraint, resolved in frame_a of revolute joint
        position_b = RealInput(nin=3)
    end
    @parameters e[1:3] = _normalize(n) [description = "normalized axis of rotation"]
    # @parameters n[1:3]=n [description = "axis of rotation"] # Can't have this as parameter since e = _normalize(n) does not work :/
    @parameters phi_offset = phi_offset, [description = "offset of the joint in animations"]
    @parameters length_constraint = length_constraint, [description = "Fixed length of length constraint"]
    pars = @parameters begin
        radius = radius, [description = "radius of the joint in animations"]
        length = length, [description = "length of the joint in animations"]
        color[1:4] = color, [description = "color of the joint in animations (RGBA)"]
    end
    @variables tau(t) [
        connect = Flow,
        description = "Driving torque in direction of axis of rotation",
    ]
    @variables phi(t) [
        state_priority = 1,
        description = "Relative rotation angle from frame_a to frame_b",
    ]
    @variables angle(t) [
        state_priority = -1,
        description = "= phi + phi_offset (relative rotation angle between frame_a and frame_b)",
    ]
    @variables r_a(t)[1:3], [description = "Position vector from frame_a to frame_a side of length constraint, resolved in frame_a of revolute joint"]
    @variables r_b(t)[1:3], [description = "Position vector from frame_b to frame_b side of length constraint, resolved in frame_b of revolute joint"]

    # vars = [tau; phi; angle; r_a; r_b]
    # pars = [collect(e); phi_offset; length_constraint]

    # @parameters positive_branch::Bool=false
    # NOTE: final parameters in modelica can be implemented by parameter_dependencies = [final_parameter => expression with other parameters]

    Rrel = planar_rotation(e, angle, D(angle))
    Rb = absolute_rotation(ori(frame_a), Rrel)

    # IR = JuliaSimCompiler
    # e_array = IR.make_array((3,), e...)
    # r_a_array = IR.make_array((3,), r_a...)
    # r_b_array = IR.make_array((3,), r_b...)

    eqs = [
        r_a ~ position_a.u
        r_b ~ position_b.u

        axis.tau ~ tau
        axis.phi ~ phi
        bearing.phi ~ 0
        angle ~ phi + phi_offset
        frame_b.r_0 ~ frame_a.r_0

        ori(frame_b) ~ Rb
        # ori(frame_b).w ~ Rb.w

        zeros(3) ~ frame_a.f + resolve1(Rrel, frame_b.f)
        zeros(3) ~ frame_a.tau + resolve1(Rrel, frame_b.tau)

        if use_arrays
            angle ~ compute_angle2(length_constraint, e, r_a, r_b, positive_branch)[1]
            # angle ~ Symbolics.term(compute_angle2, length_constraint, e, r_a, r_b, positive_branch, type=Real)
        else
            # angle ~ Symbolics.term(compute_angle, length_constraint, e..., r_a..., r_b..., positive_branch, type=Real)
            angle ~ compute_angle(length_constraint, e, r_a, r_b, positive_branch)
        end
    ]

    sys = System(eqs, t; name=:nothing, systems)#, parameter_dependencies = [positive_branch => select_branch(length_constraint, e, phi_offset + phi_guess, r_a, r_b)])  # JuliaSimCompiler ignores parameter dependencies, the user has to provide it instead
    
    add_params(sys, pars; name)
end


"""
    JointUSR(;
        name,
        n1_a = [0, 0, 1],
        n_b = [0, 0, 1],
        rRod1_ia = [1, 0, 0],
        rRod1_ib = [-1, 0, 0],
        phi_offset = 0,
        phi_guess = 0,
    )

This component consists of a universal joint at `frame_a`, a revolute joint at `frame_b` and a spherical joint which is connected via rod1 to the universal and via rod2 to the revolute joint.

This joint aggregation has no mass and no inertia and introduces neither constraints nor potential state variables. It should be used in kinematic loops whenever possible since the non-linear system of equations introduced by this joint aggregation is solved analytically (i.e., a solution is always computed, if a unique solution exists).

The universal joint is defined in the following way:

- The rotation axis of revolute joint 1 is along parameter vector `n1_a` which is fixed in `frame_a`.
- The rotation axis of revolute joint 2 is perpendicular to axis 1 and to the line connecting the universal and the spherical joint (= rod 1).

The definition of axis 2 of the universal joint is performed according to the most often occurring case for the sake of simplicity. Otherwise, the treatment is much more complicated and the number of operations is considerably higher, if axis 2 is not orthogonal to axis 1 and to the connecting rod.

Note, there is a singularity when axis 1 and the connecting rod are parallel to each other. Therefore, if possible `n1_a` should be selected in such a way that it is perpendicular to rRod1_ia in the initial configuration (i.e., the distance to the singularity is as large as possible).

The rest of this joint aggregation is defined by the following parameters:

- `positive_branch`: The positive branch of the revolute joint is selected (cf. elbow up vs. elbow down).
- The position of the spherical joint with respect to the universal joint is defined by vector `rRod1_ia`. This vector is directed from `frame_a` to the spherical joint and is resolved in `frame_ia` (it is most simple to select `frame_ia` such that it is parallel to `frame_a` in the reference or initial configuration).
- The position of the spherical joint with respect to the revolute joint is defined by vector `rRod2_ib`. This vector is directed from the inner frame of the revolute joint (`frame_ib` or `revolute.frame_a`) to the spherical joint and is resolved in `frame_ib` (note, that `frame_ib` and `frame_b` are parallel to each other).
- The axis of rotation of the revolute joint is defined by axis vector `n_b`. It is fixed and resolved in `frame_b`.
- When specifying this joint aggregation with the definitions above, two different configurations are possible. Via parameter `phi_guess` a guess value for `revolute.phi(t0)` at the initial time `t0` is given. The configuration is selected that is closest to `phi_guess` (`|revolute.phi - phi_guess|` is minimal).

# Connectors
- `frame_a`: Frame for the universal joint
- `frame_b`: Frame for the revolute joint
- An additional `frame_ia` is present. It is fixed in the rod connecting the universal and the spherical joint at the origin of `frame_a`. The placement of `frame_ia` on the rod is implicitly defined by the universal joint (`frame_a` and `frame_ia` coincide when the angles of the two revolute joints of the universal joint are zero) and by parameter vector rRod1_ia, the position vector from the origin of `frame_a` to the spherical joint, resolved in `frame_ia`.
- An additional `frame_ib` is present. It is fixed in the rod connecting the revolute and the spherical joint at the side of the revolute joint that is connected to this rod (`= rod2.frame_a = revolute.frame_a`).
- An additional `frame_im` is present. It is fixed in the rod connecting the revolute and the spherical joint at the side of the spherical joint that is connected to this rod (`= rod2.frame_b`). It is always parallel to `frame_ib`.
"""
@component function JointUSR(;
    name,
    n1_a = [0, 0, 1],
    n_b = [0, 0, 1],
    rRod1_ia = [1, 0, 0],
    rRod2_ib = [-1, 0, 0],
    rod_color = purple,
    rod_radius = 0.05,
    phi_offset = 0,
    phi_guess = 0,
    positive_branch,
    use_arrays = false,
    render = true,
)
    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
        frame_ia = Frame()
        frame_ib = Frame()
        frame_im = Frame()
        axis = Rotational.Flange()
        bearing = Rotational.Flange()
    end

    @parameters begin
        # n1_a[1:3] = n1_a, [description = "Axis 1 of universal joint fixed and resolved in frame_a (axis 2 is orthogonal to axis 1 and to rod 1)"]
        # n_b[1:3] = n_b, [description = "Axis of revolute joint fixed and resolved in frame_b"]
        rRod1_ia[1:3] = rRod1_ia, [description = "Vector from origin of frame_a to spherical joint, resolved in frame_ia"]
        # rRod2_ib[1:3] = rRod2_ib, [description = "Vector from origin of frame_ib to spherical joint, resolved in frame_ib"]
        # phi_offset = phi_offset, [description = "Relative angle offset of revolute joint (angle = phi(t) + from_deg(phi_offset))"]
        phi_guess = phi_guess, [description = "Select the configuration such that at initial time |phi(t0) - from_deg(phi_guess)| is minimal"]
        # render = render, [description = "Whether or not to render the joint in animations"]
    end



    @variables begin
        aux(t), [description = "Denominator used to compute force in rod connecting universal and spherical joint"]
        f_rod(t), [description = "Constraint force in direction of the rod (positive, if rod is pressed)"]
    end

    more_systems = @named begin
        rod1 = UniversalSpherical(;
            n1_a,
            rRod_ia = rRod1_ia,
            kinematic_constraint = false,
            constraint_residue = :external,
            sphere_color = [0,0,0,0],
            rod_color = rod_color,
            cylinder_color = rod_color,
            cylinder_diameter = 2*rod_radius,
            rod_width = 2*rod_radius,
            rod_height = 2*rod_radius,
            render,
        )
        rod2 = FixedTranslation(;
            r = rRod2_ib,
            color = rod_color,
            radius = rod_radius,
            render,
        )
        position_b = Constant3(k = rRod2_ib)
        relative_position = RelativePosition(resolve_frame = :frame_a)
    end

    @named revolute = RevoluteWithLengthConstraint(;
        n = n_b,
        length_constraint = rod1.rodLength,
        phi_offset,
        phi_guess,
        positive_branch,
        use_arrays,
    )
    push!(more_systems, revolute)

    eqs = [
        aux ~ cross(revolute.e, rRod2_ib)'resolve_relative(rod1.eRod_a, ori(rod1.frame_a), ori(rod1.frame_b))
        f_rod ~ (
            -revolute.tau - dot(revolute.e, (frame_ib.tau + frame_im.tau +
                cross(rRod2_ib, frame_im.f) -
                cross(rRod2_ib, resolve_relative(rod1.f_b_a1, ori(rod1.frame_a), ori(rod1.frame_b)))
            ))
        )/max(abs(aux), 1e-10)

        rod1.constraint_residue ~ rod1.f_rod - f_rod # Externally provided residue
        connect(revolute.frame_b, rod2.frame_a)
        connect(rod2.frame_b, rod1.frame_b)
        connect(revolute.frame_a, frame_b)
        connect(rod2.frame_a, frame_ib)
        connect(rod1.frame_a, frame_a)
        connect(relative_position.frame_b, frame_a)
        connect(relative_position.frame_a, frame_b)
        connect(position_b.output, revolute.position_b)
        connect(rod2.frame_b, frame_im)
        connect(rod1.frame_ia, frame_ia)
        connect(revolute.axis, axis)
        connect(relative_position.r_rel, revolute.position_a)
        connect(revolute.bearing, bearing)
    ]
    System(eqs, t; name, systems=[systems; more_systems])

end

@component function Constant3(; name, k = zeros(3))
    pars = @parameters begin
        k[1:3] = k, [description = "Constant output value of block"]
    end

    systems = @named begin
        output = Blocks.RealOutput(nout=3)
    end

    vars = @variables begin
    end

    equations = Equation[
        output.u[1] ~ k[1]
        output.u[2] ~ k[2]
        output.u[3] ~ k[3]
    ]

    return System(equations, t; name, systems)
end

"""
    JointRRR(;
        name,
        n_a = [0,0,1],
        n_b = [0,0,1],
        rRod1_ia = [1,0,0],
        rRod2_ib = [-1,0,0],
        phi_offset = 0, 
        phi_guess = 0,

    )   

This component consists of 3 revolute joints with parallel axes of rotation that are connected together by two rods.

This joint aggregation introduces neither constraints nor state variables and should therefore be used in kinematic loops whenever possible to avoid non-linear systems of equations. It is only meaningful to use this component in planar loops. Basically, the position and orientation of the 3 revolute joints as well as of `frame_ia`, `frame_ib`, and `frame_im` are calculated by solving analytically a non-linear equation, given the position and orientation at `frame_a` and at `frame_b`.

Connector `frame_a` is the "left" side of the first revolute joint whereas `frame_ia` is the "right side of this revolute joint, fixed in rod 1. Connector `frame_b` is the "right" side of the third revolute joint whereas `frame_ib` is the "left" side of this revolute joint, fixed in rod 2. Finally, connector `frame_im` is the connector at the "right" side of the revolute joint in the middle, fixed in rod 2.

The easiest way to define the parameters of this joint is by moving the MultiBody system in a reference configuration where all frames of all components are parallel to each other (alternatively, at least `frame_a`, `frame_ia`, `frame_im`, `frame_ib`, `frame_b` of the JointRRR joint should be parallel to each other when defining an instance of this component).

Basically, the JointRRR model internally consists of a universal-spherical-revolute joint aggregation (= JointUSR). In a planar loop this will behave as if 3 revolute joints with parallel axes are connected by rigid rods.

# Arguments
- `n_a` Axis of revolute joints resolved in `frame_a` (all axes are parallel to each other)
- `n_b` Axis of revolute joint fixed and resolved in `frame_b`
- `rRod1_ia` Vector from origin of `frame_a` to revolute joint in the middle, resolved in `frame_ia`
- `rRod2_ib` Vector from origin of `frame_ib` to revolute joint in the middle, resolved in `frame_ib`
- `phi_offset` Relative angle offset of revolute joint at `frame_b` `(angle = phi(t) + phi_offset)`

# Connectors
- `frame_a`: Coordinate system fixed to the component with one cut-force and cut-torque
- `frame_b`: Coordinate system fixed to the component with one cut-force and cut-torque
- `frame_ia`: Coordinate system at origin of `frame_a` fixed at connecting rod of left and middle revolute joint
- `frame_ib`: Coordinate system at origin of `frame_ib` fixed at connecting rod of middle and right revolute joint
- `frame_im`: Coordinate system at origin of revolute joint in the middle fixed at connecting rod of middle and right revolute joint
- `axis`: 1-dim. rotational flange that drives the right revolute joint at `frame_b`
- `bearing`: 1-dim. rotational flange of the drive bearing of the right revolute joint at `frame_b`
"""
@component function JointRRR(;
    name,
    n_a = [0,0,1],
    n_b = [0,0,1],
    rRod1_ia = [1,0,0],
    rRod2_ib = [-1,0,0],
    phi_offset = 0, 
    phi_guess = 0,
    positive_branch = true,
    kwargs...
)

    @parameters begin
        # n_a[1:3] = n_a, [description = "Axes of revolute joints resolved in frame_a (all axes are parallel to each other)"]
        # n_b[1:3] = n_b, [description = "Axis of revolute joint fixed and resolved in frame_b"]
        rRod1_ia[1:3] = rRod1_ia, [description = "Vector from origin of frame_a to revolute joint in the middle, resolved in frame_ia"]
        # rRod2_ib[1:3] = rRod2_ib, [description = "Vector from origin of frame_ib to revolute joint in the middle, resolved in frame_ib"]
        # phi_offset = phi_offset, [description = "Relative angle offset of revolute joint at frame_b (angle = phi(t) + from_deg(phi_offset))"]
        phi_guess = phi_guess, [description = "Select the configuration such that at initial time |phi(t0) - from_deg(phi_guess)| is minimal"]

    end
    systems = @named begin
        frame_a = Frame()
        frame_b = Frame()
        frame_ia = Frame()
        frame_ib = Frame()
        frame_im = Frame()
        axis = Rotational.Flange()
        bearing = Rotational.Flange()
        jointUSR = JointUSR(;
            n1_a = n_a,
            n_b = n_b,
            phi_offset,
            phi_guess,
            rRod2_ib,
            rRod1_ia,
            positive_branch,
            kwargs...
        )
    end

    # e_a = _normalize(n_a)
    # e_ia = jointUSR.e2_ia
    # e_b = jointUSR.revolute.e

    eqs = [
        connect(jointUSR.frame_a, frame_a)
        connect(jointUSR.frame_b, frame_b)
        connect(jointUSR.frame_ia, frame_ia)
        connect(jointUSR.frame_im, frame_im)
        connect(jointUSR.frame_ib, frame_ib)
        connect(jointUSR.axis, axis)
        connect(jointUSR.bearing, bearing)
    ]

    System(eqs, t; name, systems)
end


function select_branch(L::Real, e::AbstractVector, angle_guess::Real, r_a::AbstractVector, r_b::AbstractVector)
    e_r_a = e'r_a
    e_r_b = e'r_b
    A = -2*(r_b'r_a - e_r_b'e_r_a)
    B = 2*r_b'cross(e, r_a)
    C = r_a'r_a + r_b'r_b - L^2 - 2*e_r_b'e_r_a
    k1 = A^2 + B^2
    k1a = k1 - C^2
    k1a > 1e-10 || error("Singular position of loop (either no or two analytic solutions; the mechanism has lost one-degree-of freedom in this position). Try to use another joint-assembly component. In most cases it is best to let joints outside of the JointXXX component be revolute and _not_ prismatic joints. If this also leads to singular positions, it could be that this kinematic loop cannot be solved analytically. In this case you have to build up the loop with basic joints (_no_ aggregation JointXXX components) and rely on dynamic state selection, i.e., during simulation the states will be dynamically selected in such a way that in no position a degree of freedom is lost.")
    k1b = max(k1a, 1.0e-12)
    k2 = sqrt(k1b)
    kcos1 = -A*C + B*k2
    ksin1 = -B*C - A*k2
    angle1 = atan(ksin1, kcos1)
    kcos2 = -A*C - B*k2
    ksin2 = -B*C + A*k2
    angle2 = atan(ksin2, kcos2)
    positive_branch = abs(angle1 - angle_guess) <= abs(angle2 - angle_guess)
    [positive_branch]
end

@register_array_symbolic select_branch(L::Real, e::AbstractVector, angle_guess::Real, r_a::AbstractVector, r_b::AbstractVector) begin
    size = (1, )
    eltype = Bool
end

function compute_angle(L::Real, e, r_a, r_b, positive_branch)
    e_r_a = e'r_a
    e_r_b = e'r_b
    A = -2*(r_b'r_a - e_r_b'e_r_a)
    B = 2*r_b'cross(e, r_a)
    C = r_a'r_a + r_b'r_b - L^2 - 2*e_r_b'e_r_a
    k1 = A^2 + B^2
    k1a = k1 - C^2
    if k1a isa AbstractFloat
        # k1a > 1e-10 || error("Singular position of loop (either no or two analytic solutions; the mechanism has lost one-degree-of freedom in this position). Try to use another joint-assembly component. In most cases it is best to let joints outside of the JointXXX component be revolute and _not_ prismatic joints. If this also leads to singular positions, it could be that this kinematic loop cannot be solved analytically. In this case you have to build up the loop with basic joints (_no_ aggregation JointXXX components) and rely on dynamic state selection, i.e., during simulation the states will be dynamically selected in such a way that in no position a degree of freedom is lost.")
    end
    k1b = max(k1a, 1.0e-12)
    k2 = sqrt(k1b)
    kcos1 = -A*C + B*k2*ifelse(positive_branch == true, 1, -1)
    ksin1 = -B*C + A*k2*ifelse(positive_branch == true, -1, 1)
    atan(ksin1, kcos1)
end

function compute_angle2(L::Real, e::AbstractVector, r_a::AbstractVector, r_b::AbstractVector, positive_branch)
    e_r_a = e'r_a
    e_r_b = e'r_b
    A = -2*(r_b'r_a - e_r_b'e_r_a)
    B = 2*r_b'cross(e, r_a)
    C = r_a'r_a + r_b'r_b - L^2 - 2*e_r_b'e_r_a
    k1 = A^2 + B^2
    k1a = k1 - C^2
    # k1a > 1e-10 || error("Singular position of loop (either no or two analytic solutions; the mechanism has lost one-degree-of freedom in this position). Try to use another joint-assembly component. In most cases it is best to let joints outside of the JointXXX component be revolute and _not_ prismatic joints. If this also leads to singular positions, it could be that this kinematic loop cannot be solved analytically. In this case you have to build up the loop with basic joints (_no_ aggregation JointXXX components) and rely on dynamic state selection, i.e., during simulation the states will be dynamically selected in such a way that in no position a degree of freedom is lost.")
    k1b = max(k1a, 1.0e-12)
    k2 = sqrt(k1b)
    kcos1 = -A*C + B*k2*ifelse(positive_branch == true, 1, -1)
    ksin1 = -B*C + A*k2*ifelse(positive_branch == true, -1, 1)
    [atan(ksin1, kcos1)]
end

@register_array_symbolic compute_angle2(L::Real, e::AbstractVector, r_a::AbstractVector, r_b::AbstractVector, positive_branch::Bool) begin
    size = (1, )
    eltype = Float64
end

# @register_symbolic compute_angle(L::Real, e::AbstractVector, r_a::AbstractVector, r_b::AbstractVector, positive_branch::Bool)::Real

# @register_symbolic compute_angle(L::Num, e1::Num, e1::Num, e2::Num, r_a1::Num, r_a2::Num, r_a3::Num, r_b1::Num, r_b2::Num, r_b3::Num, positive_branch)

function compute_angle(L::Real, e1::Real, e2::Real, e3::Real, r_a1::Real, r_a2::Real, r_a3::Real, r_b1::Real, r_b2::Real, r_b3::Real, positive_branch)
    e = SA[e1, e2, e3]
    r_a = SA[r_a1, r_a2, r_a3]
    r_b = SA[r_b1, r_b2, r_b3]
    compute_angle(L, e, r_a, r_b, positive_branch)
end
