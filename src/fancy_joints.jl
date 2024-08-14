
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
                    r_0 = [0,0,0],
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

    rRod_0 = collect(rRod_0)
    rRod_a = collect(rRod_a)
    eRod_a = collect(eRod_a)
    r_cm_0 = collect(r_cm_0)
    v_cm_0 = collect(v_cm_0)
    f_cm_a = collect(f_cm_a)
    f_cm_e = collect(f_cm_e)
    f_b_a1 = collect(f_b_a1)
    rodlength = _norm(r_0)
    constraint_residue = rRod_0'rRod_0 - rodlength^2


    eqs = [
        # Determine relative position vector between the two frames
  if kinematic_constraint
    rRod_0 .~ transpose(ori(frame_b).R)*(ori(frame_b)*collect(frame_b.r_0)) - transpose(ori(frame_a).R)*(ori(frame_a)*collect(frame_a.r_0))
  else
    rRod_0 .~ frame_b.r_0 - frame_a.r_0
  end

  #rRod_0 = frame_b.r_0 - frame_a.r_0;
  rRod_a .~ resolve2(ori(frame_a), rRod_0)
  eRod_a .~ rRod_a/rodlength

  # Constraint equation
  constraint_residue ~ 0

  # Cut-torques at frame_a and frame_b
  frame_a.tau .~ zeros(3)
  frame_b.tau .~ zeros(3)

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
    [r_cm_0 .~ frame_a.r_0 + rRod_0/2;
    v_cm_0 .~ D.(r_cm_0);
    f_cm_a .~ m*resolve2(ori(frame_a), D.(v_cm_0) - gravity_acceleration(r_cm_0))
    f_cm_e .~ (f_cm_a'eRod_a)*eRod_a
    frame_a.f .~ (f_cm_a - f_cm_e)./2 + f_rod*eRod_a
    f_b_a1 .~ (f_cm_a + f_cm_e)./2
    frame_b.f .~ resolve_relative(f_b_a1 - f_rod*eRod_a, ori(frame_a),
      ori(frame_b));
    ]
  else
    [r_cm_0 .~ zeros(3);
    v_cm_0 .~ zeros(3);
    f_cm_a .~ zeros(3);
    f_cm_e .~ zeros(3);
    f_b_a1 .~ zeros(3);
    frame_a.f .~ f_rod*eRod_a;
    frame_b.f .~ -resolve_relative(frame_a.f, ori(frame_a), ori(frame_b));
    ]
  end
    ]


    sys = extend(ODESystem(eqs, t; name=:nothing), ptf)
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
        r_rel_a .~ resolve2(ori(frame_a), frame_b.r_0 - frame_a.r_0);
        zeros(3) .~ collect(frame_b.tau);
        collect(frame_b.f) .~ -resolve2(Rrel, frame_a.f);
        zeros(3) .~ collect(frame_a.tau) + resolve1(Rrel, frame_b.tau) - cross(r_rel_a, frame_a.f);
    ]

    sys = extend(ODESystem(eqs, t; name=:nothing), ptf)
    add_params(sys, pars; name)
end

"""
    UniversalSpherical(; name, n1_a, rRod_ia, sphere_diameter = 0.1, sphere_color, rod_width = 0.1, rod_height = 0.1, rod_color, cylinder_length = 0.1, cylinder_diameter = 0.1, cylinder_color, kinematic_constraint = true)

Universal - spherical joint aggregation (1 constraint, no potential states)

This component consists of a universal joint at `frame_a` and a spherical joint at `frame_b` that are connected together with a rigid rod.

This joint aggregation has no mass and no inertia and introduces the constraint that the distance between the origin of `frame_a` and the origin of `frame_b` is constant (= `length(rRod_ia)`). The universal joint is defined in the following way:
- The rotation axis of revolute joint 1 is along parameter vector `n1_a` which is fixed in `frame_a`.
- The rotation axis of revolute joint 2 is perpendicular to axis 1 and to the line connecting the universal and the spherical joint.

Note, there is a singularity when axis 1 and the connecting rod are parallel to each other. Therefore, if possible `n1_a` should be selected in such a way that it is perpendicular to `rRod_ia` in the initial configuration (i.e., the distance to the singularity is as large as possible).

An additional `frame_ia` is present. It is fixed in the connecting rod at the origin of `frame_a`. The placement of `frame_ia` on the rod is implicitly defined by the universal joint (frame_a and `frame_ia` coincide when the angles of the two revolute joints of the universal joint are zero) and by parameter vector `rRod_ia`, the position vector from the origin of `frame_a` to the origin of `frame_b`, resolved in `frame_ia`.

This joint aggregation can be used in cases where in reality a rod with spherical joints at end are present. Such a system has an additional degree of freedom to rotate the rod along its axis. In practice this rotation is usually of no interest and is mathematically removed by replacing one of the spherical joints by a universal joint. Still, in most cases the [`SphericalSpherical`](@ref) joint aggregation can be used instead of the UniversalSpherical joint since the rod is animated and its mass properties are approximated by a point mass in the middle of the rod. The [`SphericalSpherical`](@ref) joint has the advantage that it does not have a singular configuration.

# Arguments
- `n1_a` Axis 1 of universal joint resolved in frame_a (axis 2 is orthogonal to axis 1 and to rod)
- `rRod_ia` Vector from origin of frame_a to origin of frame_b, resolved in `frame_ia` (if computeRodLength=true, rRod_ia is only an axis vector along the connecting rod)
- `kinematic_constraint = true` Set to false if no constraint shall be defined, due to analytically solving a kinematic loop

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
                    cylinder_length = 0.1,
                    cylinder_diameter = 0.1,
                    cylinder_color = [1, 0.2, 0, 1], 
                    kinematic_constraint = true,
)

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
    end

    vars = @variables begin
        f_rod(t), [description="Constraint force in direction of the rod (positive on frame_a, when directed from frame_a to frame_b)"]
        # rodLength(t), [description="Length of rod (distance between origin of frame_a and origin of frame_b)"]
        # eRod_ia(t)[1:3], [description="Unit vector from origin of frame_a to origin of frame_b, resolved in frame_ia"]
        # e2_ia(t)[1:3], [description="Unit vector in direction of axis 2 of universal joint, resolved in frame_ia (orthogonal to n1_a and eRod_ia; note: frame_ia is parallel to frame_a when the universal joint angles are zero)"]
        # e3_ia(t)[1:3], [description="Unit vector perpendicular to eRod_ia and e2_ia, resolved in frame_ia"]
        f_b_a1(t)[1:3], [description="frame_b.f without f_rod part, resolved in frame_a (needed for analytic loop handling)"]
        eRod_a(t)[1:3], [description="Unit vector in direction of rRod_a, resolved in frame_a (needed for analytic loop handling)"]
        rRod_0(t)[1:3] = rRod_ia, [description="Position vector from frame_a to frame_b resolved in world frame"]
        rRod_a(t)[1:3] = rRod_ia, [description="Position vector from frame_a to frame_b resolved in frame_a"]
        (constraintResidue(t) = 0), [description="Constraint equation of joint in residue form: Either length constraint (= default) or equation to compute rod force (for analytic solution of loops in combination with RevoluteWithLengthConstraint)"]
        f_b_a(t)[1:3], [description="frame_b.f resolved in frame_a"]
        f_ia_a(t)[1:3], [description="frame_ia.f resolved in frame_a"]
        t_ia_a(t)[1:3], [description="frame_ia.t resolved in frame_a"]
        n2_a(t)[1:3], [description="Vector in direction of axis 2 of the universal joint (e2_ia), resolved in frame_a"]
        (length2_n2_a(t) = 1), [description="Square of length of vector n2_a"]
        length_n2_a(t), [description="Length of vector n2_a"]
        e2_a(t)[1:3], [description="Unit vector in direction of axis 2 of the universal joint (e2_ia), resolved in frame_a"]
        e3_a(t)[1:3], [description="Unit vector perpendicular to eRod_ia and e2_a, resolved in frame_a"]
        der_rRod_a_L(t)[1:3], [description="= der(rRod_a)/rodLength"]
        w_rel_ia1(t)[1:3], [description="Angular velocity of frame_ia1 with respect to frame_a"]
        # R_rel_ia1(t), [description="Rotation object from frame_a to frame_ia1 (frame that is fixed in frame_ia such that x-axis is along the rod axis)"]
        # R_rel_ia2(t), [description="Fixed rotation object from frame_ia1 to frame_ia"]
        # R_rel_ia(t), [description="Rotation from frame_a to frame_ia"]
    end

    n1_a, rRod_ia = collect.((n1_a, rRod_ia))

    eRod_ia = normalize(rRod_ia)
    e2_ia = cross(n1_a, eRod_ia)
    e3_ia = cross(eRod_ia, e2_ia)
    rodLength = _norm(rRod_ia)

    eRod_ia,e2_ia,e3_ia,f_b_a1,eRod_a,rRod_0,rRod_a,f_b_a,f_ia_a,t_ia_a,n2_a,e2_a,e3_a,der_rRod_a_L,w_rel_ia1 = collect.((eRod_ia,e2_ia,e3_ia,f_b_a1,eRod_a,rRod_0,rRod_a,f_b_a,f_ia_a,t_ia_a,n2_a,e2_a,e3_a,der_rRod_a_L,w_rel_ia1))


    R_rel_ia1 = RotationMatrix(transpose([eRod_a e2_a e3_a]), w_rel_ia1)
    R_rel_ia2 = RotationMatrix([eRod_ia e2_ia e3_ia], zeros(3))
    # R_rel_ia = ori(R_rel_ia_frame, true)
    # R_rel_ia = ori(R_rel_ia_frame)
    R_rel_ia = absolute_rotation(R_rel_ia1, R_rel_ia2)

    Ra = ori(frame_a)
    
    Ria = absolute_rotation(Ra, R_rel_ia)

    eqs = [

        if kinematic_constraint
            rRod_0 .~ ori(frame_b).R.mat'*(ori(frame_b).R.mat*collect(frame_b.r_0)) - Ra.R.mat'*(Ra.R.mat*collect(frame_a.r_0))
        else
            rRod_0 .~ frame_b.r_0 - frame_a.r_0
        end

        rRod_a .~ resolve2(Ra, rRod_0)

        constraintResidue ~ rRod_0'rRod_0 - rodLength'rodLength
        constraintResidue ~ 0

        eRod_a .~ rRod_a/rodLength
        n2_a .~ cross(n1_a, eRod_a)
        length2_n2_a ~ n2_a'n2_a

        # assert(length2_n2_a > 1e-10, "A MultiBody.Joints.UniversalSpherical joint (consisting of a universal joint and a spherical joint connected together by a rigid rod) is in the singular configuration of the universal joint. This means that axis 1 of the universal joint defined via parameter \"n1_a\" is parallel to vector \"rRod_ia\" that is directed from the origin of frame_a to the origin of frame_b. You may try to use another \"n1_a\" vector. If this fails, use instead MultiBody.Joints.SphericalSpherical, if this is possible, because this joint aggregation does not have a singular configuration.")

        length_n2_a ~ sqrt(length2_n2_a)
        e2_a .~ n2_a/length_n2_a
        e3_a .~ cross(eRod_a, e2_a)

        der_rRod_a_L .~ (resolve2(Ra, D.(rRod_0)) - cross(Ra.w, rRod_a))/rodLength
        w_rel_ia1 .~ [e3_a'cross(n1_a, der_rRod_a_L)/length_n2_a, -(e3_a'der_rRod_a_L), e2_a'der_rRod_a_L]

        # R_rel_ia ~ Rrelia# absolute_rotation(R_rel_ia1, R_rel_ia2)
        # R_rel_ia.w ~ Rrelia.w

        frame_ia.r_0 .~ frame_a.r_0
        ori(frame_ia) ~ Ria # absolute_rotation(frame_a, R_rel_ia)
        # ori(frame_ia).w .~ Ria.w

        f_ia_a .~ resolve1(R_rel_ia, frame_ia.f)
        t_ia_a .~ resolve1(R_rel_ia, frame_ia.tau)

        f_b_a1 .~ -e2_a*((n1_a't_ia_a)/(rodLength*(n1_a'e3_a))) + e3_a*((e2_a't_ia_a)/rodLength)
        f_b_a .~ -f_rod*eRod_a + f_b_a1
        frame_b.f .~ resolve_relative(f_b_a, Ra, ori(frame_b))
        frame_b.tau .~ 0
        0 .~ collect(frame_a.f) + f_b_a + f_ia_a
        0 .~ collect(frame_a.tau) + t_ia_a + cross(rRod_a, f_b_a)
    ]

    sys = ODESystem(eqs, t; name=:nothing, systems)
    add_params(sys, pars; name)
end

