import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica as TP

function LineForceBase(; name, length=0, s_small = 1e-10)
    @named frame_a = Frame()
    @named frame_b = Frame()

    @variables length(t) [
        description = "Distance between the origin of frame_a and the origin of frame_b",
    ]
    @variables s(t) [
        description = "(Guarded) distance between the origin of frame_a and the origin of frame_b (>= s_small))",
    ]
    @variables r_rel_0(t)[1:3] [
        description = "Position vector from frame_a to frame_b resolved in world frame",
    ]
    @variables e_rel_0(t)[1:3] [
        description = "Unit vector in direction from frame_a to frame_b, resolved in world frame",
    ]

    eqs = [r_rel_0 .~ frame_b.r_0 .- frame_a.r_0;
           length ~ norm(r_rel_0)
           s ~ max(length, s_small)
           e_rel_0 .~ r_rel_0 ./ s]

    compose(ODESystem(eqs, t; name), frame_a, frame_b)
end

function LineForceWithMass(; name, length=0, m = 1.0, lengthFraction = 0.5)
    m0 = m
    @named lfb = LineForceBase(; length)
    @unpack length, s, r_rel_0, e_rel_0, frame_a, frame_b = lfb
    @named flange_a = TP.Flange()
    @named flange_b = TP.Flange()
    @parameters m=m [description = "mass", bounds = (0, Inf)]
    @variables fa(t)=0 [description = "scalar force from flange_a"]
    @variables fb(t)=0 [description = "scalar force from flange_b"]
    @variables r_CM_0(t)[1:3]=zeros(3) [
        description = "Position vector from world frame to point mass, resolved in world frame",
    ]
    @variables v_CM_0(t)[1:3]=zeros(3) [description = "First derivative of r_CM_0"]
    @variables ag_CM_0(t)[1:3]=zeros(3) [description = "D(v_CM_0) - gravityAcceleration"]
    @parameters lengthFraction=lengthFraction [
        description = "Location of point mass with respect to frame_a as a fraction of the distance from frame_a to frame_b",
        bounds = (0, 1),
    ]

    r_rel_0, e_rel_0, r_CM_0, v_CM_0, ag_CM_0 = collect.((r_rel_0, e_rel_0, r_CM_0, v_CM_0, ag_CM_0))

    eqs = [flange_a.s ~ 0
           flange_b.s ~ length]

    # Determine translational flange forces
    # if cardinality(flange_a) > 0 && cardinality(flange_b) > 0
    # the cardinality of a flange is the number of connections to a flange
    eqs = [eqs
           fa ~ flange_a.f
           fb ~ flange_b.f]
    # elseif cardinality(flange_a) > 0 && cardinality(flange_b) == 0
    #   fa ~ flange_a.f
    #   fb ~ -fa
    # elseif cardinality(flange_a) == 0 && cardinality(flange_b) > 0
    #   fa ~ -fb
    #   fb ~ flange_b.f
    # else
    #   fa ~ 0
    #   fb ~ 0
    # end

    #= Force and torque balance of point mass
       - Kinematics for center of mass CM of point mass including gravity
         r_CM_0 = frame_a.r0 + r_rel_CM_0;
         v_CM_0 = D.(r_CM_0);
         ag_CM_0 = D.(v_CM_0) - world.gravityAcceleration(r_CM_0);
       - Power balance for the connection line
         (f1=force on frame_a side, f2=force on frame_b side, h=lengthFraction)
         0 = f1*va - m*ag_CM*(va+(vb-va)*h) + f2*vb
           = (f1 - m*ag_CM*(1-h))*va + (f2 - m*ag_CM*h)*vb
         since va and vb are completely independent from other
         the parenthesis must vanish:
           f1 := m*ag_CM*(1-h)
           f2 := m*ag_CM*h
       - Force balance on frame_a and frame_b finally results in
           0 = frame_a.f + e_rel_a*fa - f1_a
           0 = frame_b.f + e_rel_b*fb - f2_b
         and therefore
           frame_a.f = -e_rel_a*fa + m*ag_CM_a*(1-h)
           frame_b.f = -e_rel_b*fb + m*ag_CM_b*h
    =#
    if m0 > 0
        eqs = [eqs
               r_CM_0 .~ frame_a.r_0 + r_rel_0 * lengthFraction
               v_CM_0 .~ D.(r_CM_0)
               ag_CM_0 .~ D.(v_CM_0) - gravity_acceleration(r_CM_0)
               frame_a.f .~ resolve2(ori(frame_a),
                                     (m * (1 - lengthFraction)) * ag_CM_0 - e_rel_0 * fa)
               frame_b.f .~ resolve2(ori(frame_b),
                                     (m * lengthFraction) * ag_CM_0 - e_rel_0 * fb)]
    else
        eqs = [eqs
               r_CM_0 .~ zeros(3)
               v_CM_0 .~ zeros(3)
               ag_CM_0 .~ zeros(3)
               frame_a.f .~ -resolve2(ori(frame_a), e_rel_0 * fa)
               frame_b.f .~ -resolve2(ori(frame_b), e_rel_0 * fb)]
    end

    extend(ODESystem(eqs, t; name, systems = [flange_a, flange_b]), lfb)
end

function Spring(c; name, m = 0, lengthFraction = 0.5, s_unstretched = 0)
    @named ptf = PartialTwoFrames()
    @unpack frame_a, frame_b = ptf
    @named lineForce = LineForceWithMass(; length = s_unstretched, m, lengthFraction)
    @parameters c=c [description = "spring constant", bounds = (0, Inf)]
    @parameters s_unstretched=s_unstretched [
        description = "unstretched length of spring",
        bounds = (0, Inf),
    ]

    @variables r_rel_a(t)[1:3] [
        description = "Position vector from origin of frame_a to origin of frame_b, resolved in frame_a",
    ]
    @variables e_a(t)[1:3] [
        description = "Unit vector on the line connecting the origin of frame_a with the origin of frame_b resolved in frame_a (directed from frame_a to frame_b)",
    ]
    @variables f(t) [
        description = "Line force acting on frame_a and on frame_b (positive, if acting on frame_b and directed from frame_a to frame_b)",
    ]
    @variables length(t) [
        description = "Distance between the origin of frame_a and the origin of frame_b",
    ]
    @variables s(t) [
        description = "(Guarded) distance between the origin of frame_a and the origin of frame_b (>= s_small))",
    ]
    @variables r_rel_0(t)[1:3] [
        description = "Position vector from frame_a to frame_b resolved in world frame",
    ]
    @variables e_rel_0(t)[1:3] [
        description = "Unit vector in direction from frame_a to frame_b, resolved in world frame",
    ]
    @named spring2d = TP.Spring(c; s_rel0=s_unstretched)

    eqs = [
        r_rel_a .~ resolve2(ori(frame_a), r_rel_0)
        e_a .~ r_rel_a/s
        f ~ spring2d.f
        length ~ lineForce.length
        s ~ lineForce.s
        r_rel_0 .~ lineForce.r_rel_0
        e_rel_0 .~ lineForce.e_rel_0
    
        connect(lineForce.frame_a, frame_a)
        connect(lineForce.frame_b, frame_b)
        connect(spring2d.flange_b, lineForce.flange_b)
        connect(spring2d.flange_a, lineForce.flange_a)
    ]

    extend(ODESystem(eqs, t; name, systems = [lineForce, spring2d]), ptf)
end
