using Test
import ModelingToolkitStandardLibrary.Mechanical.Rotational
@mtkmodel FurutaPendulum begin
    @components begin
        world = W()
        shoulder_joint = Revolute(n = [0, 1, 0], isroot = true, axisflange = true)
        elbow_joint    = Revolute(n = [0, 0, 1], isroot = true, axisflange = true, phi0=0.1)
        upper_arm = BodyShape(; m = 0.1, isroot = false, r = [0, 0, 0.6], radius=0.04)
        lower_arm = BodyShape(; m = 0.1, isroot = false, r = [0, 0.6, 0], radius=0.04)
        tip = Body(; m = 0.3, isroot = false)

        damper1 = Rotational.Damper(d = 0.07)
        damper2 = Rotational.Damper(d = 0.07)
    end
    @equations begin
        connect(world.frame_b, shoulder_joint.frame_a)
        connect(shoulder_joint.frame_b, upper_arm.frame_a)
        connect(upper_arm.frame_b, elbow_joint.frame_a)
        connect(elbow_joint.frame_b, lower_arm.frame_a)
        connect(lower_arm.frame_b, tip.frame_a)

        connect(shoulder_joint.axis, damper1.flange_a)
        connect(shoulder_joint.support, damper1.flange_b)

        connect(elbow_joint.axis, damper2.flange_a)
        connect(elbow_joint.support, damper2.flange_b)

    end
end

@named model = FurutaPendulum()
model = complete(model)
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys, [model.shoulder_joint.phi => 0.0, model.elbow_joint.phi => 0.1], (0, 12))
sol = solve(prob, Rodas4())





get_rot(sol, model.shoulder_joint.frame_b, 0)

# we see that at time $t = 0$, we have no rotation of `frame_b` around the $y$ axis of the world (frames are always resolved in the world frame), but a second into the simulation, we have:
R1 = get_rot(sol, model.shoulder_joint.frame_b, 1)

# Here, the `frame_b` has rotated around the $y$ axis of the world (if you are not familiar with rotation matrices, we can ask for the rotation axis and angle)
using Multibody.Rotations
(rotation_axis(R1), rotation_angle(R1))


# The next body is the upper arm. This body has an extent of `0.6` in the $z$ direction, as measured in its local `frame_a`
get_trans(sol, model.upper_arm.frame_b, 0)

# One second into the simulation, the upper arm has rotated around the $y$ axis of the world
rb1 = get_trans(sol, model.upper_arm.frame_b, 1)
@test rb1 ≈ [-0.2836231361058395, 0.0, 0.528732367711197]

# If we look at the variable `model.upper_arm.r`, we do not see this rotation!
arm_r = sol(1, idxs=collect(model.upper_arm.r))
@test arm_r == [0, 0, 0.6]

# The reason is that this variable is resolved in the local `frame_a` and not in the world frame. To transform this variable to the world frame, we may multiply with the rotation matrix of `frame_a` which is always resolved in the world frame:
@test get_rot(sol, model.upper_arm.frame_a, 1)*arm_r ≈ rb1

# We now get the same result has when we asked for the translation vector of `frame_b` above.

# Slightly more formally, let $R_A^B$ denote the rotation matrix that rotates a vector expressed in a frame $A$ into one that is expressed in frame $B$, i.e., $r_B = R_B^A r_A$. We have then just performed the transformation $r_W = R_W^A r_A$, where $W$ denotes the world frame, and $A$ denotes `body.frame_a`.

# The next joint, the elbow joint, has the rotational axis `n = [0, 0, 1]`. This indicates that the joint rotates around the $z$-axis of its `frame_a`. Since the upper arm was oriented along the $z$ direction, the joint is rotating around the axis that coincides with the upper arm. 

# The lower arm is finally having an extent along the $y$-axis. At the final time when the pendulum motion has been fully damped, we see that the second frame of this body ends up with an $y$-coordinate of `-0.6`:
t1 = get_trans(sol, model.lower_arm.frame_b, 12)
@test t1 ≈ [-0.009040487302666853, -0.59999996727278, 0.599931920189277]


# If we rotate the vector of extent of the lower arm to the world frame, we indeed see that the only coordinate that is nonzero is the $y$ coordinate:
t1rot = get_rot(sol, model.lower_arm.frame_a, 12)*sol(12, idxs=collect(model.lower_arm.r))

@test abs(t1rot[1]) < 1e-2
@test abs(t1rot[3]) < 1e-2


# The reason that the latter vector differs from `get_trans(sol, model.lower_arm.frame_b, 12)` above is that `get_trans(sol, model.lower_arm.frame_b, 12)` has been _translated_ as well. To both translate and rotate `model.lower_arm.r` into the world frame, we must use the full transformation matrix $T_W^A \in SE(3)$:

r_A = sol(12, idxs=collect(model.lower_arm.r))
r_A = [r_A; 1] # Homogeneous coordinates

@test (get_frame(sol, model.lower_arm.frame_a, 12)*r_A)[1:3] ≈ t1

# the vector is now coinciding with `get_trans(sol, model.lower_arm.frame_b, 12)`.

@test get_trans(sol, model.tip.frame_a, 12) ≈ t1
@test get_frame(sol, model.tip.frame_a, 12)[1:3, 4] ≈ t1