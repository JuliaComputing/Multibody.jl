using LinearAlgebra


function isroot(sys)
    sys.metadata isa Dict || return false
    get(sys.metadata, :isroot, false)
end

function ori(sys)
    if sys.metadata isa Dict && (O = get(sys.metadata, :orientation, nothing)) !== nothing
        R,w = collect.((O.R, O.w))
        # Since we are using this function instead of sys.ori, we need to handle namespacing properly as well
        ns = nameof(sys)
        R = ModelingToolkit.renamespace.(ns, R) .|> Num
        w = get_w(R)
        # w = ModelingToolkit.renamespace.(ns, w) .|> Num
        RotationMatrix(R, w)
    else
        error("System $(sys.name) does not have an orientation object.")
    end
end

function World(; name)
    # World should have
    # 3+3+9+3 // r_0+f+R.R+τ
    # - (3+3) // (f+t)
    # = 12 equations 
    @named frame_b = Frame()
    @parameters n[1:3]=[0, 1, 0] [description="gravity direction of world"]
    @parameters g=9.81 [description="gravitational acceleration of world"]
    O = ori(frame_b)
    eqs = Equation[collect(frame_b.r_0) .~ 0;
                O ~ nullrotation();
                # vec(D(O).R .~ 0); # QUESTION: not sure if I should have to add this, should only have 12 equations according to modelica paper
                ]
    ODESystem(eqs, t, [], [n; g]; name, systems=[frame_b])
end

const world = World(; name=:world)

"Function to compute the gravity acceleration, resolved in world frame"
gravity_acceleration(r) = world.g*world.n # NOTE: This is hard coded for now to use the the standard, parallel gravity model

function FixedTranslation(; name)
    @named frame_a = Frame(name = name)
    @named frame_b = Frame(name = name)
    @variables r_ab(t)[1:3] [description="position vector from frame_a to frame_b, resolved in frame_a"]
    fa = frame_a.f
    fb = frame_b.f
    taua = frame_a.tau
    taub = frame_b.tau
    eqs = Equation[
        frame_b.r_0 .~ frame_a.r_0 + resolve1(ori(frame_a), r_ab)
        ori(frame_b) ~ ori(frame_a)
        0 .~ fa + fb
        0 .~ taua + taub + cross(r_ab, fb)
    ]
    compose(ODESystem(eqs, t; name), frame_b)
end

function Revolute(; name, ϕ0=0, ω0=0, n=[0, 0, 1])
    @named frame_a = Frame(name = name)
    @named frame_b = Frame(name = name)
    @parameters n[1:3]=n [description="axis of rotation"]
    @variables ϕ(t)=ϕ0 [description="angle of rotation (rad)"]
    @variables ω(t)=ω0 [description="angular velocity (rad/s)"]
    Rrel = planar_rotation(n, ϕ0, ω0)

    
    eqs = Equation[
        frame_a.r_0 .~ frame_b.r_0
        Rrel ~ planar_rotation(n, ϕ, ω)
        ori(frame_b) ~ abs_rotation(ori(frame_a), Rrel)
        D(ϕ) ~ ω

        0 .~ frame_a.f   + resolve1(Rrel, frame_b.f)
        0 .~ frame_a.tau + resolve1(Rrel, frame_b.tau)
        0 .~ n'frame_b.tau # no torque through joint
    ]
    compose(ODESystem(eqs, t; name), frame_b)
end


function Body(; name, m=1, r_cm=[0, 0, 0], I=collect(0.001LinearAlgebra.I(3)), isroot=false)
    @named frame_a = Frame()
    f = frame_a.f |> collect
    tau = frame_a.tau |> collect
    @variables r_0(t)[1:3]=0 [description="Position vector from origin of world frame to origin of frame_a"]
    @variables v_0(t)[1:3]=0 [description="Absolute velocity of frame_a, resolved in world frame (= D(r_0))"]
    @variables a_0(t)[1:3]=0 [description="Absolute acceleration of frame_a resolved in world frame (= D(v_0))"]
    @variables g_0(t)[1:3]=0 [description="gravity acceleration"]
    @variables w_a(t)[1:3]=0 [description="Absolute angular velocity of frame_a resolved in frame_a"]
    @variables z_a(t)[1:3]=0 [description="Absolute angular acceleration of frame_a resolved in frame_a"]
    # 6*3 potential variables + Frame: 2*3 flow + 3 potential + 3 residual = 24 equations + 2*3 flow
    @parameters m=m [description="mass"]
    @parameters r_cm[1:3]=r_cm [description="center of mass"]
    @parameters I[1:3, 1:3]=I [description="inertia tensor"]


    @parameters I_11 = I[1,1] [description="Element (1,1) of inertia tensor"]
    @parameters I_22 = I[2,2] [description="Element (2,2) of inertia tensor"]
    @parameters I_33 = I[3,3] [description="Element (3,3) of inertia tensor"]
    @parameters I_21 = I[2,1] [description="Element (2,1) of inertia tensor"]
    @parameters I_31 = I[3,1] [description="Element (3,1) of inertia tensor"]
    @parameters I_32 = I[3,2] [description="Element (3,2) of inertia tensor"]

    I = [I_11 I_21 I_31; I_21 I_22 I_32; I_31 I_32 I_33]

    r_0,v_0,a_0,g_0,w_a,z_a,r_cm = collect.((r_0,v_0,a_0,g_0,w_a,z_a,r_cm))

    Ra = ori(frame_a)
    DRa = D(Ra)

    eqs = if isroot # isRoot
        Equation[
            0 .~ orientation_constraint(Ra); # TODO: Replace with non-unit quaternion
            # collect(w_a .~ DRa.w);
            # q̇ .~ D.(q)
            # w_a .~ 0;# angular_velocity2(q, q̇)
        ]
    else
        # This equation is defined here and not in the Rotation component since the branch above might use another equation
        [
            collect(w_a .~ Ra.w);
            # collect(w_a .~ DRa.w); # angular_velocity2(R, D.(R)): skew(R.w) = R.T*der(transpose(R.T))
            # vec(DRa.R .~ 0)
        ]
    end

    eqs = [
        eqs;
        collect(r_0 .~ frame_a.r_0)
        collect(g_0 .~ [0, -9.82, 0])# NOTE: should be gravity_acceleration(frame_a.r_0 .+ resolve1(Ra, r_cm)))
        collect(v_0 .~ D.(r_0))
        collect(a_0 .~ D.(v_0))
        collect(z_a .~ D.(w_a))

        collect(f .~ m*(resolve2(Ra, a_0 - g_0) + cross(z_a, r_cm) + cross(w_a, cross(w_a, r_cm))))
        collect(tau .~ I*z_a + cross(w_a, I*w_a) + cross(r_cm, f))
    ]

    ODESystem(eqs, t; name, metadata = Dict(:isroot => isroot), systems=[frame_a])
end