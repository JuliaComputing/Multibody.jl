using LinearAlgebra

function World(; name)
    @named frame_b = Frame()
    @parameters n[1:3]=[0, 1, 0] [description="gravity direction of world"]
    @parameters g=9.81 [description="gravitational acceleration of world"]
    eqs = Equation[collect(frame_b.r_0) .~ 0
                orientation_equal(frame_b.R, nullrotation())]
    ODESystem(eqs, t, [], [n; g]; name, systems=[frame_b])
end

const world = World(; name=:world)

"Function to compute the gravity acceleration, resolved in world frame"
gravity_acceleration(r) = world.g*world.n # NOTE: This is hard coded to use the the standard, parallel gravity model




function FixedTranslation(; name)
    @named frame_a = Frame(name = name)
    @named frame_b = Frame(name = name)
    @variables r_ab(t)[1:3] [description="position vector from frame_a to frame_b, resolved in frame_a"]
    fa = frame_a.f
    fb = frame_b.f
    taua = frame_a.tau
    taub = frame_b.tau
    eqs = Equation[
        frame_b.r_0 .~ frame_a.r_0 + resolve1(frame_a.R, r_ab)
        orientation_equal(frame_b.R, frame_a.R)
        0 .~ fa + fb
        0 .~ taua + taub + cross(r_ab, fb)
    ]
    compose(ODESystem(eqs, t; name), frame_b)
end

function Revolute(; name, ϕ0=0, ω0=0)
    @named frame_a = Frame(name = name)
    @named frame_b = Frame(name = name)
    @parameters n[1:3]=[0, 0, 1] [description="axis of rotation"]
    @variables ϕ0(t)=ϕ0 [description="angle of rotation (rad)"]
    @variables ω0(t)=ω0 [description="angular velocity (rad/s)"]
    @variables Rrel(t)[1:3, 1:3]=planar_rotation(n, ϕ0) [description="relative rotation matrix"]

    
    eqs = Equation[
        frame_a.r_0 .~ frame_b.r_0
        orientation_equal(Rrel, planar_rotation(n, ϕ))
        orientation_equal(frame_b.R, abs_rotation(frame_a.R, Rrel))
        D(ϕ) ~ ω

        0 .~ frame_a.f   + resolve1(Rrel, frame_b.f)
        0 .~ frame_a.tau + resolve1(Rrel, frame_b.tau)
        0 .~ n'frame_b.tau # no torque through joint
    ]
    compose(ODESystem(eqs, t; name), frame_b)
end


function Body(; name, m=1, r_cm=[0, 0, 0], I=collect(0.001LinearAlgebra.I(3)), isroot=false)
    @named frame_a = Frame()
    R = frame_a.R.R |> collect
    w = frame_a.R.w |> collect
    f = frame_a.f
    tau = frame_a.tau
    r = frame_a.r_0
    @variables v(t)[1:3]=0 [description="linear velocity"]
    @variables a(t)[1:3]=0 [description="linear acceleration"]
    @variables w(t)[1:3]=0 
    @variables g(t)[1:3]=0 [description="gravity acceleration"]
    @variables q(t)[1:4]=[0,0,0,1] [description="quaternion orientation"]
    @variables q̇(t)[1:4]=[0,0,0,0] [description="quaternion time derivative"]
    @parameters m=m [description="mass"]
    @parameters r_cm[1:3]=r_cm [description="center of mass"]
    @parameters I[1:3, 1:3]=I [description="inertia tensor"]

    v,a,w,g,q,q̇,r_cm,I,tau = collect.((v,a,w,g,q,q̇,r_cm,I,tau))

    eqs = if isroot # isRoot
        Equation[
            0 .~ orientation_constraint(q) # TODO: Replace with non-unit quaternion
            q̇ .~ D.(q)
            w .~ angular_velocity2(q, q̇)
        ]
    else
        # This equation is defined here and not in the Rotation component since the branch above might use another equaiton
        w .~ skewcoords(R*D.(R')) # angular_velocity2(R, D.(R)): skew(R.w) = R.T*der(transpose(R.T))
            # Quaternion not used in this branch
    end

    eqs = [
        eqs
        collect(D.(r) .~ v)
        collect(g .~ gravity_acceleration(r + resolve1(R, r_cm)))
        collect(a .~ resolve2(R, D.(v) - g))
        collect(f/m .~ a + cross(D.(w), r_cm) + cross(w, cross(w, r_cm)))
        collect(tau .~ I * D.(w) + cross(w, I * w) + cross(r_cm, f))
    ]

    ODESystem(eqs, t; name, systems=[frame_a])
end