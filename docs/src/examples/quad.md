# Quadrotor with cable-suspended load


```@example QUAD
using Multibody
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using LinearAlgebra
using Plots
using JuliaSimCompiler
using OrdinaryDiffEq
using Test
t = Multibody.t
D = Differential(t)

world = Multibody.world

num_arms = 4 # Number of arms of the rotor craft.
angle_between_arms = 2*pi/num_arms # Angle between the arms, assuming they are evenly spaced.
arm_length = 0.2 # Length of each arm.
arm_outer_diameter = 0.03
arm_inner_diameter = 0.02
arm_density = 800 # Density of the arm [kg/m³].
body_mass = 0.2 # Mass of the body.
load_mass = 0.1 # Mass of the load.
cable_length = 1 # Length of the cable.
cable_mass = 0.1 # Mass of the cable.
cable_diameter = 0.01 # Diameter of the cable.
number_of_links = 5 # Number of links in the cable.
kalt = 1
Tialt = 3
Tdalt = 3

kroll = 0.02
Tiroll = 100
Tdroll = 1

kpitch = 0.02
Tipitch = 100
Tdpitch = 1

@mtkmodel Thruster begin
    @components begin
        frame_b = Frame()
        thrust3d = WorldForce(resolve_frame = :frame_b, scale=0.1, radius=0.02)
        thrust = RealInput()
    end
    @variables begin
        u(t), [state_priority=1000]
    end
    @equations begin
        thrust3d.force.u[1] ~ 0
        thrust3d.force.u[2] ~ thrust.u
        thrust3d.force.u[3] ~ 0
        thrust.u ~ u
        connect(frame_b, thrust3d.frame_b)
    end
end

rou(x) = round(x, digits=3)

function RotorCraft(; cl = true, addload=true)
    arms = [
        BodyCylinder(
            r = rou.([arm_length*cos(angle_between_arms*(i-1)), 0, arm_length*sin(angle_between_arms*(i-1))]),
            diameter = arm_outer_diameter,
            inner_diameter = arm_inner_diameter,
            density = arm_density,
            name=Symbol("arm$i"),
        ) for i = 1:num_arms
    ]

    thrusters = [Thruster(name = Symbol("thruster$i")) for i = 1:num_arms]

    # @named feedback_gain = MatrixGain(K = -[kp*I(4) kd*I(4) ki*I(4)])

    @parameters Galt[1:4] = ones(4)
    @parameters Groll[1:4] = [1,0,-1,0]
    @parameters Gpitch[1:4] = [0,1,0,-1]

    @named Calt = PID(; k=kalt, Ti=Tialt, Td=Tdalt)
    @named Croll = PID(; k=kroll, Ti=Tiroll, Td=Tdroll)
    @named Cpitch = PID(; k=kpitch, Ti=Tipitch, Td=Tdpitch)
    # @named Calt = PI(; k=kalt, T=Tialt)
    # @named Croll = PI(; k=kroll, T=Tiroll)
    # @named Cpitch = PI(; k=kpitch, T=Tipitch)

    @named body = Body(m = body_mass, state_priority = 0, I_11=0.01, I_22=0.01, I_33=0.01, air_resistance=1)
    @named load = Body(m = load_mass, air_resistance=1)
    @named freemotion = FreeMotion(state=true, isroot=true, quat=false, state_priority=1000, neg_w=false)

    @named cable = Rope(
        l = cable_length,
        m = cable_mass,
        n = number_of_links,
        c = 0,
        d = 0,
        air_resistance = 0.5,
        d_joint = 0.1,
        radius = cable_diameter/2,
        color = [0.5, 0.4, 0.4, 1],
        dir = [0.0, -1, 0]
    )

    connections = [
        connect(world.frame_b, freemotion.frame_a)
        connect(freemotion.frame_b, body.frame_a)
        [connect(body.frame_a, arms[i].frame_a) for i = 1:num_arms]
        [connect(arms[i].frame_b, thrusters[i].frame_b) for i = 1:num_arms]
    ]
    systems = [world; arms; body; thrusters; freemotion]
    if addload
        push!(systems, load)
        push!(systems, cable)
        
        push!(connections, connect(body.frame_a, cable.frame_a))
        push!(connections, connect(cable.frame_b, load.frame_a))
    end
    if cl

        uc = Galt*Calt.ctr_output.u + Groll*Croll.ctr_output.u + Gpitch*Cpitch.ctr_output.u
        uc = collect(uc)
        append!(connections, [thrusters[i].u ~ uc[i] for i = 1:num_arms])

        append!(connections, [
            Calt.err_input.u ~ -body.r_0[2]
            Croll.err_input.u ~ freemotion.phi[3]
            Cpitch.err_input.u ~ -freemotion.phi[1]
        ])
        append!(systems, [Calt; Croll; Cpitch])

        #=
        # append!(connections, [thrusters[i].thrust.u ~ feedback_gain.output.u[i] for i = 1:num_arms])
        # append!(connections, [feedback_gain.input.u[i] ~ arms[i].frame_b.r_0[2] for i = 1:num_arms ]) # Connect positions to controller
        # append!(connections, [feedback_gain.input.u[i+num_arms] ~ D(arms[i].frame_b.r_0[2]) for i = 1:num_arms]) # Connect velocities to controller
        # append!(connections, [feedback_gain.input.u[i+2num_arms] ~ Ie[i] for i = 1:num_arms]) #
        # append!(connections, [feedback_gain.input.u[i] ~ freemotion.phi[[1,3][i]] for i = 1:2 ]) # Connect positions to controller
        # append!(connections, [feedback_gain.input.u[i+2] ~ freemotion.phid[[1,3][i]] for i = 1:2]) # Connect velocities to controller
        # push!(systems, feedback_gain)
        =#
    end
    @named model = ODESystem(connections, t; systems)
    complete(model)
end
model = RotorCraft(cl=true, addload=false)
ssys = structural_simplify(IRSystem(model))


op = [
    # model.load.v_0[1] => 0.0;
    model.body.v_0[1] => 0;
    # collect(model.freemotion.phi) .=> 0.1;
    # collect(model.cable.joint_2.phi) .=> 0.03;
    model.world.g => 2;
    model.body.frame_a.render => true
    model.body.frame_a.radius => 0.01
    model.body.frame_a.length => 0.1
    # collect(model.body.a_0) .=> 0;
    ]
# ModelingToolkit.generate_initializesystem(ssys; u0map=op)


# init_sys_ir = generate_initializesystem(ssys, u0map = op)
# isys_ir = structural_simplify(init_sys_ir, fully_determined = false)
# init_prob_ir = NonlinearProblem(isys_ir, [model.Calt.int.y => 1], [t => 0.0])
# sol_ir = solve(init_prob_ir)


prob = ODEProblem(ssys, op, (0, 30))
sol = solve(prob, FBDF(autodiff=false), reltol=1e-8, abstol=1e-8)
# sol = solve(prob, Tsit5(), abstol=1e-5, reltol=1e-5)
@test SciMLBase.successful_retcode(sol)
# plot(sol) |> display
plot(sol, idxs=[model.freemotion.phi;], layout=1) |> display

plot(sol, idxs=[model.arm1.frame_b.r_0[2], model.arm2.frame_b.r_0[2], model.arm3.frame_b.r_0[2], model.arm4.frame_b.r_0[2]], layout=4, framestyle=:zerolines) |> display
# plot(sol, idxs=[model.arm1.frame_b.f[2], model.arm2.frame_b.f[2], model.arm3.frame_b.f[2], model.arm4.frame_b.f[2]], layout=1) |> display
# plot(sol, idxs=[model.Calt.ctr_output.u, model.Croll.ctr_output.u, model.Cpitch.ctr_output.u], layout=3) |> display
# plot(sol, idxs=[model.thruster1.u, model.thruster2.u, model.thruster3.u, model.thruster4.u], layout=1) |> display
# plot(sol, idxs=[model.Calt.ctr_output.u], layout=1) |> display

# import GLMakie
# first(render(model, sol, 0, show_axis=true)) # Interactive plot
# Multibody.render(model, sol, filename = "quad.gif")
# nothing # hide
```


![quadrotor animation](quad.gif)

```@example QUAD
outputs = []# [model.freemotion.phi; model.freemotion.phid]
inputs = [model.thruster1.u; model.thruster2.u; model.thruster3.u; model.thruster4.u]
irsys = IRSystem(RotorCraft(cl=false))
ssys, diff_idxs, alge_idxs, input_idxs = JuliaSimCompiler.io_preprocessing(irsys, inputs, [])

ModelingToolkit.linearization_function(irsys, inputs, outputs)

```