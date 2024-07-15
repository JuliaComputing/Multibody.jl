# Quadrotor with cable-suspended load


```@example QUAD
using Multibody
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using Plots
using JuliaSimCompiler
using OrdinaryDiffEq
using Test
t = Multibody.t

world = Multibody.world

num_arms = 4 # Number of arms of the rotor craft.
angle_between_arms = 2*pi/num_arms # Angle between the arms, assuming they are evenly spaced.
arm_length = 0.2 # Length of each arm.
arm_outer_diameter = 0.03
arm_inner_diameter = 0.02
arm_density = 800 # Density of the arm [kg/m³].
body_mass = 0.2 # Mass of the body.
load_mass = 1.0 # Mass of the load.
cable_length = 1 # Length of the cable.
cable_mass = 0.1 # Mass of the cable.
cable_diameter = 0.01 # Diameter of the cable.
number_of_links = 3 # Number of links in the cable.
thrust = 13

arms = [
    BodyCylinder(
        r = [arm_length*cos(angle_between_arms*(i-1)), 0, arm_length*sin(angle_between_arms*(i-1))],
        diameter = arm_outer_diameter,
        density = arm_density,
        name=Symbol("arm$i")
    ) for i = 1:num_arms
]

@mtkmodel Thruster begin
    @components begin
        frame_b = Frame()
        thrust3d = WorldForce(resolve_frame = :frame_b)
        thrust = RealInput()
    end
    @equations begin
        thrust3d.force.u[1] ~ 0
        thrust3d.force.u[2] ~ thrust.u
        thrust3d.force.u[3] ~ 0
        connect(frame_b, thrust3d.frame_b)
    end
end

thrusters = [Thruster(name = Symbol("thruster$i")) for i = 1:num_arms]

@named body = Body(m = body_mass)
@named load = Body(m = load_mass)

@named cable = Rope(
    l = 1,
    m = cable_mass,
    n = number_of_links,
    c = 0,
    d = 0,
    air_resistance = 0.1,
    d_joint = 0.1,
    radius = cable_diameter/2,
    color = [0.5, 0.4, 0.4, 1],
    dir = [0.0, -1, 0]
)

@named freemotion = FreeMotion( # This connects the rotorcraft to the world, and is needed to give the rotorcraft a state.
    state = true,
    isroot = true,
    quat = true,
)

connections = [
    connect(world.frame_b, freemotion.frame_a)
    connect(freemotion.frame_b, body.frame_a)
    [connect(body.frame_a, arms[i].frame_a) for i = 1:num_arms]
    connect(body.frame_a, cable.frame_a)
    connect(cable.frame_b, load.frame_a)
    [connect(arms[i].frame_b, thrusters[i].frame_b) for i = 1:num_arms]
    [thrusters[i].thrust.u ~ thrust * t for i = 1:num_arms]
]

@named model = ODESystem(connections, t, systems = [world; body; load; arms; thrusters; cable; freemotion])
model = complete(model)
ssys = structural_simplify(IRSystem(model))


# ModelingToolkit.generate_initializesystem(IRSystem(model); u0map)

prob = ODEProblem(ssys, [
    load.v_0[1] => 2;
    body.v_0[1] => 0;
    collect(cable.joint_2.phi_d) .=> 1;
    ], (0, 13))

sol = solve(prob, Tsit5())
@test SciMLBase.successful_retcode(sol)
plot(sol) |> display

import GLMakie
# first(render(model, sol, 0, show_axis=true)) # Interactive plot
Multibody.render(model, sol, filename = "quad.gif")
# nothing # hide
```


![quadrotor animation](quad.gif)