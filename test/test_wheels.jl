using Test
# ==============================================================================
## Rolling wheel ===============================================================
# ==============================================================================
using LinearAlgebra
# The wheel does not need the world
# @testset "Rolling wheel" begin
@mtkmodel WheelInWorld begin
    @components begin
        # world = World(n=[0,0,-1])
        world = World()
        wheel = RollingWheel(radius = 0.3, m = 2, I_axis = 0.06,
                            I_long = 0.12,
                            x0 = 0.2,
                            z0 = 0.2,
                            der_angles = [0, 5, 1])
    end
end

@named worldwheel = WheelInWorld()
worldwheel = complete(worldwheel)

# pars = collect(worldwheel.world.n) .=> [0,0,-1];
defs = Dict([
    # collect(worldwheel.world.n) .=> [0,0,-1];
    worldwheel.wheel.body.r_0[1] => 0.2;
    worldwheel.wheel.body.r_0[2] => 0.3;
    worldwheel.wheel.body.r_0[3] => 0.2;
    # collect(D.(cwheel.wheel.angles)) .=> [0, 5, 1]
])

ssys = structural_simplify(IRSystem(worldwheel))
prob = ODEProblem(ssys, defs, (0, 4))
length(filter(x->occursin("world₊n", string(x)), parameters(worldwheel))) == 3
# @test prob[collect(worldwheel.world.n)] == [0,0,-1]
@test prob[collect(worldwheel.wheel.wheeljoint.der_angles)] == prob[collect(worldwheel.wheel.wheeljoint.der_angles)]

sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8)
@test SciMLBase.successful_retcode(sol)
@test sol(4, idxs=[worldwheel.wheel.x; worldwheel.wheel.z]) ≈ [0.162547, -2.23778] atol=1e-3

# plot(sol, idxs=[
#     worldwheel.wheel.x;
#     worldwheel.wheel.z;
#     # worldwheel.wheel.body.r_0[1];
#     # worldwheel.wheel.body.r_0[3];
# ])


# ==============================================================================
## Rolling wheel on interesting surface ===============================================================
# ==============================================================================
using LinearAlgebra
# The wheel does not need the world
# @testset "Rolling wheel" begin
@mtkmodel WheelInWorld begin
    @structural_parameters begin
        surface
    end
    @components begin
        # world = World(n=[0,0,-1])
        world = World()
        wheel = RollingWheel(; radius = 0.3, m = 2, I_axis = 0.06,
                            I_long = 0.12,
                            x0 = 0.2,
                            z0 = 0.2,
                            surface,
                            der_angles = [0, 0, 0])
    end
end

@named worldwheel = WheelInWorld(surface = (x,z)->0)
worldwheel = complete(worldwheel)

# pars = collect(worldwheel.world.n) .=> [0,0,-1];
defs = Dict([
    # collect(worldwheel.world.n) .=> [0,0,-1];
    worldwheel.wheel.body.r_0[1] => 0.2;
    worldwheel.wheel.body.r_0[2] => 0.3;
    worldwheel.wheel.body.r_0[3] => 0.2;
    collect(worldwheel.wheel.wheeljoint.der_angles) .=> [0, 5, 1];
    # collect(D.(cwheel.wheel.angles)) .=> [0, 5, 1]
])

ssys = structural_simplify(IRSystem(worldwheel))
prob = ODEProblem(ssys, defs, (0, 4))
sol = solve(prob, FBDF(autodiff=false), abstol=1e-8, reltol=1e-8)
@test SciMLBase.successful_retcode(sol)
# first(Multibody.render(worldwheel, sol, 0, show_axis=true))
@test sol(4, idxs=[worldwheel.wheel.x; worldwheel.wheel.z]) ≈ [0.162547, -2.23778] atol=1e-3

@test all(norm.(sol[collect(worldwheel.wheel.wheeljoint.e_lat_0)]) .≈ 1)

@named worldwheel = WheelInWorld(surface = (x,z)->x)
worldwheel = complete(worldwheel)

defs = Dict([
    # collect(worldwheel.world.n) .=> [0,0,-1];
    worldwheel.wheel.body.r_0[1] => 0.0;
    worldwheel.wheel.body.r_0[2] => 0.3/sqrt(2);
    worldwheel.wheel.body.r_0[3] => 0.0;
    collect(worldwheel.wheel.wheeljoint.der_angles) .=> [0, 0, 0];
])

ssys = structural_simplify(IRSystem(worldwheel))
prob = ODEProblem(ssys, defs, (0, 4))
sol = solve(prob, FBDF(autodiff=false), abstol=1e-8, reltol=1e-8)
# plot(sol)
tv = 0:0.5:4
@test sol(tv, idxs=worldwheel.wheel.body.r_0[1]) ≈ sol(tv, idxs=worldwheel.wheel.body.r_0[2]) .- 0.3*sqrt(2) rtol=1e-6 # The sqrt(2) is to account for the shifted contact point at a 45 degree plane

dd = diff(sol(tv, idxs=worldwheel.wheel.wheeljoint.der_angles[2]).u) # angular acceleration
@test norm(dd .- dd[1]) < 1e-10 # constant acceleration
@test abs(dd[1]) < 9.81
@test abs(dd[1]) > 5
@test all(norm.(sol[collect(worldwheel.wheel.wheeljoint.e_lat_0)]) .≈ 1)
@test all(norm.(sol[collect(worldwheel.wheel.wheeljoint.e_long_0)]) .≈ 1)




## Rolling wheel with axis ====================================================
import ModelingToolkitStandardLibrary.Blocks
@mtkmodel WheelWithAxis begin
    @components begin
        world = World()
        prismatic = Prismatic(n = [0,1,0])
        world_axis = Revolute(n = [0,1,0], iscut=false, state_priority=100, w0=10)
        # world_axis = RevolutePlanarLoopConstraint(n = [0,1,0])
        spin_axis = Revolute(n = [0,0,1], iscut=false, state_priority=100)
        # spin_axis = RevolutePlanarLoopConstraint(n = [0,0,1])
        bar = FixedTranslation(r = [0, 0.3, 1])
        wheel = RollingWheel(radius = 0.3, m = 2, I_axis = 0.06,
                            I_long = 0.12,
                            x0 = 1,
                            z0 = 1,
                            iscut=true,
                            der_angles = [0, 5, 0])
    end
    @equations begin
        connect(world.frame_b, prismatic.frame_a)
        connect(prismatic.frame_b, world_axis.frame_a)
        connect(world_axis.frame_b, bar.frame_a)
        connect(bar.frame_b, spin_axis.frame_a)
        connect(spin_axis.frame_b, wheel.frame_a)
    end

end
@named model = WheelWithAxis()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, [], (0, 4))
@test_skip begin # Singular linear system
    sol = solve(prob, Rodas4(autodiff=false), abstol=1e-8, reltol=1e-8)
    @test_broken !all(iszero, sol.u)
end
# first(Multibody.render(model, sol, 0, show_axis=true))

# ==============================================================================
## RollingWheelSet
# ==============================================================================
@mtkmodel DrivingWheelSet begin
    @components begin
        sine1 = Blocks.Sine(frequency=1, amplitude=2)
        sine2 = Blocks.Sine(frequency=1, amplitude=2, phase=pi/2)
        torque1 = Rotational.Torque()
        torque2 = Rotational.Torque()
        wheels = RollingWheelSet(radius=0.1, m_wheel=0.5, I_axis=0.01, I_long=0.02, track=0.5, state_priority=100)
        bar = FixedTranslation(r = [0.2, 0, 0])
        body = Body(m=0.01, state_priority=1)
        world = World()
    end
    @equations begin
        connect(sine1.output, torque1.tau)
        connect(sine2.output, torque2.tau)
        connect(torque1.flange, wheels.axis1)
        connect(torque2.flange, wheels.axis2)
        connect(wheels.frame_middle, bar.frame_a)
        connect(bar.frame_b, body.frame_a)
    end
end

@named model = DrivingWheelSet()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
# display(unknowns(ssys))
prob = ODEProblem(ssys, [
    model.wheels.wheelSetJoint.prismatic1.s => 0.1
    model.wheels.wheelSetJoint.prismatic2.s => 0.1
], (0, 3))
sol = solve(prob, Tsit5())
@test SciMLBase.successful_retcode(sol)

@test sol(0:0.1:1, idxs=[model.wheels.x, model.wheels.z]) ≈ [
    0.1  0.0611068  -0.0658862  -0.270161  -0.519216  -0.741042  -0.834663   -0.820189   -0.821855  -0.862654  -0.888457
    0.1  0.101768    0.122277    0.169771   0.193664   0.120891  -0.0183936  -0.0968069  -0.100183  -0.102072  -0.10982
] atol=1e-3

# plot(sol)
# plot(sol, idxs=[model.wheels.x, model.wheels.z])
# first(Multibody.render(model, sol, 0, show_axis=true))
