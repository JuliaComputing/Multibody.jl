using ModelingToolkit
using Multibody
using Test
using JuliaSimCompiler
t = Multibody.t
D = Differential(t)
using OrdinaryDiffEq
using OrdinaryDiffEq.SciMLBase: successful_retcode

using Multibody: AccSensor,AxisType2,AxisType1,Controller,GearType2,BearingFriction,GearType1,Motor,MechanicalStructure,RobotAxis,Robot6DOF,PathPlanning1,PathPlanning6

doplot() = false

world = Multibody.world
@named structure = MechanicalStructure()



@named motor = Motor()
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks

myfixed(args...; kwargs...) = ModelingToolkitStandardLibrary.Mechanical.Rotational.Fixed(args...; kwargs...)

myinertia(args...; kwargs...) = ModelingToolkitStandardLibrary.Mechanical.Rotational.Inertia(args...; kwargs...)

mytorque(args...; kwargs...) = ModelingToolkitStandardLibrary.Mechanical.Rotational.ConstantTorque(args...; kwargs...)

##

@mtkmodel MotorTest begin
    @components begin
        motor = Motor()
        # fixed = myfixed() # bug in @mtkmodel
        inertia = myinertia(J=1, phi=0, w=0)
        constant = Constant(k=1)
    end
    @equations begin
        # connect(motor.flange_motor, fixed.flange)
        connect(motor.flange_motor, inertia.flange_a)
        constant.output.u ~ motor.axisControlBus.current_ref
    end
end

@named motorTest = MotorTest()
m = structural_simplify(IRSystem(motorTest))
# @test length(unknowns(m)) == 3
    # D(motorTest.motor.gear.bearingFriction.w) => 0
cm = complete(motorTest)

prob = ODEProblem(m, [
    D(D(cm.motor.Jmotor.phi)) => 0,
], (0.0, 5.0))
sol = solve(prob, Rodas4())
@test successful_retcode(sol)
doplot() && plot(sol, idxs=cm.motor.phi.phi.u)

##

@mtkmodel GearTest begin
    @components begin
        motor = Motor()
        inertia = myinertia(J=1, phi=0, w=0)
        gear = GearType2()
        # fixed = myfixed() # bug in @mtkmodel
        constant = Constant(k=1)
    end
    @equations begin
        constant.output.u ~ motor.axisControlBus.current_ref
        connect(motor.flange_motor, gear.flange_a)
        connect(gear.flange_b, inertia.flange_a)
        # connect(gear2.flange_a, fixed.flange)
    end
end

@named gearTest = GearTest()
m = structural_simplify(IRSystem(gearTest))
cm = complete(gearTest)

prob = ODEProblem(m, [
    cm.gear.gear.phi_b => 0,
    D(cm.gear.gear.phi_b) => 0,
], (0.0, 5.0))
sol = solve(prob, Rodas4())
@test successful_retcode(sol)
doplot() && plot(sol, idxs=cm.motor.phi.phi.u)


@mtkmodel GearTest begin
    @components begin
        gear2 = GearType1()
        fixed = myfixed() # bug in @mtkmodel
        constant = Constant(k=1)
    end
    @equations begin
        connect(gear2.flange_a, fixed.flange)
    end
end

@named gearTest = GearTest()
m = structural_simplify(IRSystem(gearTest))

##

@mtkmodel ControllerTest begin
    @components begin
        controller = Controller()
        constant1 = Constant(k=1)
        constant2 = Constant(k=1)
        constant3 = Constant(k=1)
        constant4 = Constant(k=1)
    end
    @equations begin
        constant1.output.u ~ controller.axisControlBus.motorAngle
        constant2.output.u ~ controller.axisControlBus.speed_ref
        constant3.output.u ~ controller.axisControlBus.angle_ref
        constant4.output.u ~ controller.axisControlBus.motorSpeed
    end
end

@named controllerTest = ControllerTest()
m = structural_simplify(controllerTest)
m = structural_simplify(IRSystem(controllerTest))


## Test Axis


@mtkmodel AxisTest2 begin
    @components begin
        axis2 = AxisType2()
        # fixed = myfixed() # bug in @mtkmodel
        # fixed = mytorque(tau_constant=1, use_support=true) # bug in @mtkmodel
        inertia = myinertia(J=1, phi=0, w=0)
        # constant1 = Constant(k=1)
        constant2 = Constant(k=0)
        constant3 = Constant(k=2)
        constant4 = Constant(k=1)
        constant5 = Constant(k=0)
    end
    @equations begin
        # connect(axis2.flange, fixed.flange)
        connect(axis2.flange, inertia.flange_a)
        # constant1.output.u ~ motor.axisControlBus.current_ref # This is connected to the controller output
        constant2.output.u ~ axis2.axisControlBus.speed_ref
        constant3.output.u ~ axis2.axisControlBus.angle_ref
        constant4.output.u ~ axis2.axisControlBus.motion_ref
        constant5.output.u ~ axis2.axisControlBus.acceleration_ref
        # axis2.motor.emf.support.phi ~ 0
    end
end

@named axisTest = AxisTest2()
# m = structural_simplify(axisTest)
m = structural_simplify(IRSystem(axisTest)) # Yingbo: solution unstable with IRSystem simplification

cm = complete(axisTest)
tspan = (0.0, 5.0)
prob = ODEProblem(m, [
    # ModelingToolkit.missing_variable_defaults(m);
    # D(cm.axis2.gear.bearingFriction.w) => 0
    cm.axis2.motor.flange_motor.phi => deg2rad(20) *  0,
    D(cm.axis2.motor.flange_motor.phi) => 0,
    cm.axis2.motor.Jmotor.phi => deg2rad(20) *  0,
    cm.axis2.gear.gear.phi_b => 0,
    D(cm.axis2.gear.gear.phi_b) => 0,
], tspan)
sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)

@test sol(0.0, idxs=cm.axis2.motor.emf.phi) == 0
# @test sol(tspan[2], idxs=cm.axis2.motor.emf.phi) == 0

doplot() && plot(sol, layout=length(unknowns(m)))
doplot() && plot(sol, idxs=[
    cm.axis2.gear.gear.phi_a
    cm.axis2.gear.gear.phi_b
    cm.axis2.gear.gear.flange_b.phi
    # cm.axis2.gear.bearingFriction.flange_a.phi
    cm.axis2.gear.flange_b.phi
    cm.axis2.gear.gear.phi_support
    cm.axis2.angleSensor.phi.u
    cm.axis2.motor.phi.phi.u
], layout=8, size=(800, 800))
u = cm.axis2.controller.PI.ctr_output.u
@test abs(sol(prob.tspan[2], idxs=u)) < 1e-6 # test control output is zero at the end of simulation

##

@named controller = Controller()
@named axis2 = AxisType2()
@named gear2 = GearType2()
# @named axis1 = AxisType1()
@named gear1 = GearType1()

@named pp = PathPlanning1(;)
@named pp6 = PathPlanning6(;)


# ==============================================================================
## Test OneAxis
# ==============================================================================

@testset "one axis" begin
    @info "Testing one axis"
    @named oneaxis = RobotAxis()
    oneaxis = complete(oneaxis)
    op = Dict([
        oneaxis.axis.flange.phi => 0
        D(oneaxis.axis.flange.phi) => 0
        D(D(oneaxis.axis.flange.phi)) => 0
        D(D(oneaxis.load.phi)) => 0
        D(oneaxis.axis.gear.gear.phi_b) => 0
        oneaxis.axis.controller.PI.T => 0.01
        oneaxis.axis.controller.PI.gainPI.k => 1
        oneaxis.axis.controller.P.k => 10
        oneaxis.load.J => 1.3*15
        oneaxis.load.phi => 0
    ])
    # matrices_S, simplified_sys = Blocks.get_sensitivity(oneaxis, :axis₊controller_e; op)


    # using ControlSystemsBase 
    # S = ss(matrices_S...) |> minreal
    # @test isstable(S)
    # bodeplot(S)


    ssys = structural_simplify(IRSystem(oneaxis)) # Yingbo: IRSystem does not handle the DataInterpolations.CubicSpline
    # ssys = structural_simplify(oneaxis)
    # cm = oneaxis
    # prob = ODEProblem(ssys, [
    #     cm.axis.flange.phi => 0
    #     D(cm.axis.flange.phi) => 0
    # ], (0.0, 5.0))

    zdd = ModelingToolkit.missing_variable_defaults(oneaxis); op = merge(Dict(zdd), op)

    prob = ODEProblem(ssys, collect(op), (0.0, 3),)
    sol = solve(prob, Rodas4());
    if doplot()
        plot(sol, layout=length(unknowns(ssys)), size=(1900, 1200))
        plot!(sol, idxs=oneaxis.pathPlanning.controlBus.axisControlBus1.angle_ref)
        display(current())
    end
    @test SciMLBase.successful_retcode(sol)
    # @test sol(10, idxs=oneaxis.axis.controller.PI.err_input.u) ≈ 0 atol=1e-8

    tv = 0:0.1:10.0
    control_error = sol(tv, idxs=oneaxis.pathPlanning.controlBus.axisControlBus1.angle_ref-oneaxis.load.phi)

    @test sol(tv[1], idxs=oneaxis.pathPlanning.controlBus.axisControlBus1.angle_ref) ≈ deg2rad(0) atol=1e-8
    @test sol(tv[end], idxs=oneaxis.pathPlanning.controlBus.axisControlBus1.angle_ref) ≈ deg2rad(120)
    @test maximum(abs, control_error) < 1e-3
end


##
@testset "full robot" begin
    @info "Testing full robot"

    @named robot = Robot6DOF()
    robot = complete(robot)

    @time "full robot" begin 
        @time "structural_simplify" ssys = structural_simplify(IRSystem(robot))
        @time "ODEProblem creation" prob = ODEProblem(ssys, [
            robot.mechanics.r1.phi => deg2rad(-60)
            robot.mechanics.r2.phi => deg2rad(20)
            robot.mechanics.r3.phi => deg2rad(90)
            robot.mechanics.r4.phi => deg2rad(0)
            robot.mechanics.r5.phi => deg2rad(-110)
            robot.mechanics.r6.phi => deg2rad(0)
        
            robot.axis1.motor.Jmotor.phi => deg2rad(-60) * (-105) # Multiply by gear ratio
            robot.axis2.motor.Jmotor.phi => deg2rad(20) * (210)
            robot.axis3.motor.Jmotor.phi => deg2rad(90) * (60)
        ], (0.0, 2.0))
        @time "simulation (solve)" sol = solve(prob, Rodas5P(autodiff=false));
        @test SciMLBase.successful_retcode(sol)
    end

    if doplot()
        # plot(sol, layout=30, size=(1900,1200), legend=false)
        @time "Plotting ref" plot(sol, idxs = [
            robot.pathPlanning.controlBus.axisControlBus1.angle_ref# * (-105)
            robot.pathPlanning.controlBus.axisControlBus2.angle_ref# * (210)
            robot.pathPlanning.controlBus.axisControlBus3.angle_ref# * (60)
            robot.pathPlanning.controlBus.axisControlBus4.angle_ref# * (-99)
            robot.pathPlanning.controlBus.axisControlBus5.angle_ref# * (79.2)
            robot.pathPlanning.controlBus.axisControlBus6.angle_ref# * (-99)
        ], layout=9, size=(800,800), l=(:black, :dash), legend=false)
        @time "Plotting ang." plot!(sol, idxs = [
            robot.pathPlanning.controlBus.axisControlBus1.angle
            robot.pathPlanning.controlBus.axisControlBus2.angle
            robot.pathPlanning.controlBus.axisControlBus3.angle
            robot.pathPlanning.controlBus.axisControlBus4.angle
            robot.pathPlanning.controlBus.axisControlBus5.angle
            robot.pathPlanning.controlBus.axisControlBus6.angle
        ], sp=(1:6)')

        plot!(sol, idxs = [
            robot.axis1.motor.Jmotor.phi / ( -105) - robot.pathPlanning.controlBus.axisControlBus1.angle_ref
            robot.axis2.motor.Jmotor.phi / (210) - robot.pathPlanning.controlBus.axisControlBus2.angle_ref
            robot.axis3.motor.Jmotor.phi / (60) - robot.pathPlanning.controlBus.axisControlBus3.angle_ref
        ], sp=(7:9)')
        display(current())

    end

    tv = 0:0.1:2
    angle_ref = sol(tv, idxs=robot.pathPlanning.controlBus.axisControlBus1.angle_ref)
    @test !all(iszero, angle_ref)

    control_error = sol(tv, idxs=robot.pathPlanning.controlBus.axisControlBus1.angle_ref-robot.pathPlanning.controlBus.axisControlBus1.angle)
    @test maximum(abs, control_error) < 0.01
end



@testset "subs constants" begin
    @info "Testing subs constants"
    @named robot = Robot6DOF()
    robot = complete(robot)
    ssys = structural_simplify(IRSystem(robot))
    ssys = Multibody.subs_constants(robot; ssys)
    prob = ODEProblem(ssys, [
        robot.mechanics.r1.phi => deg2rad(-60)
        robot.mechanics.r2.phi => deg2rad(20)
        robot.mechanics.r3.phi => deg2rad(90)
        robot.mechanics.r4.phi => deg2rad(0)
        robot.mechanics.r5.phi => deg2rad(-110)
        robot.mechanics.r6.phi => deg2rad(0)

        robot.axis1.motor.Jmotor.phi => deg2rad(-60) * (-105) # Multiply by gear ratio
        robot.axis2.motor.Jmotor.phi => deg2rad(20) * (210)
        robot.axis3.motor.Jmotor.phi => deg2rad(90) * (60)
    ], (0.0, 2.0))
    sol2 = solve(prob, Rodas5P(autodiff=false))

    @test SciMLBase.successful_retcode(sol2)

    tv = 0:0.1:2
    # control_error = sol2(tv, idxs=robot.pathPlanning.controlBus.axisControlBus1.angle_ref-robot.pathPlanning.controlBus.axisControlBus1.angle)
    # @test maximum(abs, control_error) < 0.002

    # using BenchmarkTools
    # @btime solve(prob, Rodas5P(autodiff=false));
    # 152.225 ms (2272926 allocations: 40.08 MiB)
    # 114.113 ms (1833534 allocations: 33.38 MiB) # sub 0, 1
end
