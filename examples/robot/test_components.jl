using ModelingToolkit
using Multibody
using Test
using JuliaSimCompiler
t = Multibody.t

cd(@__DIR__)
world = Multibody.world
include("OneAxis.jl")
include("FullRobot.jl")
@named structure = MechanicalStructure()



@named motor = Motor()
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Mechanical.Rotational: Fixed
using ModelingToolkitStandardLibrary.Blocks

myfixed(args...; kwargs...) = ModelingToolkitStandardLibrary.Mechanical.Rotational.Fixed(args...; kwargs...)

myinertia(args...; kwargs...) = ModelingToolkitStandardLibrary.Mechanical.Rotational.Inertia(args...; kwargs...)

mytorque(args...; kwargs...) = ModelingToolkitStandardLibrary.Mechanical.Rotational.ConstantTorque(args...; kwargs...)

##

@mtkmodel MotorTest begin
    @components begin
        motor = Motor()
        # fixed = myfixed() # bug in @mtkmodel
        inertia = myinertia(J=1)
        constant = Constant(k=1)
    end
    @equations begin
        # connect(motor.flange_motor, fixed.flange)
        connect(motor.flange_motor, inertia.flange_a)
        constant.output.u ~ motor.axisControlBus.current_ref
    end
end

@named motorTest = MotorTest()
m = structural_simplify(motorTest)
# @test length(states(m)) == 3
    # D(motorTest.motor.gear.bearingFriction.w) => 0
cm = complete(motorTest)

prob = ODEProblem(m, [
    D(D(cm.motor.Jmotor.phi)) => 0
], (0.0, 5.0))
sol = solve(prob, Rodas4())
plot(sol, idxs=cm.motor.phi.phi.u)

##

@mtkmodel GearTest begin
    @components begin
        motor = Motor()
        inertia = myinertia(J=1)
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
m = structural_simplify(gearTest)
cm = complete(gearTest)

prob = ODEProblem(m, [
    ModelingToolkit.missing_variable_defaults(m);
], (0.0, 5.0))
sol = solve(prob, Rodas4())
plot(sol, idxs=cm.motor.phi.phi.u)


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
m = structural_simplify(gearTest)

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
        inertia = myinertia(J=1)
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
m = structural_simplify(axisTest)
# m = structural_simplify(IRSystem(axisTest)) # Yingbo: solution unstable with IRSystem simplification

cm = complete(axisTest)
tspan = (0.0, 5.0)
prob = ODEProblem(m, [
    ModelingToolkit.missing_variable_defaults(m);
    # D(cm.axis2.gear.bearingFriction.w) => 0
], tspan)
sol = solve(prob, Rodas4())
@test SciMLBase.successful_retcode(sol)

@test sol(0.0, idxs=cm.axis2.motor.emf.phi) == 0
# @test sol(tspan[2], idxs=cm.axis2.motor.emf.phi) == 0

isinteractive() && plot(sol, layout=length(states(m)))
isinteractive() && plot(sol, idxs=[
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

@named oneaxis = OneAxis()
oneaxis = complete(oneaxis)
op = Dict([
    oneaxis.axis.flange.phi => 0
    D(oneaxis.axis.flange.phi) => 0
    D(D(oneaxis.axis.flange.phi)) => 0
    D(D(oneaxis.load.phi)) => 0
    oneaxis.axis.controller.PI.T => 0.01
    oneaxis.axis.controller.PI.gainPI.k => 1
    oneaxis.axis.controller.P.k => 10
    oneaxis.load.J => 1.3*15
])
matrices_S, simplified_sys = Blocks.get_sensitivity(oneaxis, :axis₊controller_e; op)


# using ControlSystemsBase 
# S = ss(matrices_S...) |> minreal
# @test isstable(S)
# bodeplot(S)


ssys = structural_simplify(IRSystem(oneaxis)) # Yingbo: solution unstable with IRSystem simplification
ssys = structural_simplify(oneaxis)
# cm = oneaxis
# prob = ODEProblem(ssys, [
#     cm.axis.flange.phi => 0
#     D(cm.axis.flange.phi) => 0
# ], (0.0, 5.0))

prob = ODEProblem(ssys, op, (0.0, 10),)
sol = solve(prob, Rodas4());
isinteractive() && plot(sol, layout=length(states(ssys)))
@test SciMLBase.successful_retcode(sol)
@test sol(10, idxs=oneaxis.axis.controller.PI.err_input.u) ≈ 0 atol=1e-8




@named robot = FullRobot()

# ssys = structural_simplify(robot, allow_parameters = false)
ssys = structural_simplify(IRSystem(robot))



dummyder = setdiff(states(ssys), states(oneaxis))
op = merge(ModelingToolkit.defaults(oneaxis), Dict(x => 0.0 for x in dummyder))


prob = ODEProblem(ssys, op, (0.0, 1.0))

using OrdinaryDiffEq
sol = solve(prob, Rodas4(), u0=prob.u0 .+ 0.0.*randn.())


