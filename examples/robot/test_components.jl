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


##

@mtkmodel MotorTest begin
    @components begin
        motor = Motor()
        fixed = myfixed() # bug in @mtkmodel
        constant = Constant(k=1)
    end
    @equations begin
        connect(motor.flange_motor, fixed.flange)
        constant.output.u ~ motor.axisControlBus.current_ref
    end
end

@named motorTest = MotorTest()
m = structural_simplify(motorTest)
@test length(states(m)) == 3

##

@mtkmodel GearTest begin
    @components begin
        gear2 = GearType2()
        fixed = myfixed() # bug in @mtkmodel
        constant = Constant(k=1)
    end
    @equations begin
        connect(gear2.flange_a, fixed.flange)
    end
end

@named gearTest = GearTest()
m = structural_simplify(gearTest)


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


##
myinertia(args...; kwargs...) = ModelingToolkitStandardLibrary.Mechanical.Rotational.Inertia(args...; kwargs...)

mytorque(args...; kwargs...) = ModelingToolkitStandardLibrary.Mechanical.Rotational.ConstantTorque(args...; kwargs...)

@mtkmodel AxisTest begin
    @components begin
        axis2 = AxisType2()
        # fixed = myfixed() # bug in @mtkmodel
        fixed = mytorque(tau_constant=1, use_support=false) # bug in @mtkmodel
        # inertia = myinertia(J=1)
        # constant1 = Constant(k=1)
        constant2 = Constant(k=1)
        constant3 = Constant(k=1)
        constant4 = Constant(k=1)
        constant5 = Constant(k=1)

        constant6 = Constant(k=1)
        constant7 = Constant(k=1)
        constant8 = Constant(k=1)
        constant9 = Constant(k=1)
        constant10 = Constant(k=1)

    end
    @equations begin
        connect(axis2.flange, fixed.flange)
        # connect(axis2.flange, inertia.flange_a)
        # constant1.output.u ~ motor.axisControlBus.current_ref # This is connected to the controller output
        constant2.output.u ~ motor.axisControlBus.speed_ref
        constant3.output.u ~ motor.axisControlBus.angle_ref
        constant4.output.u ~ motor.axisControlBus.motion_ref
        constant5.output.u ~ motor.axisControlBus.acceleration_ref
    end
end

@named axisTest = AxisTest()
m = structural_simplify(axisTest)
m = structural_simplify(IRSystem(axisTest))

##

@named controller = Controller()
@named axis2 = AxisType2()
@named gear2 = GearType2()
# @named axis1 = AxisType1()
@named gear1 = GearType1()

@named pp = PathPlanning1(;)
@named pp6 = PathPlanning6(;)

@named oneaxis = OneAxis()
@named oneaxis = OneAxis()

ssys = structural_simplify(IRSystem(oneaxis))
ssys = structural_simplify(oneaxis)#, allow_parameters = false)

@named robot = FullRobot()

ssys = structural_simplify(robot, allow_parameters = false)
ssys = structural_simplify(IRSystem(robot))



dummyder = setdiff(states(ssys), states(oneaxis))
op = merge(ModelingToolkit.defaults(oneaxis), Dict(x => 0.0 for x in dummyder))
prob = ODEProblem(ssys, op, (0.0, 10.0))

using OrdinaryDiffEq
sol = solve(prob, Rodas4(), u0=prob.u0 .+ 0.01.*randn.())


