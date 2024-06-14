using ModelingToolkit
import ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
t = Multibody.t
D = Differential(t)

"""
    @connector AxisControlBus(; name)    

- `motion_ref(t) = 0`: = true, if reference motion is not in rest
- `angle_ref(t) = 0`: Reference angle of axis flange
- `angle(t) = 0`: Angle of axis flange
- `speed_ref(t) = 0`: Reference speed of axis flange
- `speed(t) = 0`: Speed of axis flange
- `acceleration_ref(t) = 0`: Reference acceleration of axis flange
- `acceleration(t) = 0`: Acceleration of axis flange
- `current_ref(t) = 0`: Reference current of motor
- `current(t) = 0`: Current of motor
- `motorAngle(t) = 0`: Angle of motor flange
- `motorSpeed(t) = 0`: Speed of motor flange
"""
@connector function AxisControlBus(; name)
    vars = @variables begin
        (motion_ref(t)), [guess=0.0,input = true, description = "= true, if reference motion is not in rest"]
        (angle_ref(t)), [guess=0.0,input = true, description = "Reference angle of axis flange"]
        (angle(t)), [guess=0.0,output = true, description = "Angle of axis flange"]
        (speed_ref(t)), [guess=0.0,input = true, description = "Reference speed of axis flange"]
        (speed(t)), [guess=0.0,output = true, description = "Speed of axis flange"]
        (acceleration_ref(t)), [guess=0.0,input = true, description = "Reference acceleration of axis flange"]
        (acceleration(t)), [guess=0.0,output = true, description = "Acceleration of axis flange"]
        (current_ref(t)), [guess=0.0,input = true, description = "Reference current of motor"]
        (current(t)), [guess=0.0,output = true, description = "Current of motor"]
        (motorAngle(t)), [guess=0.0,output = true, description = "Angle of motor flange"]
        (motorSpeed(t)), [guess=0.0,output = true, description = "Speed of motor flange"]
    end
    ODESystem(Equation[], t, vars, []; name)
end

@connector function ControlBus(; name)
    systems = @named begin
        axisControlBus1 = AxisControlBus()
        axisControlBus2 = AxisControlBus()
        axisControlBus3 = AxisControlBus()
        axisControlBus4 = AxisControlBus()
        axisControlBus5 = AxisControlBus()
        axisControlBus6 = AxisControlBus()
    end
    ODESystem(Equation[], t; systems, name)
end

"""
    AccSensor(;name)

Ideal sensor to measure the absolute flange angular acceleration

# Connectors:

  - `flange`: [Flange](@ref) Flange of shaft from which sensor information shall be measured
  - `a`: [RealOutput](@ref) Absolute angular acceleration of flange
"""
@component function AccSensor(; name)
    @named flange = Rotational.Flange()
    @variables w(t) [description = "Absolute angular velocity of flange"]
    @named a = Blocks.RealOutput() #[description = "Absolute angular acceleration of flange"]
    eqs = [D(flange.phi) ~ w
           a.u ~ D(w)
           flange.tau ~ 0
           ]
    return ODESystem(eqs, t, [], []; name = name, systems = [flange, a])
end

RotationalFlange = Rotational.Flange
RotationalAngleSensor = Rotational.AngleSensor
RotationalSpeedSensor = Rotational.SpeedSensor

# @mtkmodel AxisType2 begin #(; name, kp = 10, ks = 1, Ts = 0.01, k = 1.1616, w = 4590, D = 0.6,
#     #J = 0.0013, ratio = -105, Rv0 = 0.4, Rv1 = 0.13 / 160, peak = 1)
# @parameters begin
# kp = 10, [description = "Gain of position controller"]
# ks = 1, [description = "Gain of speed controller"]
# Ts = 0.01, [description = "Time constant of integrator of speed controller"]
# k = 1.1616, [description = "Gain of motor"]
# w = 4590, [description = "Time constant of motor"]
# D = 0.6, [description = "Damping constant of motor"]
# J = 0.0013, [description = "Moment of inertia of motor"]
# ratio = -105, [description = "Gear ratio"]
# Rv0 = 0.4, [description = "Viscous friction torque at zero velocity"]
# Rv1 = 0.13 / 160, [description = "Viscous friction coefficient"]
# peak = 1,
# [description = "Maximum static friction torque is peak*Rv0 (peak >= 1)"]
# end

"""
    AxisType2(; name)

Axis model of the r3 joints 4,5,6
"""
function AxisType2(; name, kp = 10, ks = 1, Ts = 0.01, k = 1.1616, w = 4590, D = 0.6,
                   J = 0.0013, ratio = -105, Rv0 = 0.4, Rv1 = 0.13 / 160, peak = 1)
    # pars = @parameters begin
    #     kp = kp, [description = "Gain of position controller"]
    #     ks = ks, [description = "Gain of speed controller"]
    #     Ts = Ts, [description = "Time constant of integrator of speed controller"]
    #     k = k, [description = "Gain of motor"]
    #     w = w, [description = "Time constant of motor"]
    #     D = D, [description = "Damping constant of motor"]
    #     J = J, [description = "Moment of inertia of motor"]
    #     # ratio = ratio, [description = "Gear ratio"]
    #     Rv0 = Rv0, [description = "Viscous friction torque at zero velocity"]
    #     Rv1 = Rv1, [description = "Viscous friction coefficient"]
    #     peak = peak,
    #            [description = "Maximum static friction torque is peak*Rv0 (peak >= 1)"]
    # end

    systems = @named begin
        flange = Rotational.Flange()
        gear = GearType2(; Rv0, Rv1, peak, i = ratio)
        motor = Motor(; J, k, w, D)
        controller = Controller(; kp, ks, Ts, ratio)
        angleSensor = Rotational.AngleSensor()
        speedSensor = Rotational.SpeedSensor()
        accSensor = AccSensor() # TODO: shift to MTKstdlib version when merged
        # Const = Blocks.Constant(k = 0)
        axisControlBus = AxisControlBus()
    end

    eqs = [
           connect(flange, gear.flange_b)
           connect(motor.flange_motor, gear.flange_a, angleSensor.flange, speedSensor.flange, accSensor.flange)


           connect(motor.axisControlBus, axisControlBus)
           (angleSensor.phi.u/ratio ~ axisControlBus.angle)
           (speedSensor.w.u/ratio ~ axisControlBus.speed)
           (accSensor.a.u/ratio ~ axisControlBus.acceleration)
           connect(controller.axisControlBus, axisControlBus)]

    ODESystem(eqs, t; name, systems)
end

function AxisType1(; name, c = 43, cd = 0.005, kp = 10, ks = 1, Ts = 0.01, k = 1.1616, w = 4590, D = 0.6,
    J = 0.0013, ratio = -105, Rv0 = 0.4, Rv1 = 0.13 / 160, peak = 1)
    # @parameters begin
    #     c = c, [description = "Spring constant"]
    #     cd = cd, [description = "Damper constant"]
    # end

    systems = @named begin
        flange = Rotational.Flange()
        spring = Rotational.SpringDamper(c=c, d=cd)
        gear = GearType2(; Rv0, Rv1, peak, i = ratio)
        motor = Motor(; J, k, w, D)
        controller = Controller(; kp, ks, Ts, ratio)
        angleSensor = Rotational.AngleSensor()
        speedSensor = Rotational.SpeedSensor()
        accSensor = AccSensor() # TODO: shift to MTKstdlib version when merged
        # Const = Blocks.Constant(k = 0)
        axisControlBus = AxisControlBus()
    end

    eqs = [
        connect(flange, gear.flange_b)
        connect(gear.flange_a, spring.flange_b)
        connect(motor.flange_motor, spring.flange_a, angleSensor.flange, speedSensor.flange, accSensor.flange)
        connect(motor.axisControlBus, axisControlBus)
        (angleSensor.phi.u/ratio ~ axisControlBus.angle)
        (speedSensor.w.u/ratio ~ axisControlBus.speed)
        (accSensor.a.u/ratio ~ axisControlBus.acceleration)
        connect(controller.axisControlBus, axisControlBus)
    ]

    ODESystem(eqs, t; name, systems)
end

function Controller(; name, kp = 10, ks = 1, Ts = 0.01, ratio = 1)
    # pars = @parameters begin
    #     kp = kp, [description = "Gain of position controller"]
    #     ks = ks, [description = "Gain of speed controller"]
    #     Ts = Ts, [description = "Time constant of integrator of speed controller"]
    #     ratio = ratio, [description = "Gear ratio of gearbox"]
    # end
    systems = @named begin
        gain1 = Blocks.Gain(ratio)
        PI = Blocks.PI(gainPI.k = ks, T = Ts)
        feedback1 = Blocks.Feedback()
        P = Blocks.Gain(kp)
        add3 = Blocks.Add3(k3 = -1)
        gain2 = Blocks.Gain(ratio)
        axisControlBus = AxisControlBus()
    end

    eqs = [connect(gain1.output, feedback1.input1)
           connect(feedback1.output, P.input)
           connect(P.output, add3.input2)
           connect(gain2.output, add3.input1)
           connect(add3.output, :e, PI.err_input)
           (gain2.input.u ~ axisControlBus.speed_ref)
           (gain1.input.u ~ axisControlBus.angle_ref)
           (feedback1.input2.u ~ axisControlBus.motorAngle)
           (add3.input3.u ~ axisControlBus.motorSpeed)
           (PI.ctr_output.u ~ axisControlBus.current_ref)]

    ODESystem(eqs, t; name, systems)
end

function GearType2(; name, i = -99,
                   Rv0 = 21.8,
                   Rv1 = 9.8,
                   peak = (26.7 / 21.8))
    #     Modelica.Mechanics.Rotational.Components.BearingFriction bearingFriction(
    #       tau_pos=[0,
    #            Rv0/unitTorque; 1, (Rv0 + Rv1*unitAngularVelocity)/unitTorque], peak=peak,
    #       useSupport=false)

    unitAngularVelocity = 1

    # @parameters begin
    #     i = i, [description = "Gear ratio"]
    #     Rv0 = Rv0, [description = "Viscous friction torque at zero velocity"]
    #     Rv1 = Rv1, [description = "Viscous friction coefficient (R=Rv0+Rv1*abs(qd))"]
    #     peak = peak,
    #            [description = "Maximum static friction torque is peak*Rv0 (peak >= 1)"]
    # end
    systems = @named begin
        flange_a = Rotational.Flange()
        flange_b = Rotational.Flange()
        gear = Rotational.IdealGear(; ratio = i, use_support = false)
        # bearingFriction = Rotational.BearingFriction(; tau_pos=[0, Rv0; 1, (Rv0 + Rv1*unitAngularVelocity)], peak=peak, useSupport=false) # Not yet supported
        # bearingFriction = BearingFriction(; f = Rv1, tau_brk = peak * Rv0, tau_c = Rv0, w_brk = 0.1) # NOTE: poorly chosen w_brk            
        bearingFriction = BearingFriction(;)
    end
    #=
    NOTE: We do not yet have the bearingFriction component, bearingFriction this component extends PartialElementaryTwoFlangesAndSupport2 which is implicitly grounded when use_support=false. This component has a relative angle state. Instead, we use a RotationalFriction component, which extends PartialCompliantWithRelativeStates that does not have implicit grounding. We therefore add the explicit grounding using a fixed component
    =#
    eqs = [
        connect(flange_a, gear.flange_a)
        connect(gear.flange_b, bearingFriction.flange_a)
        connect(bearingFriction.flange_b, flange_b)


        # Equations below are the save as above, but without the gear
        # connect(bearingFriction.flange_b, flange_b)
        # connect(bearingFriction.flange_a, flange_a)
    ]
    ODESystem(eqs, t; name, systems)
end

using ModelingToolkitStandardLibrary.Mechanical.Rotational: PartialElementaryTwoFlangesAndSupport2
import ModelingToolkitStandardLibrary.Mechanical.Rotational.Flange as fl
@mtkmodel BearingFriction begin
    # @extend flange_a, flange_b, phi_support = partial_comp = PartialElementaryTwoFlangesAndSupport2(;use_support = false)
    # @parameters begin
    #     f, [description = "Viscous friction coefficient"]
    #     tau_c, [description = "Coulomb friction torque"]
    #     w_brk, [description = "Breakaway friction velocity"]
    #     tau_brk, [description = "Breakaway friction torque"]
    # end
    # @variables begin
    #     phi(t) = 0.0, [description = "Angle between shaft flanges (flange_a, flange_b) and support"]
    #     tau(t) = 0.0, [description = "Torque between flanges"]
    #     w(t) = 0.0
    #     a(t) = 0.0
    # end

    # begin
    #     str_scale = sqrt(2 * exp(1)) * (tau_brk - tau_c)
    #     w_st = w_brk * sqrt(2)
    #     w_coul = w_brk / 10
    # end
    # @equations begin
    #     tau ~ str_scale * (exp(-(w / w_st)^2) * w / w_st) +
    #           tau_c * tanh(w / w_coul) + f * w # Stribeck friction + Coulomb friction + Viscous friction

    #     phi ~ flange_a.phi - phi_support;
    #     flange_b.phi ~ flange_a.phi;
    
    #     # Angular velocity and angular acceleration of flanges
    #     w ~ D(phi)
    #     a ~ D(w)

    #     flange_a.tau + flange_b.tau - tau ~ 0
    # end
    @components begin
        flange_a = fl()
        flange_b = fl()
    end
    @equations begin
        connect(flange_a, flange_b)
    end
end


function GearType1(; name, i = -105, c = 43, d = 0.005,
                   Rv0 = 0.4,
                   Rv1 = (0.13 / 160),
                   peak = 1)
    unitAngularVelocity = 1
    unitTorque = 1
    # pars = @parameters begin
    #     i = i, [description = "Gear ratio"]
    #     c = c, [description = "Spring constant"]
    #     d = d, [description = "Damper constant"]
    #     Rv0 = Rv0, [description = "Viscous friction torque at zero velocity"]
    #     Rv1 = Rv1,
    #           [description = "Viscous friction coefficient (R=Rv0+Rv1*abs(qd))"]
    #     peak = peak,
    #            [description = "Maximum static friction torque is peak*Rv0 (peak >= 1)"]
    # end

    #   Modelica.Mechanics.Rotational.Components.BearingFriction bearingFriction(
    #     tau_pos=[0,
    #          Rv0/unitTorque; 1, (Rv0 + Rv1*unitAngularVelocity)/unitTorque],
    #       useSupport=false) 

    systems = @named begin
        flange_a = Rotational.Flange()
        flange_b = Rotational.Flange()
        gear = Rotational.IdealGear(; ratio = i, use_support = false)
        spring = Rotational.SpringDamper(; c, d)
        # bearingFriction = Rotational.BearingFriction(; tau_pos=[0, Rv0; 1, (Rv0 + Rv1*unitAngularVelocity)], useSupport=false) # Not yet supported
        bearingFriction = Rotational.RotationalFriction(; f = Rv1, tau_brk = peak * Rv0,
                                                        tau_c = Rv0, w_brk = 0.1) # NOTE: poorly chosen w_brk
    end
    # vars = @variables a_rel(t)=D(spring.w_rel) [ # This is only used inside "initial equation" block
    #     description = "Relative angular acceleration of spring",
    # ]

    eqs = [connect(spring.flange_b, gear.flange_a)
           connect(bearingFriction.flange_b, spring.flange_a)
           connect(gear.flange_b, flange_b)
           connect(bearingFriction.flange_a, flange_a)]
    ODESystem(eqs, t; name, systems)
end

function Motor(; name, J = 0.0013, k = 1.1616, w = 4590, D = 0.6, w_max = 315, i_max = 9)
    # @parameters begin
    #     J = J, [description = "Moment of inertia of motor"]
    #     k = k, [description = "Gain of motor"]
    #     w = w, [description = "Time constant of motor"]
    #     D = D, [description = "Damping constant of motor"]
    #     w_max = w_max, [description = "Maximum speed of motor"]
    #     i_max = i_max, [description = "Maximum current of motor"]
    # end

    #   Electrical.Analog.Basic.RotationalEMF emf(k=k, useSupport=false)

    systems = @named begin
        flange_motor = Rotational.Flange()
        Vs = Voltage() # was SignalVoltage
        power = IdealOpAmp()
        diff = IdealOpAmp()
        fixed = Rotational.Fixed() # NOTE: this was added to account for the fixed frame that is added to RotationalEMF when useSupport=false
        emf = Electrical.EMF(; k = k) # Was RotationalEMF(k, useSupport=false), which differs from EMF in that it can be instantiated without support
        La = Inductor(; L = (250 / (2 * D * w)), i=0)
        Ra = Resistor(; R = 250)
        Rd2 = Resistor(; R = 100)
        C = Capacitor(C = 0.004 * D / w, v=0)
        OpI = IdealOpAmp()
        Ri = Resistor(; R = 10)
        Rd1 = Resistor(; R = 100)
        Rp1 = Resistor(; R = 200)
        Rp2 = Resistor(; R = 50)
        Rd4 = Resistor(; R = 100)
        hall2 = Voltage() # was SignalVoltage
        Rd3 = Resistor(; R = 100)
        g1 = Ground()
        g2 = Ground()
        g3 = Ground()
        hall1 = CurrentSensor()
        g4 = Ground()
        g5 = Ground()
        phi = Rotational.AngleSensor()
        speed = Rotational.SpeedSensor()
        Jmotor = Rotational.Inertia(; J = J, w=0)
        axisControlBus = AxisControlBus()
        convert1 = Blocks.Gain(1)
        convert2 = Blocks.Gain(1)
    end

    eqs = [connect(fixed.flange, emf.support) # NOTE: extra equation added 
           connect(La.n, emf.p)
           connect(Ra.n, La.p)
           connect(Rd2.n, diff.n1)
           connect(C.n, OpI.p2)
           connect(OpI.p2, power.p1)
           connect(Vs.p, Rd2.p)
           connect(diff.n1, Rd1.p)
           connect(Rd1.n, diff.p2)
           connect(diff.p2, Ri.p)
           connect(Ri.n, OpI.n1)
           connect(OpI.n1, C.p)
           connect(power.n1, Rp1.p)
           connect(power.p2, Rp1.n)
           connect(Rp1.p, Rp2.p)
           connect(power.p2, Ra.p)
           connect(Rd3.p, hall2.p)
           connect(Rd3.n, diff.p1)
           connect(Rd3.n, Rd4.p)
           connect(Vs.n, g1.g)
           connect(g2.g, hall2.n)
           connect(Rd4.n, g3.g)
           connect(g3.g, OpI.p1)
           connect(g5.g, Rp2.n)
           connect(emf.n, hall1.p)
           connect(hall1.n, g4.g)
           connect(emf.flange, phi.flange)
           connect(emf.flange, speed.flange)
           connect(OpI.n2, power.n2)
           connect(OpI.p1, OpI.n2)
           connect(OpI.p1, diff.n2)
           connect(Jmotor.flange_b, flange_motor)
           (phi.phi.u ~ axisControlBus.motorAngle)
           (speed.w.u ~ axisControlBus.motorSpeed)
           (hall1.i ~ axisControlBus.current)
           (hall1.i ~ convert1.u)
           (convert1.y ~ hall2.v)
           (convert2.u ~ axisControlBus.current_ref)
           (convert2.y ~ Vs.v)
           connect(emf.flange, Jmotor.flange_a)]

    compose(ODESystem(eqs, t; name), systems)
end

robot_orange = [1, 0.51, 0, 1]

function MechanicalStructure(; name, mLoad = 15, rLoad = [0, 0.25, 0], g = 9.81)
    # @parameters begin
    #     mLoad = mLoad, [description = "Mass of load"]
    #     rLoad[1:3] = rLoad, [description = "Distance from last flange to load mass"]
    #     g = g, [description = "Gravity acceleration"]
    # end

    @variables begin
        (q(t)[1:6]), [guess = 0, state_priority = typemax(Int), description = "Joint angles"]
        (qd(t)[1:6]), [guess = 0, state_priority = typemax(Int), description = "Joint speeds"]
        (qdd(t)[1:6]),
        [guess = 0, state_priority = typemax(Int), description = "Joint accelerations"]
        (tau(t)[1:6]),
        [guess = 0, state_priority = typemax(Int), description = "Joint driving torques"]
    end
    path = @__DIR__()

    systems = @named begin
        axis1 = Rotational.Flange()
        axis2 = Rotational.Flange()
        axis3 = Rotational.Flange()
        axis4 = Rotational.Flange()
        axis5 = Rotational.Flange()
        axis6 = Rotational.Flange()
        r1 = Revolute(n = [0, 1, 0], axisflange = true, isroot = false, radius=0.12, color=robot_orange)
        r2 = Revolute(n = [1, 0, 0], axisflange = true, isroot = false, radius=0.1, color=robot_orange)
        r3 = Revolute(n = [1, 0, 0], axisflange = true, isroot = false, radius=0.075, color=robot_orange)
        r4 = Revolute(n = [0, 1, 0], axisflange = true, isroot = false, radius=0.06, color=robot_orange)
        r5 = Revolute(n = [1, 0, 0], axisflange = true, isroot = false, radius=0.05, color=robot_orange)
        r6 = Revolute(n = [0, 1, 0], axisflange = true, isroot = false, radius=0.02, color=[0.5, 0.5, 0.5, 1])
        b0 = BodyShape(r = [0, 0.351, 0],
                        shapefile = joinpath(path, "../../examples/resources/b0.stl"),
                       #    r_shape = [0, 0, 0],
                       #    lengthDirection = [1, 0, 0],
                       #    widthDirection = [0, 1, 0],
                       #    length = 0.225,
                       #    width = 0.3,
                       #    height = 0.3,
                       radius = 0.3/2,
                       color = [0.5, 0.5, 0.5, 1],
                       m = 1)
        b1 = BodyShape(r = [0, 0.324, 0.3],
                        shapefile = joinpath(path, "../../examples/resources/b1.stl"),
                       I_22 = 1.16,
                       #    lengthDirection = [1, 0, 0],
                       #    widthDirection = [0, 1, 0],
                       #    length = 0.25,
                       #    width = 0.15,
                       #    height = 0.2,
                       radius = 0.2/2,
                       color = robot_orange,
                       m = 1)
        b2 = BodyShape(r = [0, 0.65, 0],
                        shapefile = joinpath(path, "../../examples/resources/b2.stl"),
                       r_cm = [0.172, 0.205, 0],
                       m = 56.5,
                       I_11 = 2.58,
                       I_22 = 0.64,
                       I_33 = 2.73,
                       I_21 = -0.46,
                       #    lengthDirection = [1, 0, 0],
                       #    widthDirection = [0, 1, 0],
                       #    length = 0.5,
                       #    width = 0.2,
                       #    height = 0.15,
                       radius = 0.2/2,
                       color = robot_orange,
                       )
        b3 = BodyShape(r = [0, 0.414, -0.155],
                        shapefile = joinpath(path, "../../examples/resources/b3.stl"),
                       r_cm = [0.064, -0.034, 0],
                       m = 26.4,
                       I_11 = 0.279,
                       I_22 = 0.245,
                       I_33 = 0.413,
                       I_21 = -0.070,
                       #    lengthDirection = [1, 0, 0],
                       #    widthDirection = [0, 1, 0],
                       #    length = 0.15,
                       #    width = 0.15,
                       #    height = 0.15,
                       radius = 0.15/2,
                       color = robot_orange,
                       )
        b4 = BodyShape(r = [0, 0.186, 0],
                        shapefile = joinpath(path, "../../examples/resources/b4.stl"),
                       m = 28.7,
                       I_11 = 1.67,
                       I_22 = 0.081,
                       I_33 = 1.67,
                       #   lengthDirection = [1, 0, 0],
                       #   widthDirection = [0, 1, 0],
                       #   length = 0.73,
                       #   width = 0.1,
                       #   height = 0.1,
                       radius = 0.1/2,
                       color = robot_orange,
                       )
        b5 = BodyShape(r = [0, 0.125, 0],
                        shapefile = joinpath(path, "../../examples/resources/b5.stl"),
                       m = 5.2,
                       I_11 = 1.25,
                       I_22 = 0.81,
                       I_33 = 1.53,
                       # lengthDirection = [1, 0, 0],
                       # widthDirection = [0, 1, 0],
                       # length = 0.225,
                       # width = 0.075,
                       # height = 0.1,
                        radius = 0.1/2,
                        color = [0.5, 0.5, 0.5, 1],
                       )
        b6 = BodyShape(r = [0, 0, 0],
                        shapefile = joinpath(path, "../../examples/resources/b6.stl"),
                       r_cm = [0.05, 0.05, 0.05],
                       m = 0.5,
                       # lengthDirection = [1, 0, 0],
                       # widthDirection = [0, 1, 0],
                       color = [0.5, 0.5, 0.5, 1],
                       )
        load = BodyShape(r = [0, 0, 0],
                         r_cm = rLoad,
                         m = mLoad,
                         # widthDirection = [1, 0, 0],
                         # width = 0.05,
                         # height = 0.05,
                         color = [1, 0, 0, 1],
                         # lengthDirection = to_unit1(rLoad),
                         # length = length(rLoad),
                         )
    end
    eqs = [q .~ [r1.phi, r2.phi, r3.phi, r4.phi, r5.phi, r6.phi]
           qd .~ D.(q)
           qdd .~ D.(qd)
           tau .~ [r1.tau, r2.tau, r3.tau, r4.tau, r5.tau, r6.tau]
           connect(load.frame_a, b6.frame_b)
           connect(world.frame_b, b0.frame_a)
           connect(b0.frame_b, r1.frame_a)
           connect(b1.frame_b, r2.frame_a)
           connect(r1.frame_b, b1.frame_a)
           connect(r2.frame_b, b2.frame_a)
           connect(b2.frame_b, r3.frame_a)
           connect(r2.axis, axis2)
           connect(r1.axis, axis1)
           connect(r3.frame_b, b3.frame_a)
           connect(b3.frame_b, r4.frame_a)
           connect(r3.axis, axis3)
           connect(r4.axis, axis4)
           connect(r4.frame_b, b4.frame_a)
           connect(b4.frame_b, r5.frame_a)
           connect(r5.axis, axis5)
           connect(r5.frame_b, b5.frame_a)
           connect(b5.frame_b, r6.frame_a)
           connect(r6.axis, axis6)
           connect(r6.frame_b, b6.frame_a)]

    compose(ODESystem(eqs, t; name), [world; systems])
end
