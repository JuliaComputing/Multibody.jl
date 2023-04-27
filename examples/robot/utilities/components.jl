using ModelingToolkit
using Multibody
import ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
t = Multibody.t
D = Differential(t)

function AxisControlBus(; name)
    vars = @variables begin
        (motion_ref(t) = 0), [description = "= true, if reference motion is not in rest"]
        (angle_ref(t) = 0), [description = "Reference angle of axis flange"]
        (angle(t) = 0), [description = "Angle of axis flange"]
        (speed_ref(t) = 0), [description = "Reference speed of axis flange"]
        (speed(t) = 0), [description = "Speed of axis flange"]
        (acceleration_ref(t) = 0), [description = "Reference acceleration of axis flange"]
        (acceleration(t) = 0), [description = "Acceleration of axis flange"]
        (current_ref(t) = 0), [description = "Reference current of motor"]
        (current(t) = 0), [description = "Current of motor"]
        (motorAngle(t) = 0), [description = "Angle of motor flange"]
        (motorSpeed(t) = 0), [description = "Speed of motor flange"]
    end
    ODESystem(Equation[], t, vars, []; name)
end

"""
    AxisType2(; name)

Axis model of the r3 joints 4,5,6
"""
function AxisType2(; name)
    @parameters begin
        kp = 10, [description = "Gain of position controller"]
        ks = 1, [description = "Gain of speed controller"]
        Ts = 0.01, [description = "Time constant of integrator of speed controller"]
        k = 1.1616, [description = "Gain of motor"]
        w = 4590, [description = "Time constant of motor"]
        D = 0.6, [description = "Damping constant of motor"]
        J = 0.0013, [description = "Moment of inertia of motor"]
        ratio = -105, [description = "Gear ratio"]
        Rv0 = 0.4, [description = "Viscous friction torque at zero velocity"]
        Rv1 = 0.13 / 160, [description = "Viscous friction coefficient"]
        peak = 1, [description = "Maximum static friction torque is peak*Rv0 (peak >= 1)"]
    end

    systems = @named begin
        flange_b = Flange()
        gear = GearType2(; Rv0, Rv1, peak, i)
        motor = Motor(; J, k, w, D)
        controller = Controller(; kp, ks, Ts, ratio)
        angleSensor = Rotational.AngleSensor()
        speedSensor = Rotational.SpeedSensor()
        accSensor = Rotational.AccSensor()
        Const = Blocks.Constant(k = 0)
        axisControlBus = AxisControlBus()
    end

    eqs = [connect(gear.flange_b, flange)
           connect(gear.flange_b, angleSensor.flange)
           connect(gear.flange_b, speedSensor.flange)
           connect(motor.flange_motor, gear.flange_a)
           connect(gear.flange_b, accSensor.flange)
           connect(motor.axisControlBus, axisControlBus)
           connect(angleSensor.phi, axisControlBus.angle)
           connect(speedSensor.w, axisControlBus.speed)
           connect(accSensor.a, axisControlBus.acceleration)
           #    connect(axisControlBus.angle_ref, initializeFlange.phi_start)
           #    connect(axisControlBus.speed_ref, initializeFlange.w_start)
           #    connect(initializeFlange.flange, flange)
           #    connect(Const.y, initializeFlange.a_start)
           connect(controller.axisControlBus, axisControlBus)]

    compose(ODESystem(eqs, t; name), systems)
end

function AxisType1(; name)
    @parameters begin
        c = 43, [description = "Spring constant"]
        cd = 0.005, [description = "Damper constant"]
    end
    error("Not implemented")
    # @named axisType2 = AxisType2(redeclare GearType1 gear(c=c, d=cd)) # TODO: Figure out how to handle the redeclare directive https://github.com/SciML/ModelingToolkit.jl/issues/2038
end

function Controller(; name)
    @parameters begin
        kp = 10, [description = "Gain of position controller"]
        ks = 1, [description = "Gain of speed controller"]
        Ts = 0.01, [description = "Time constant of integrator of speed controller"]
        ratio = 1, [description = "Gear ratio of gearbox"]
    end
    systems = @named begin
        gain1 = Blocks.Gain(ratio)
        PI = Blocks.PI(k = ks, T = Ts)
        feedback1 = Blocks.Feedback()
        P = Blocks.Gain(kp)
        add3 = Blocks.Add3(k3 = -1)
        gain2 = Blocks.Gain(ratio)
        axisControlBus = AxisControlBus()
    end

    eqs = [connect(gain1.output, feedback1.input1)
           connect(feedback1.output, P.input)
           connect(gain2.output, add3.input1)
           connect(add3.output, PI.err_input)
           (gain2.input.u ~ axisControlBus.speed_ref)
           (gain1.input.u ~ axisControlBus.angle_ref)
           (feedback1.input2.u ~ axisControlBus.motorAngle)
           (add3.input3.u ~ axisControlBus.motorSpeed)
           (PI.ctr_output.u ~ axisControlBus.current_ref)]

    compose(ODESystem(eqs, t; name), systems)
end

function GearType2(; name)
    #     Modelica.Mechanics.Rotational.Components.BearingFriction bearingFriction(
    #       tau_pos=[0,
    #            Rv0/unitTorque; 1, (Rv0 + Rv1*unitAngularVelocity)/unitTorque], peak=peak,
    #       useSupport=false)

    unitAngularVelocity = 1

    @parameters begin
        i = -99, [description = "Gear ratio"]
        Rv0 = 21.8, [description = "Viscous friction torque at zero velocity"]
        Rv1 = 9.8, [description = "Viscous friction coefficient (R=Rv0+Rv1*abs(qd))"]
        peak = (26.7 / 21.8),
               [description = "Maximum static friction torque is peak*Rv0 (peak >= 1)"]
    end
    systems = @named begin
        flange_a = Flange()
        flange_b = Flange()
        gear = Rotational.IdealGear(; ratio = i, useSupport = false)
        # bearingFriction = Rotational.BearingFriction(; tau_pos=[0, Rv0; 1, (Rv0 + Rv1*unitAngularVelocity)], peak=peak, useSupport=false) # Not yet supported
        bearingFriction = Rotational.RotationalFriction(; f = Rv1, tau_brk = peak * Rv0,
                                                        tau_c = Rv0, w_brk = 0.1) # NOTE: poorly chosen w_brk
    end
    eqs = [connect(gear.flange_b, bearingFriction.flange_a)
           connect(bearingFriction.flange_b, flange_b)
           connect(gear.flange_a, flange_a)]
    compose(ODESystem(eqs, t; name), systems)
end

function GearType1(; name)
    unitAngularVelocity = 1
    unitTorque = 1
    @parameters begin
        i = -105, [description = "Gear ratio"]
        c = 43, [description = "Spring constant"]
        d = 0.005, [description = "Damper constant"]
        Rv0 = 0.4, [description = "Viscous friction torque at zero velocity"]
        Rv1 = (0.13 / 160),
              [description = "Viscous friction coefficient (R=Rv0+Rv1*abs(qd))"]
        peak = 1, [description = "Maximum static friction torque is peak*Rv0 (peak >= 1)"]
    end
    @variables a_rel=D(spring.w_rel) [
        description = "Relative angular acceleration of spring",
    ]

    #   Modelica.Mechanics.Rotational.Components.BearingFriction bearingFriction(
    #     tau_pos=[0,
    #          Rv0/unitTorque; 1, (Rv0 + Rv1*unitAngularVelocity)/unitTorque],
    #       useSupport=false) 

    systems = @named begin
        flange_a = Flange()
        flange_b = Flange()
        gear = Rotational.IdealGear(; ratio = i, useSupport = false)
        spring = Rotational.SpringDamper(; c, d)
        # bearingFriction = Rotational.BearingFriction(; tau_pos=[0, Rv0; 1, (Rv0 + Rv1*unitAngularVelocity)], useSupport=false) # Not yet supported
        bearingFriction = Rotational.RotationalFriction(; f = Rv1, tau_brk = peak * Rv0,
                                                        tau_c = Rv0, w_brk = 0.1) # NOTE: poorly chosen w_brk
    end
    eqs = [connect(spring.flange_b, gear.flange_a)
           connect(bearingFriction.flange_b, spring.flange_a)
           connect(gear.flange_b, flange_b)
           connect(bearingFriction.flange_a, flange_a)]
    compose(ODESystem(eqs, t; name), systems)
end


function Motor(; name)
    @parameters begin
        J = 0.0013, [description = "Moment of inertia of motor"]
        k = 1.1616, [description = "Gain of motor"]
        w = 4590, [description = "Time constant of motor"]
        D = 0.6, [description = "Damping constant of motor"]
        w_max = 315, [description = "Maximum speed of motor"]
        i_max = 9, [description = "Maximum current of motor"]
    end

    #   Modelica.Electrical.Analog.Sources.SignalVoltage Vs
    #   Electrical.Analog.Basic.RotationalEMF emf(k=k, useSupport=false)

    systems = @named begin
        flange_motor = Rotational.Flange()
        Vs = Voltage() # was SignalVoltage
        power = IdealOpAmp()
        diff = IdealOpAmp()
        emf = Electrical.EMF(; k = k) # Was RotationalEMF(k, useSupport=false), which differs from EMF in that it can be instantiated without support
        La = Inductor(; L = (250 / (2 * D * w)))
        Ra = Resistor(; R = 250)
        Rd2 = Resistor(; R = 100)
        C = Capacitor(C=0.004*D/w)
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
        Jmotor = Rotational.Inertia(; J = J)
        axisControlBus = AxisControlBus()
        convert1 = Blocks.Gain(1)
        convert2 = Blocks.Gain(1)
    end

    eqs = [connect(La.n, emf.p)
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

function MechanicalStructure(; name)
    # parameter Boolean animation=true "= true, if animation shall be enabled";
    # parameter SI.Mass mLoad(min=0)=15 "Mass of load";
    # parameter SI.Position rLoad[3]={0,0.25,0}
    #   "Distance from last flange to load mass";
    # parameter SI.Acceleration g=9.81 "Gravity acceleration";

    @parameters begin
        mLoad = 15, [description = "Mass of load"]
        rLoad[1:3] = [0, 0.25, 0], [description = "Distance from last flange to load mass"]
        g = 9.81, [description = "Gravity acceleration"]
    end

    @variables begin
        (q(t)[1:6] = 0), [description = "Joint angles"]
        (qd(t)[1:6] = 0), [description = "Joint speeds"]
        (qdd(t)[1:6] = 0), [description = "Joint accelerations"]
        (tau(t)[1:6] = 0), [description = "Joint driving torques"]
    end

    systems = @named begin
        axis1 = Rotational.Flange()
        axis2 = Rotational.Flange()
        axis3 = Rotational.Flange()
        axis4 = Rotational.Flange()
        axis5 = Rotational.Flange()
        axis6 = Rotational.Flange()
        r1 = Revolute(n = [0, 1, 0], useAxisFlange = true)
        r2 = Revolute(n = [1, 0, 0], useAxisFlange = true)
        r3 = Revolute(n = [1, 0, 0], useAxisFlange = true)
        r4 = Revolute(n = [0, 1, 0], useAxisFlange = true)
        r5 = Revolute(n = [1, 0, 0], useAxisFlange = true)
        r6 = Revolute(n = [0, 1, 0], useAxisFlange = true)
        b0 = BodyShape(r = [0, 0.351, 0],
                    #    r_shape = [0, 0, 0],
                       #    lengthDirection = [1, 0, 0],
                       #    widthDirection = [0, 1, 0],
                       #    length = 0.225,
                       #    width = 0.3,
                       #    height = 0.3,
                       r_cm = [0, 0, 0],
                       m = 1)
        b1 = BodyShape(r = [0, 0.324, 0.3],
                       I_22 = 1.16,
                       #    lengthDirection = [1, 0, 0],
                       #    widthDirection = [0, 1, 0],
                       #    length = 0.25,
                       #    width = 0.15,
                       #    height = 0.2,
                       #    color = [255, 0, 0],
                       r_cm = [0, 0, 0],
                       m = 1)
        b2 = BodyShape(r = [0, 0.65, 0],
                       r_cm = [0.172, 0.205, 0],
                       m = 56.5,
                       I_11 = 2.58,
                       I_22 = 0.64,
                       I_33 = 2.73,
                       I_21 = -0.46
                       #    lengthDirection = [1, 0, 0],
                       #    widthDirection = [0, 1, 0],
                       #    length = 0.5,
                       #    width = 0.2,
                       #    height = 0.15,
                       #    color = [255, 178, 0],
                       )
        b3 = BodyShape(r = [0, 0.414, -0.155],
                       r_cm = [0.064, -0.034, 0],
                       m = 26.4,
                       I_11 = 0.279,
                       I_22 = 0.245,
                       I_33 = 0.413,
                       I_21 = -0.070
                       #    lengthDirection = [1, 0, 0],
                       #    widthDirection = [0, 1, 0],
                       #    length = 0.15,
                       #    width = 0.15,
                       #    height = 0.15,
                       #    color = [255, 0, 0],
                       )
        b4 = BodyShape(r = [0, 0.186, 0],
                       r_cm = [0, 0, 0],
                       m = 28.7,
                       I_11 = 1.67,
                       I_22 = 0.081,
                       I_33 = 1.67
                       #   lengthDirection = [1, 0, 0],
                       #   widthDirection = [0, 1, 0],
                       #   length = 0.73,
                       #   width = 0.1,
                       #   height = 0.1,
                       #   color = [255, 178, 0],
                       )
        b5 = BodyShape(r = [0, 0.125, 0],
                       r_cm = [0, 0, 0],
                       m = 5.2,
                       I_11 = 1.25,
                       I_22 = 0.81,
                       I_33 = 1.53
                       # lengthDirection = [1, 0, 0],
                       # widthDirection = [0, 1, 0],
                       # length = 0.225,
                       # width = 0.075,
                       # height = 0.1,
                       # color = [0, 0, 255],
                       )
        b6 = BodyShape(r = [0, 0, 0],
                       r_cm = [0.05, 0.05, 0.05],
                       m = 0.5
                       # lengthDirection = [1, 0, 0],
                       # widthDirection = [0, 1, 0],
                       # color = [0, 0, 255],
                       )
        load = BodyShape(r = [0, 0, 0],
                         r_cm = rLoad,
                         m = mLoad
                         # widthDirection = [1, 0, 0],
                         # width = 0.05,
                         # height = 0.05,
                         # color = [255, 0, 0],
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

    compose(ODESystem(eqs, t; name), systems)
end





@named structure = MechanicalStructure()
@named motor = Motor()
@named controller = Controller()