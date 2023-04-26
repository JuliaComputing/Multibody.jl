using ModelingToolkit
using Multibody
import ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical
t = Multibody.t

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
        gain1 = Blocks.Gain(k = ratio)
        PI = Blocks.PI(k = ks, T = Ts)
        feedback1 = Blocks.Feedback()
        P = Blocks.Gain(k = kp)
        add3 = Blocks.Add3(k3 = -1)
        gain2 = Blocks.Gain(k = ratio)
        axisControlBus = AxisControlBus()
    end

    eqs = [connect(gain1.y, feedback1.u1)
           connect(feedback1.y, P.u)
           connect(gain2.y, add3.u1)
           connect(add3.y, PI.u)
           connect(gain2.u, axisControlBus.speed_ref)
           connect(gain1.u, axisControlBus.angle_ref)
           connect(feedback1.u2, axisControlBus.motorAngle)
           connect(add3.u3, axisControlBus.motorSpeed)
           connect(PI.y, axisControlBus.current_ref)]

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

# TODO: MechanicalStructure

# TODO: Motor

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
        Vs = Voltage() # NOTE: should be SignalVoltage
        power = IdealOpAmp()
        diff = IdealOpAmp()
        emf = Electrical.RotationalEMF(; k = k, useSupport = false) # NOTE: not implemented
        La = Inductor(; L = (250 / (2 * D * w)))
        Ra = Resistor(; R = 250)
        Rd2 = Resistor(; R = 100)
        Ri = Resistor(; R = 10)
        Rd1 = Resistor(; R = 100)
        Rp1 = Resistor(; R = 200)
        Rp2 = Resistor(; R = 50)
        Rd4 = Resistor(; R = 100)
        hall2 = Voltage() # NOTE: should be SignalVoltage
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
        convert1 = Gain(k = 1)
        convert2 = Gain(k = 1)
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
           connect(Vs.n, g1.p)
           connect(g2.p, hall2.n)
           connect(Rd4.n, g3.p)
           connect(g3.p, OpI.p1)
           connect(g5.p, Rp2.n)
           connect(emf.n, hall1.p)
           connect(hall1.n, g4.p)
           connect(emf.flange, phi.flange)
           connect(emf.flange, speed.flange)
           connect(OpI.n2, power.n2)
           connect(OpI.p1, OpI.n2)
           connect(OpI.p1, diff.n2)
           connect(Jmotor.flange_b, flange_motor)
           connect(phi.phi, axisControlBus.motorAngle)
           connect(speed.w, axisControlBus.motorSpeed)
           connect(hall1.i, axisControlBus.current)
           connect(hall1.i, convert1.u)
           connect(convert1.y, hall2.v)
           connect(convert2.u, axisControlBus.current_ref)
           connect(convert2.y, Vs.v)
           connect(emf.flange, Jmotor.flange_a)]

    compose(ODESystem(eqs, t; name), systems)
end
