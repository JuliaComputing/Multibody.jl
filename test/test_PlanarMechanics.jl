# using Revise
# using Plots
using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkitStandardLibrary.Blocks
import Multibody.PlanarMechanics as Pl
using JuliaSimCompiler

tspan = (0.0, 3.0)
g = -9.807

@testset "Free body" begin
    m = 2
    I = 1
    @named body = Pl.Body(; m, I)
    @named model = ODESystem(Equation[],
        t,
        [],
        [],
        systems = [body])
    sys = structural_simplify(IRSystem(model))
    unset_vars = setdiff(unknowns(sys), keys(ModelingToolkit.defaults(sys)))
    prob = ODEProblem(sys, unset_vars .=> 0.0, tspan)

    sol = solve(prob, Rodas5P(), initializealg=BrownFullBasicInit())
    @test SciMLBase.successful_retcode(sol)

    free_falling_displacement = 0.5 * g * tspan[end]^2  # 0.5 * g * t^2
    @test sol[body.r[2]][end] ≈ free_falling_displacement
    @test sol[body.r[1]][end] == 0  # no horizontal displacement
    @test all(sol[body.phi] .== 0)
    # plot(sol, idxs = [body.r[1], body.r[2]])
end

@testset "Pendulum" begin
    @named ceiling = Pl.Fixed()
    @named rod = Pl.FixedTranslation(r = [1.0, 0.0])
    @named body = Pl.Body(m = 1, I = 0.1)
    @named revolute = Pl.Revolute()

    connections = [
        connect(ceiling.frame_b, revolute.frame_a),
        connect(revolute.frame_b, rod.frame_a),
        connect(rod.frame_b, body.frame_a)
    ]

    @named model = ODESystem(connections,
        t,
        systems = [body, revolute, rod, ceiling])
    model = complete(model)
    ssys = structural_simplify(IRSystem(model))

    @test length(unknowns(ssys)) == 2
    unset_vars = setdiff(unknowns(ssys), keys(ModelingToolkit.defaults(ssys)))
    prob = ODEProblem(ssys, unset_vars .=> 0.0, tspan)

    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)
    @test sol(1, idxs=model.rod.frame_a.phi) ≈ -2.881383661312169 atol=1e-2
    @test sol(2, idxs=model.rod.frame_a.phi) ≈ -1 atol=1e-2
end

@testset "Pendulum with body shape" begin
    @named ceiling = Pl.Fixed()
    @named rod = Pl.BodyShape(r = [1.0, 0.0], m=1, I=0.1)
    @named revolute = Pl.Revolute()

    connections = [
        connect(ceiling.frame_b, revolute.frame_a),
        connect(revolute.frame_b, rod.frame_a),
    ]

    @named model = ODESystem(connections,
        t,
        systems = [revolute, rod, ceiling])
    model = complete(model)
    ssys = structural_simplify(IRSystem(model))

    @test length(unknowns(ssys)) == 2
    unset_vars = setdiff(unknowns(ssys), keys(ModelingToolkit.defaults(ssys)))
    prob = ODEProblem(ssys, unset_vars .=> 0.0, tspan)

    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)
    @test sol(1, idxs=model.rod.frame_a.phi) ≈ -pi atol=1e-2
    @test sol(2, idxs=model.rod.frame_a.phi) ≈ 0 atol=1e-2
end

@testset "Prismatic" begin
    # just testing instantiation
    @test_nowarn @named prismatic = Pl.Prismatic(r = [1.0, 0.0])
end

@testset "AbsoluteAccCentrifugal" begin
    m = 1
    I = 0.1
    w = 10
    resolve_in_frame = :world

    # components
    @named body = Pl.Body(; m, I, gy = 0.0)
    @named fixed_translation = Pl.FixedTranslation(; r = [10, 0])
    @named fixed = Pl.Fixed()
    @named revolute = Pl.Revolute()#constant_w = w)

    # sensors
    @named abs_v_sensor = Pl.AbsoluteVelocity(; resolve_in_frame)

    eqs = [
        connect(fixed.frame_b, revolute.frame_a),
        connect(revolute.frame_b, fixed_translation.frame_a),
        connect(fixed_translation.frame_b, body.frame_a),
        Pl.connect_sensor(body.frame_a, abs_v_sensor.frame_a)... # QUESTION: why?
        # connect(body.frame_a, abs_v_sensor.frame_a)
    ]

    @named model = ODESystem(eqs,
        t,
        [],
        [],
        systems = [
            body,
            fixed_translation,
            fixed,
            revolute,
            abs_v_sensor
        ])
    model = complete(model)
    @test_skip begin # Yingbo: BoundsError: attempt to access 137-element Vector{Vector{Int64}} at index [138]
        ssys = structural_simplify(IRSystem(model))
        prob = ODEProblem(ssys, [model.body.w => w], tspan)
        sol = solve(prob, Rodas5P(), initializealg=BrownFullBasicInit())

        # phi 
        @test sol[body.phi][end] ≈ tspan[end] * w
        @test all(sol[body.w] .≈ w)

        test_points = [i / w for i in 0:0.1:10]

        # instantaneous linear velocity
        v_signal(t) = -w^2 * sin.(w .* t)
        @test all(v_signal.(test_points) .≈ sol.(test_points; idxs = abs_v_sensor.v_x.u))

        # instantaneous linear acceleration
        a_signal(t) = -w^3 * cos.(w .* t)
        @test all(a_signal.(test_points) .≈ sol.(test_points; idxs = body.ax))
    end
end

@testset "Sensors (two free falling bodies)" begin
    m = 1
    I = 1
    resolve_in_frame = :world

    @named body1 = Pl.Body(; m, I)
    @named body2 = Pl.Body(; m, I)
    @named base = Pl.Fixed()

    @named abs_pos_sensor = Pl.AbsolutePosition(; resolve_in_frame)
    @named abs_v_sensor = Pl.AbsoluteVelocity(; resolve_in_frame)
    @named abs_a_sensor = Pl.AbsoluteAcceleration(; resolve_in_frame)
    @named rel_pos_sensor1 = Pl.RelativePosition(; resolve_in_frame)
    @named rel_pos_sensor2 = Pl.RelativePosition(; resolve_in_frame)
    @named rel_v_sensor1 = Pl.RelativeVelocity(; resolve_in_frame)
    @named rel_v_sensor2 = Pl.RelativeVelocity(; resolve_in_frame)
    @named rel_a_sensor1 = Pl.RelativeAcceleration(; resolve_in_frame)
    @named rel_a_sensor2 = Pl.RelativeAcceleration(; resolve_in_frame)

    connections = [
        Pl.connect_sensor(body1.frame_a, abs_pos_sensor.frame_a)...,
        Pl.connect_sensor(body1.frame_a, abs_v_sensor.frame_a)...,
        Pl.connect_sensor(body1.frame_a, abs_a_sensor.frame_a)...,
        Pl.connect_sensor(body1.frame_a, rel_pos_sensor1.frame_a)...,
        Pl.connect_sensor(base.frame_b, rel_pos_sensor1.frame_b)...,
        Pl.connect_sensor(body1.frame_a, rel_pos_sensor2.frame_a)...,
        Pl.connect_sensor(body2.frame_a, rel_pos_sensor2.frame_b)...,
        Pl.connect_sensor(base.frame_b, rel_v_sensor1.frame_a)...,
        Pl.connect_sensor(body1.frame_a, rel_v_sensor1.frame_b)...,
        Pl.connect_sensor(body1.frame_a, rel_v_sensor2.frame_a)...,
        Pl.connect_sensor(body2.frame_a, rel_v_sensor2.frame_b)...,
        Pl.connect_sensor(body1.frame_a, rel_a_sensor1.frame_a)...,
        Pl.connect_sensor(base.frame_b, rel_a_sensor1.frame_b)...,
        Pl.connect_sensor(body1.frame_a, rel_a_sensor2.frame_a)...,
        Pl.connect_sensor(body2.frame_a, rel_a_sensor2.frame_b)...
    ]

    # connections = [
    #     connect(body1.frame_a, abs_pos_sensor.frame_a),
    #     connect(body1.frame_a, abs_v_sensor.frame_a),
    #     connect(body1.frame_a, abs_a_sensor.frame_a),
    #     connect(body1.frame_a, rel_pos_sensor1.frame_a),
    #     connect(base.frame_b, rel_pos_sensor1.frame_b),
    #     connect(body1.frame_a, rel_pos_sensor2.frame_a),
    #     connect(body2.frame_a, rel_pos_sensor2.frame_b),
    #     connect(base.frame_b, rel_v_sensor1.frame_a),
    #     connect(body1.frame_a, rel_v_sensor1.frame_b),
    #     connect(body1.frame_a, rel_v_sensor2.frame_a),
    #     connect(body2.frame_a, rel_v_sensor2.frame_b),
    #     connect(body1.frame_a, rel_a_sensor1.frame_a),
    #     connect(base.frame_b, rel_a_sensor1.frame_b),
    #     connect(body1.frame_a, rel_a_sensor2.frame_a),
    #     connect(body2.frame_a, rel_a_sensor2.frame_b),
    # ]

    @named model = ODESystem(connections,
        t,
        [],
        [],
        systems = [
            body1,
            body2,
            base,
            abs_pos_sensor,
            abs_v_sensor,
            abs_a_sensor,
            rel_pos_sensor1,
            rel_pos_sensor2,
            rel_v_sensor1,
            rel_v_sensor2,
            rel_a_sensor1,
            rel_a_sensor2
        ])

    sys = structural_simplify((model)) # Yingbo: fails with JSCompiler
    unset_vars = setdiff(unknowns(sys), keys(ModelingToolkit.defaults(sys)))
    prob = ODEProblem(sys, unset_vars .=> 0.0, tspan)

    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)

    # the two bodyies falled the same distance, and so the absolute sensor attached to body1
    @test sol[abs_pos_sensor.y.u][end] ≈ sol[body1.r[2]][end] ≈ sol[body2.r[2]][end] ≈
          0.5 * g * tspan[end]^2

    # sensor1 is attached to body1, so the relative y-position between body1 and the base is
    # equal to the absolute y-position of body1
    @test sol[body1.r[2]][end] ≈ -sol[rel_pos_sensor1.rel_y.u][end]

    # the relative y-position between body1 and body2 is zero
    @test sol[rel_pos_sensor2.rel_y.u][end] == 0

    # no displacement in the x-direction
    @test sol[abs_pos_sensor.x.u][end] ≈ sol[body1.r[1]][end] ≈ sol[body2.r[1]][end]

    # velocity after t seconds v = g * t, so the relative y-velocity between body1 and the base is
    # equal to the absolute y-velocity of body1
    @test sol[abs_v_sensor.v_y.u][end] ≈ sol[rel_v_sensor1.rel_v_y.u][end] ≈ g * tspan[end]

    # the relative y-velocity between body1 and body2 is zero
    @test sol[rel_v_sensor2.rel_v_y.u][end] == 0

    # the body is under constant acclertation = g
    @test all(sol[abs_a_sensor.a_y.u] .≈ g)

    # the relative y-acceleration between body1 and the base is
    # equal to the absolute y-acceleration of body1
    @test sol[abs_a_sensor.a_y.u][end] ≈ -sol[rel_a_sensor1.rel_a_y.u][end]

    # the relative y-acceleration between body1 and body2 is zero
    @test sol[rel_a_sensor2.rel_a_y.u][end] == 0
end

@testset "Measure Demo" begin
    @named body = Pl.Body(; m = 1, I = 0.1)
    @named fixed_translation = Pl.FixedTranslation(;)
    @named fixed = Pl.Fixed()
    @named body1 = Pl.Body(; m = 0.4, I = 0.02)
    @named fixed_translation1 = Pl.FixedTranslation(; r = [0.4, 0])
    @named abs_pos_sensor = Pl.AbsolutePosition(; resolve_in_frame = :world)
    @named rel_pos_sensor = Pl.RelativePosition(; resolve_in_frame = :world)
    @named revolute1 = Pl.Revolute()
    @named abs_v_sensor = Pl.AbsoluteVelocity(; resolve_in_frame = :frame_a)
    @named rel_v_sensor = Pl.RelativeVelocity(; resolve_in_frame = :frame_b)
    @named abs_a_sensor = Pl.AbsoluteAcceleration(; resolve_in_frame = :world)
    @named rel_a_sensor = Pl.RelativeAcceleration(; resolve_in_frame = :frame_b)
    @named revolute2 = Pl.Revolute()

    connections = [
        connect(fixed_translation.frame_b, body.frame_a),
        connect(fixed_translation1.frame_b, body1.frame_a),
        connect(fixed.frame_b, revolute1.frame_a),
        connect(revolute1.frame_b, fixed_translation.frame_a),
        # connect(abs_a_sensor.frame_resolve, abs_a_sensor.frame_a),
        connect(revolute2.frame_b, fixed_translation1.frame_a),
        connect(revolute2.frame_a, fixed_translation.frame_b),
        # Pl.connect_sensor(fixed_translation.frame_b, rel_a_sensor.frame_a)...,
        # connect(fixed_translation.frame_b, rel_v_sensor.frame_a),
        # connect(fixed_translation.frame_b, rel_v_sensor.frame_a),
        # connect(rel_a_sensor.frame_b, body1.frame_a),
        # connect(rel_v_sensor.frame_b, body1.frame_a),
        # connect(rel_v_sensor.frame_b, body1.frame_a),
        Pl.connect_sensor(body1.frame_a, abs_a_sensor.frame_a)...        # Pl.connect_sensor(body1.frame, abs_v_sensor.frame_a)...,        # Pl.connect_sensor(body1.frame, abs_pos_sensor.frame_a)...,
    ]

    @named model = ODESystem(connections,
        t,
        [],
        [],
        systems = [
            fixed_translation,
            body,
            fixed,
            body1,
            fixed_translation1,
            revolute1,
            revolute2,
            abs_pos_sensor
        ])
    @test_skip begin # Yingbo: BoundsError again
        sys = structural_simplify(IRSystem(model))
        unset_vars = setdiff(unknowns(sys), keys(ModelingToolkit.defaults(sys)))
        prob = ODEProblem(sys, unset_vars .=> 0.0, (0, 5))
        sol = solve(prob, Rodas5P())
        @test SciMLBase.successful_retcode(sol)
    end
end

@testset "SpringDamper" begin
    @named spring_damper = Pl.SpringDamper(;
        s_relx0 = 0,
        d_y = 1,
        s_rely0 = 0,
        d_phi = 0,
        c_y = 5,
        c_x = 5,
        d_x = 1,
        c_phi = 0)
    @named body = Pl.Body(; I = 0.1, m = 0.5, r = [1,1], color=[0,1,0,1])
    @named fixed = Pl.Fixed()
    @named fixed_translation = Pl.FixedTranslation(; r = [-1, 0])

    connections = [
        connect(fixed.frame_b, fixed_translation.frame_a),
        connect(fixed_translation.frame_b, spring_damper.frame_a),
        connect(spring_damper.frame_b, body.frame_a)
    ]
    @named model = ODESystem(connections,
        t,
        [],
        [],
        systems = [
            spring_damper,
            body,
            fixed,
            fixed_translation
        ])
    sys = structural_simplify(IRSystem(model))
    unset_vars = setdiff(unknowns(sys), keys(ModelingToolkit.defaults(sys)))
    prob = ODEProblem(sys, unset_vars .=> 0.0, (0, 5))
    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)
end

@testset "Spring and damper demo" begin
    @named body = Pl.Body(; m = 0.5, I = 0.1)
    @named fixed = Pl.Fixed()
    @named spring = Pl.Spring(; c_y = 10, s_rely0 = -0.5, c_x = 1, c_phi = 1e5)
    @named damper = Pl.Damper(d = 1)
    @named prismatic = Pl.Prismatic(; r=[0, 1])

    connections = [
        connect(fixed.frame_b, spring.frame_a),
        connect(spring.frame_b, body.frame_a),
        connect(damper.frame_a, spring.frame_a),
        connect(damper.frame_b, spring.frame_b),
        connect(spring.frame_a, prismatic.frame_a),
        connect(prismatic.frame_b, spring.frame_b)
    ]

    @named model = ODESystem(connections,
        t,
        [],
        [],
        systems = [
            body,
            fixed,
            spring,
            damper,
            prismatic
        ])
    sys = structural_simplify(IRSystem(model)) # Yingbo: fails with JSCompiler
    unset_vars = setdiff(unknowns(sys), keys(ModelingToolkit.defaults(sys)))
    prob = ODEProblem(sys, unset_vars .=> 0.0, (0, 5), [])
    sol = solve(prob, Rodas5P(), initializealg=BrownFullBasicInit())
    @test SciMLBase.successful_retcode(sol)
end

##


@testset "SimpleWheel" begin
    @info "Testing SimpleWheel"
    gray = [0.1, 0.1, 0.1, 1]
    @mtkmodel TestWheel begin
        @components begin
            body = Pl.BodyShape(r = [1.0, 0.0], m=1, I=0.1, gy=0)
            revolute = Pl.Revolute()
            wheel1 = Pl.SimpleWheel(color=gray)
            wheel2 = Pl.SimpleWheel(color=gray, μ=.5)
            input = Blocks.Constant(k=1)
        end
        @equations begin
            connect(body.frame_a, revolute.frame_a)
            connect(revolute.frame_b, wheel1.frame_a)
            connect(input.output, wheel1.thrust)
            revolute.phi ~ deg2rad(50)sin(2pi*0.2*t)
            wheel2.thrust.u ~ 0

            connect(wheel2.frame_a, body.frame_b)
        end
    end
    @named model = TestWheel()
    model = complete(model)
    ssys = structural_simplify((model))
    defs = Dict(unknowns(ssys) .=> 0)
    prob = ODEProblem(ssys, defs, (0.0, 10.0))
    sol = solve(prob, Rodas5P(), initializealg = BrownFullBasicInit())
    @test SciMLBase.successful_retcode(sol)
    # Multibody.render(model, sol, show_axis=true, x=1, y=-1.8, z=5, lookat=[1,-1.8,0], traces=[model.wheel1.frame_a, model.wheel2.frame_a], filename="drifting.gif")
end


# import GLMakie, Multibody
# Multibody.render(model, sol, show_axis=true, x=1, y=1, z=5, traces=[model.wheel1.frame_a, model.wheel2.frame_a])


# plot(sol, idxs=[
#     model.revolute.phi,
#     model.revolute.frame_a.phi,
#     model.revolute.frame_b.phi,
#     model.wheel1.θ,
#     model.wheel1.frame.phi
# ])

##



import ModelingToolkitStandardLibrary.Mechanical.Rotational

@testset "SlipBasedWheel" begin
    @info "Testing SlipBasedWheel"

    @mtkmodel TestSlipBasedWheel begin
        @components begin
            slipBasedWheelJoint = Pl.SlipBasedWheelJoint(
                radius = 0.3,
                r = [1,0],
                mu_A = 0.8,
                mu_S = 0.4,
                N = 100,
                sAdhesion = 0.04,
                sSlide = 0.12,
                vAdhesion_min = 0.05,
                vSlide_min = 0.15,
                # w_roll = 10
            )
            prismatic = Pl.Prismatic(r = [0,1], s = 1, v = 0)
            revolute = Pl.Revolute(phi = 0, w = 0)
            fixed = Pl.Fixed()
            engineTorque = Rotational.ConstantTorque(tau_constant = 2)
            body = Pl.Body(m = 10, I = 1, gy=0, phi=0, w=0)
            inertia = Rotational.Inertia(J = 1, phi = 0, w = 0)
            constant = Blocks.Constant(k = 0)
        end
        @equations begin
            connect(prismatic.frame_a, revolute.frame_b)
            connect(revolute.frame_a, fixed.frame_b)
            connect(engineTorque.flange, inertia.flange_a)
            connect(body.frame_a, prismatic.frame_b)
            connect(slipBasedWheelJoint.frame_a, prismatic.frame_b)
            connect(slipBasedWheelJoint.flange_a, inertia.flange_b)
            connect(constant.output, slipBasedWheelJoint.dynamicLoad)
        end
    end

    @named model = TestSlipBasedWheel()
    model = complete(model)
    ssys = structural_simplify(IRSystem(model))
    display(unknowns(ssys))
    defs = ModelingToolkit.defaults(model)
    prob = ODEProblem(ssys, [
        model.inertia.w => 1e-10, # This is important, at zero velocity, the friction is ill-defined
        model.revolute.frame_b.phi => 0,
        model.body.w => 0,
        D(model.revolute.frame_b.phi) => 0,
        D(model.prismatic.r0[2]) => 0,
    ], (0.0, 20.0))
    sol = solve(prob, Rodas5Pr(autodiff=true)) # Since the friction model is not differentiable everywhere

    @test sol(15, idxs=[model.slipBasedWheelJoint.f_lat, model.slipBasedWheelJoint.f_long]) ≈ [80, -5] rtol=0.01
    @test sol(20, idxs=[model.slipBasedWheelJoint.f_lat, model.slipBasedWheelJoint.f_long]) ≈ [80, -4.95] rtol=0.01
    # plot(sol, idxs=[model.slipBasedWheelJoint.f_lat, model.slipBasedWheelJoint.f_long])
    # plot(sol, idxs=[model.revolute.w, model.prismatic.s])
end

##

@testset "TwoTrackModel" begin
    @info "Testing TwoTrackModel"
@mtkmodel TwoTrackWithDifferentialGear begin
    @components begin
        body = Pl.Body(m = 100, I = 1, gy = 0)
        body1 = Pl.Body(m = 300, I = 0.1, r = [1, 1], v = [0, 0], phi = 0, w = 0, gy = 0)
        body2 = Pl.Body(m = 100, I = 1, gy = 0)
        wheelJoint1 = Pl.SlipBasedWheelJoint(
            radius = 0.25,
            r = [0, 1],
            mu_A = 1,
            mu_S = 0.7,
            N = 1000,
            sAdhesion = 0.04,
            sSlide = 0.12,
            vAdhesion_min = 0.05,
            vSlide_min = 0.15,
            phi_roll = 0)
        wheelJoint2 = Pl.SlipBasedWheelJoint(
            radius = 0.25,
            r = [0, 1],
            mu_A = 1,
            mu_S = 0.7,
            N = 1500,
            sAdhesion = 0.04,
            sSlide = 0.12,
            vAdhesion_min = 0.05,
            vSlide_min = 0.15,
            phi_roll = 0)
        wheelJoint3 = Pl.SlipBasedWheelJoint(
            radius = 0.25,
            r = [0, 1],
            mu_A = 1,
            mu_S = 0.7,
            N = 1500,
            sAdhesion = 0.04,
            sSlide = 0.12,
            vAdhesion_min = 0.05,
            vSlide_min = 0.15,
            phi_roll = 0)
        wheelJoint4 = Pl.SlipBasedWheelJoint(
            radius = 0.25,
            r = [0, 1],
            mu_A = 1,
            mu_S = 0.7,
            N = 1000,
            sAdhesion = 0.04,
            sSlide = 0.12,
            vAdhesion_min = 0.05,
            vSlide_min = 0.15,
            phi_roll = 0)
        differentialGear = Pl.DifferentialGear()
        pulse = Blocks.Square(frequency = 1/2, offset = 0, start_time = 1, amplitude = -2)
        torque = Rotational.Torque()
        constantTorque1 = Rotational.ConstantTorque(tau_constant = 25)
        inertia = Rotational.Inertia(J = 1, phi = 0, w = 0)
        inertia1 = Rotational.Inertia(J = 1, phi = 0, w = 0)
        inertia2 = Rotational.Inertia(J = 1, phi = 0, w = 0)
        inertia3 = Rotational.Inertia(J = 1, phi = 0, w = 0)
        fixedTranslation1 = Pl.FixedTranslation(r = [0, 2])
        fixedTranslation2 = Pl.FixedTranslation(r = [0.75, 0])
        fixedTranslation3 = Pl.FixedTranslation(r = [-0.75, 0])
        fixedTranslation4 = Pl.FixedTranslation(r = [0.75, 0])
        fixedTranslation5 = Pl.FixedTranslation(r = [-0.75, 0])
        leftTrail = Pl.FixedTranslation(r = [0, -0.05])
        rightTrail = Pl.FixedTranslation(r = [0, -0.05])
        revolute = Pl.Revolute(axisflange=true)
        revolute2 = Pl.Revolute(axisflange=true, phi = -0.43633231299858, w = 0)
        dynamic_load = Blocks.Constant(k=0)
    end


    @equations begin
        connect(wheelJoint2.flange_a, inertia1.flange_b)
        connect(inertia.flange_b, wheelJoint1.flange_a)
        connect(fixedTranslation2.frame_b, fixedTranslation1.frame_a)
        connect(fixedTranslation2.frame_a, wheelJoint2.frame_a)
        connect(fixedTranslation3.frame_b, fixedTranslation1.frame_a)
        connect(wheelJoint3.frame_a, fixedTranslation3.frame_a)
        connect(inertia2.flange_b, wheelJoint3.flange_a)
        connect(body1.frame_a, fixedTranslation1.frame_a)
        connect(fixedTranslation1.frame_b, fixedTranslation4.frame_b)
        connect(fixedTranslation1.frame_b, fixedTranslation5.frame_b)
        connect(inertia3.flange_b, wheelJoint4.flange_a)
        connect(pulse.output, torque.tau)
        connect(differentialGear.flange_right, wheelJoint3.flange_a)
        connect(differentialGear.flange_left, wheelJoint2.flange_a)
        connect(constantTorque1.flange, differentialGear.flange_b)
        connect(body.frame_a, leftTrail.frame_b)
        connect(leftTrail.frame_b, wheelJoint1.frame_a)
        connect(body2.frame_a, rightTrail.frame_b)
        connect(wheelJoint4.frame_a, rightTrail.frame_b)
        connect(leftTrail.frame_a, revolute2.frame_a)
        connect(revolute2.frame_b, fixedTranslation4.frame_a)
        connect(torque.flange, revolute.flange_a, revolute2.flange_a)
        connect(revolute.frame_a, rightTrail.frame_a)
        connect(revolute.frame_b, fixedTranslation5.frame_a)
        connect(dynamic_load.output, wheelJoint1.dynamicLoad, wheelJoint2.dynamicLoad, wheelJoint3.dynamicLoad, wheelJoint4.dynamicLoad)
    end
end

@named model = TwoTrackWithDifferentialGear()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
defs = merge(
    Dict(unknowns(ssys) .=> 0),
    ModelingToolkit.defaults(model),
    Dict(model.body.w => 0),
)
prob = ODEProblem(ssys, defs, (0.0, 20.0))
sol = solve(prob, Rodas5P(autodiff=false))
@test SciMLBase.successful_retcode(sol)
# plot(sol)

end
