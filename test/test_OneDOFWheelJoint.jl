using ModelingToolkit, OrdinaryDiffEq, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkitStandardLibrary.Blocks
import ModelingToolkitStandardLibrary.Mechanical.Rotational
import Multibody
using Multibody: multibody
import Multibody.PlanarMechanics as Pl
# using JuliaSimCompiler

tspan = (0.0, 3.0)
g = -9.80665

@testset "OneDOFWheelJoint" begin
    @info "Testing OneDOFWheelJoint"

    # Wheel spinning in origin
    @component function SimpleTest(; name)
        systems = @named begin
            body = Pl.Body(m = 0.1, I = 0.1, phi=0, w=30, gy=-9.82)
            wheelJoint = Pl.OneDOFWheelJoint(
                x = 0,
                v = 0,
                radius = 1,
                mu_A = 0.95,
                mu_S = 0.7,
                sAdhesion = 0.04,
                sSlide = 0.12,
                vAdhesion_min = 0.05,
                vSlide_min = 0.15,
            )
        end

        vars = @variables begin
        end

        eqs = [
            # connect(wheelJoint.flange_a, inertia.flange_b)
            connect(wheelJoint.frame_a, body.frame_a)
        ]

        System(eqs, t, vars, []; systems, name)
    end

    @named model = SimpleTest()
    ssys = multibody(model)
    prob = ODEProblem(ssys, [ssys.wheelJoint.frame_a.render => true, ssys.wheelJoint.frame_a.length => 1.3], (0.0, 3))
    sol = solve(prob, Rodas5P())

    @test sol(3, idxs=ssys.body.v[1] + ssys.body.w) ≈ 0 atol=1e-3 # The rotational speed eventually match linear speed (radius = 1)
    @test sol(3, idxs=ssys.body.v[1]) < 0 # Goes in the left direction



    ## Planar segway
    @component function PlanarSegway(; name)
        pars = @parameters begin
            r_cm = 2
            radius = 1
        end
        systems = @named begin
            body = Pl.Body(m = 0.1, I = 0.0001, phi=0, w=0.0, radius=0.02, gy=-9.82)
            body2 = Pl.Body(m = 10, I = 0.001, radius=0.02, gy=-9.82)
            translation_cm = Pl.FixedTranslation(r = [r_cm, 0], radius=0.01)
            wheelJoint = Pl.OneDOFWheelJoint(
                radius = radius,
                x = 0,
                v = 0,
                mu_A = 1,
                mu_S = 0.7,
                sAdhesion = 0.04,
                sSlide = 0.12,
                vAdhesion_min = 0.05,
                vSlide_min = 0.15,
            )
        end
        eqs = [
            connect(wheelJoint.frame_a, translation_cm.frame_a, body2.frame_a)
            connect(translation_cm.frame_b, body.frame_a)
        ]
        System(eqs, t, [], pars; systems, name)
    end

    @named model = PlanarSegway()
    ssys = multibody(model)
    prob = ODEProblem(ssys,
        [
            ssys.wheelJoint.frame_a.render => true, ssys.wheelJoint.frame_a.length => 3, ssys.wheelJoint.frame_a.radius => 0.002,
            ssys.body.frame_a.render => true, ssys.body.frame_a.length => 3, ssys.body.frame_a.radius => 0.002,
            ssys.body2.frame_a.render => true, ssys.body2.frame_a.length => 3, ssys.body2.frame_a.radius => 0.002,
        ], (0.0, 50); 
        missing_guess_value = MissingGuessValue.Constant(0.0),
    )

    sol = solve(prob, Rodas5P(autodiff=true))

    @test sol(0.1, idxs=ssys.body.frame_a.phi - ssys.wheelJoint.frame_a.phi) == 0 # These two should be the same
    @test sol(0.1, idxs=ssys.body.w + ssys.wheelJoint.w_roll) == 0 # These two should be the same but opposite sign
    @test sol(0.1, idxs=ssys.body.w) < 0
    @test sol(0.1, idxs=ssys.wheelJoint.v) > 0 # Wheel should be pulled towards the falling body (positive x)
    @test sol(0.1, idxs=ssys.wheelJoint.frame_a.x) > 0 # Wheel should be pulled towards the falling body (positive x)
    
end