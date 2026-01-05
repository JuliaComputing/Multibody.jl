using ModelingToolkit
using Multibody
using Test
# using JuliaSimCompiler
using OrdinaryDiffEq
doplot() = false
world = Multibody.world

t = Multibody.t
D = Differential(t)

# ==============================================================================
## The simplest of them all, two WorldForces attached to the same body. The forces add to zero, resulting in no movement
# ==============================================================================
@component function TestWorldForce1(; name)
    pars = @parameters begin
        f[1:3]
    end

    systems = @named begin
        world = World()
        forcea = WorldForce(resolve_frame=:frame_b)
        forceb = WorldForce(resolve_frame=:frame_b)
        b = Body(m=1, state=true, isroot=true, quat=false)
    end

    vars = @variables begin
    end

    f = collect(f)

    equations = Equation[
        connect(b.frame_a, forcea.frame_b, forceb.frame_b)
        forcea.force.u ~ f
        forceb.force.u ~ -f
    ]

    return System(equations, t; name, systems)
end

@named testwf = TestWorldForce1()
testwf = complete(testwf)
ssys = multibody(testwf)
prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [0,1,0]], (0, 1))
sol = solve(prob, Tsit5())
# plot(sol)
@test sol(1, idxs=testwf.b.r_0) ≈ [0, 0.0, 0.0] atol=1e-3

# ==============================================================================
## A single force
# ==============================================================================
# Tests here pass for any combination of quat and neg_w
@component function TestWorldForce2(; name)
    pars = @parameters begin
        f[1:3]
    end

    systems = @named begin
        world = World()
        force = WorldForce()
        body = Body(m=1, state=true, isroot=true, quat=false)
    end

    vars = @variables begin
    end

    f = collect(f)

    equations = Equation[
        # connect(world.frame_b, body.frame_a)
        connect(force.frame_b, body.frame_a)
        force.force.u ~ f
    ]

    return System(equations, t; name, systems)
end

@named testwf = TestWorldForce2()
testwf = complete(testwf)
ssys = multibody(testwf)
prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [1,0,0]], (0, 1))
sol = solve(prob, Tsit5())
# plot(sol)
@test sol(1, idxs=testwf.body.r_0) ≈ [0.5*1^2*1, 0, 0] atol=1e-3


prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [0,1,0]], (0, 1))
sol = solve(prob, Tsit5())
# plot(sol)
@test sol(1, idxs=testwf.body.r_0) ≈ [0, 0.5*1^2*1, 0] atol=1e-3


# ==============================================================================
## Two forces on each side of a BodyShape, when pointing towards each other (or away from each other) they should cancel out
# When pointing in different directions, they should cause the body to spin around its center of mass
# ==============================================================================
@component function TestWorldForce3(; name)
    pars = @parameters begin
        f[1:3]
    end

    systems = @named begin
        world = World()
        forcea = WorldForce(resolve_frame=:world)
        forceb = WorldForce(resolve_frame=:world)
        body = BodyShape(r=[1,0,0], state=true, isroot=true, quat=false)
    end

    vars = @variables begin
    end

    f = collect(f)

    equations = Equation[
        # connect(world.frame_b, body.frame_a)
        connect(forcea.frame_b, body.frame_a)
        connect(forceb.frame_b, body.frame_b)
        forcea.force.u ~ f
        forceb.force.u ~ -f
    ]

    return System(equations, t; name, systems)
end

@named testwf = TestWorldForce3()
testwf = complete(testwf)
ssys = multibody(testwf)
prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [1,0,0]], (0, 1))
sol = solve(prob, Tsit5())
# plot(sol)
@test sol(1, idxs=testwf.body.r_0) ≈ [0, 0, 0] atol=1e-3


prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [0,1,0]], (0, 1))
sol = solve(prob, Tsit5(), reltol=1e-8)
# plot(sol)
@test sol(1, idxs=testwf.body.frame_a.r_0) ≈ [0.572800369240885, 0.4946715021692289, 0.0] atol=1e-3 # Spinning around center of mass

# ==============================================================================
## Same as above but with forces resolved in frame_b instead of world frame
# ==============================================================================
@component function TestWorldForce4(; name)
    pars = @parameters begin
        f[1:3]
    end

    systems = @named begin
        world = World()
        forcea = WorldForce(resolve_frame=:frame_b)
        forceb = WorldForce(resolve_frame=:frame_b)
        body = BodyShape(r=[1,0,0], state=true, isroot=true, quat=false)
    end

    vars = @variables begin
    end

    f = collect(f)

    equations = Equation[
        connect(forcea.frame_b, body.frame_a)
        connect(forceb.frame_b, body.frame_b)
        forcea.force.u ~ f
        forceb.force.u ~ -f
    ]

    return System(equations, t; name, systems)
end

@named testwf = TestWorldForce4()
testwf = complete(testwf)
ssys = multibody(testwf)
prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [1,0,0]], (0, 1))
sol = solve(prob, Tsit5())
# plot(sol)
@test sol(1, idxs=testwf.body.r_0) ≈ [0, 0, 0] atol=1e-3


prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [0,1,0]], (0, 1))
sol = solve(prob, Tsit5(), reltol=1e-8)
# plot(sol)
@test sol(1, idxs=testwf.body.frame_a.r_0) ≈ [0.9419246090353689, -0.23388592078659548, 0.0] atol=1e-3

# ==============================================================================
## Same as above but with BodyCylinder instead
# ==============================================================================
@component function TestWorldForce5(; name)
    pars = @parameters begin
        f[1:3]
    end

    systems = @named begin
        world = World()
        forcea = WorldForce(resolve_frame=:frame_b)
        forceb = WorldForce(resolve_frame=:frame_b)
        body = BodyCylinder(r=[1,0,0], state=true, isroot=true, quat=false, density=1, diameter=0.1)
    end

    vars = @variables begin
    end

    f = collect(f)

    equations = Equation[
        connect(forcea.frame_b, body.frame_a)
        connect(forceb.frame_b, body.frame_b)
        forcea.force.u ~ f
        forceb.force.u ~ -f
    ]

    return System(equations, t; name, systems)
end

@named testwf = TestWorldForce5()
testwf = complete(testwf)
ssys = multibody(testwf)
prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [1,0,0]], (0, 1))
sol = solve(prob, Tsit5())
# plot(sol)
@test sol(1, idxs=testwf.body.r_0) ≈ [0, 0, 0] atol=1e-3


prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [0,1,0]], (0, 0.1))
sol = solve(prob, Tsit5(), reltol=1e-8)
# plot(sol)
@test sol(0.1, idxs=testwf.body.frame_a.r_0) ≈ [0.36637355877861, 0.4818152481205165, 0.0] atol=1e-3

# ==============================================================================
## First create a body and then attach Cylinder to this. The body has (almost) zero mass, so the result should be identical to the test above
# ==============================================================================
using LinearAlgebra
@component function TestWorldForce6(; name)
    pars = @parameters begin
        f[1:3]
    end

    systems = @named begin
        world = World()
        forcea = WorldForce(resolve_frame=:frame_b, radius=0.15, scale=0.2)
        forceb = WorldForce(resolve_frame=:frame_b, radius=0.15, scale=0.2)
        b0 = Body(m=1e-32, I_11=1e-32, I_22=1e-32, I_33=1e-32, state_priority=0, radius=0.14, color=[1,0,0,0.5])
        body = BodyCylinder(r=[1,0,0], density=1, diameter=0.1, state=true, isroot=true, quat=false, color=[0,0,1,0.5])
    end

    vars = @variables begin
    end

    f = collect(f)

    equations = Equation[
        connect(forcea.frame_b, body.frame_a, b0.frame_a)
        connect(forceb.frame_b, body.frame_b)
        forcea.force.u ~ f
        forceb.force.u ~ -f
    ]

    return System(equations, t; name, systems)
end

@named testwf = TestWorldForce6()
testwf = complete(testwf)
ssys = multibody(testwf)
prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [1,0,0]], (0, 1))
sol = solve(prob, Tsit5())
# plot(sol)
@test sol(1, idxs=testwf.body.r_0) ≈ [0, 0, 0] atol=1e-3
@test sol(1, idxs=testwf.body.body.w_a[3]-testwf.b0.w_a[3]) ≈ 0 atol=1e-3 # These should be identical

prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [0,1,0]], (0, 1))
sol = solve(prob, Rodas4(), reltol=1e-8)
# plot(sol)
@test sol(0.1, idxs=testwf.body.frame_a.r_0) ≈ [0.36595999301484056, 0.48169308140123934, 0.0] atol=1e-3 # Identical to the test above without extra body b0

@test sol(1, idxs=testwf.body.body.w_a[3]-testwf.b0.w_a[3]) ≈ 0 atol=1e-3 # These should be identical
# plot(sol, idxs=[testwf.body.body.w_a[3], testwf.b0.w_a[3]])

# ==============================================================================
## Same as above but attach b0 to body.frame_a, this should yields identical results, but currently does not
# ==============================================================================
@component function TestWorldForce7(; name)
    pars = @parameters begin
        f[1:3]
    end

    systems = @named begin
        world = World()
        forcea = WorldForce(resolve_frame=:frame_b, radius=0.15, scale=0.2)
        forceb = WorldForce(resolve_frame=:frame_b, radius=0.15, scale=0.2)
        b0 = Body(m=1e-32, I_11=1e-32, I_22=1e-32, I_33=1e-32, state_priority=0, radius=0.14, color=[1,0,0,0.5])
        body = BodyCylinder(r=[1,0,0], density=1, diameter=0.1, color=[0,0,1,0.5], state=true, isroot=true, quat=false)
    end

    vars = @variables begin
    end

    f = collect(f)

    equations = Equation[
        connect(forcea.frame_b, body.frame_a, b0.frame_a)
        connect(forceb.frame_b, body.frame_b)
        forcea.force.u ~ f
        forceb.force.u ~ -f
    ]

    return System(equations, t; name, systems)
end

@named testwf = TestWorldForce7()
testwf = complete(testwf)
ssys = multibody(testwf)
prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [1,0,0]], (0, 1))
sol = solve(prob, Tsit5())
# plot(sol)
@test sol(1, idxs=testwf.body.r_0) ≈ [0, 0, 0] atol=1e-3


prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [0,1,0]], (0, 0.1))
sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8)
# plot(sol)
@test sol(0.1, idxs=testwf.body.frame_a.r_0) ≈ [0.36595999301484056, 0.48169308140123934, 0.0] atol=1e-3



# ==============================================================================
## Simplification of the model above into two bodies
# This might be the simplest test that demonstrates the problem with one force set to 0
# ==============================================================================

# =============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
@component function TestWorldForce8(; name)
    pars = @parameters begin
        f[1:3]
    end

    systems = @named begin
        world = World()
        forcea = WorldForce(resolve_frame=:frame_b)
        forceb = WorldForce(resolve_frame=:frame_b)
        b0 = Body(m=1, state_priority=0, radius=0.1, color=[0,0,1,0.2])
        b1 = Body(m=1, state=true, isroot=true, quat=false, radius=0.05, color=[1,0,0,1])
    end

    vars = @variables begin
    end

    f = collect(f)

    equations = Equation[
        connect(forceb.frame_b, forcea.frame_b, b1.frame_a, b0.frame_a)
        forcea.force.u ~ f
        forceb.force.u ~ -f
    ]

    return System(equations, t; name, systems)
end

@named testwf = TestWorldForce8()
testwf = complete(testwf)
ssys = multibody(testwf)
prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [1,0,0]], (0, 1))
sol = solve(prob, Tsit5())
# plot(sol)
@test sol(1, idxs=testwf.b1.r_0) ≈ [0, 0.0, 0.0] atol=1e-3

reshape(sol(1, idxs = [testwf.forceb.frame_b.f; testwf.forcea.frame_b.f; testwf.b1.frame_a.f; testwf.b0.frame_a.f;]), 3, :)

@test iszero(sol(1, idxs=(testwf.b0.frame_a.f)[1] + (testwf.b1.frame_a.f)[1] + (testwf.forcea.frame_b.f)[1] + (testwf.forceb.frame_b.f)[1]))

## Same with quat=true

@component function TestWorldForce9(; name)
    pars = @parameters begin
        f[1:3]
    end

    systems = @named begin
        world = World()
        forcea = WorldForce(resolve_frame=:frame_b)
        forceb = WorldForce(resolve_frame=:frame_b)
        b0 = Body(m=1, state_priority=0, radius=0.1, color=[0,0,1,0.2])
        body = Body(m=1, state=true, isroot=true, quat=true, radius=0.05, color=[1,0,0,1])
    end

    vars = @variables begin
    end

    f = collect(f)

    equations = Equation[
        connect(forcea.frame_b, body.frame_a, b0.frame_a)
        connect(forceb.frame_b, body.frame_a)
        forcea.force.u ~ f
        forceb.force.u ~ -f
    ]

    return System(equations, t; name, systems)
end

@named testwf = TestWorldForce9()
testwf = complete(testwf)
ssys = multibody(testwf)
prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [1,0,0]], (0, 1))
sol = solve(prob, Tsit5())
@test sol(1, idxs=testwf.body.r_0) ≈ [0, 0.0, 0.0] atol=1e-3

## Removing the body b0 from the model

@component function TestWorldForce10(; name)
    pars = @parameters begin
        f[1:3]
    end

    systems = @named begin
        world = World()
        forcea = WorldForce(resolve_frame=:frame_b)
        forceb = WorldForce(resolve_frame=:frame_b)
        body = Body(m=1, state=true, isroot=true, quat=false, radius=0.05, color=[1,0,0,1])
    end

    vars = @variables begin
    end

    f = collect(f)

    equations = Equation[
        connect(forcea.frame_b, body.frame_a)
        connect(forceb.frame_b, body.frame_a)
        forcea.force.u ~ f
        forceb.force.u ~ -f
    ]

    return System(equations, t; name, systems)
end

@named testwf = TestWorldForce10()
testwf = complete(testwf)
ssys = multibody(testwf)
prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [1,0,0]], (0, 1))
sol = solve(prob, Tsit5())
# plot(sol)
@test sol(1, idxs=testwf.body.r_0) ≈ [0, 0.0, 0.0] atol=1e-3

# ==============================================================================
## Almost same as above but with FixedTranslation instead of r_cm and force application at the end of the rod
# ==============================================================================
@component function TestWorldForce11(; name)
    pars = @parameters begin
        f[1:3]
    end

    systems = @named begin
        world = World()
        forcea = WorldForce(resolve_frame=:frame_b)
        forceb = WorldForce(resolve_frame=:frame_b)
        b0 = Body(m=1, state_priority=0)
        body = Body(m=1, state=true, isroot=true, quat=false)
        tr = FixedTranslation(r=[1,0,0])
    end

    vars = @variables begin
    end

    f = collect(f)

    equations = Equation[
        connect(forcea.frame_b, b0.frame_a, tr.frame_a)
        connect(forceb.frame_b, body.frame_a, tr.frame_b)
        forcea.force.u ~ f
        forceb.force.u ~ -f
    ]

    return System(equations, t; name, systems)
end

@named testwf = TestWorldForce11()
testwf = complete(testwf)
ssys = multibody(testwf)
prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [1,0,0]], (0, 1))
sol = solve(prob, Tsit5())
# plot(sol)
@test sol(1, idxs=testwf.body.r_0) ≈ [0, 0.0, 0.0] atol=1e-3


prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [0,10,0]], (0, 1))
sol = solve(prob, Tsit5())
# plot(sol)
@test !iszero(sol(1, idxs=testwf.body.frame_a.r_0))
@test sol(1, idxs=testwf.body.frame_a.r_0) ≈ [-0.9300324062366484, 0.25508128572301375, 0.0] atol=1e-3

sol(0, idxs=testwf.forcea.frame_b.f + testwf.b0.frame_a.f + testwf.tr.frame_a.f)
sol(0, idxs=testwf.forceb.frame_b.f + testwf.body.frame_a.f + testwf.tr.frame_b.f)

sol(1, idxs=testwf.body.frame_a.f)
sol(1, idxs=testwf.tr.frame_b.f)

sol(1, idxs=testwf.body.w_a)
sol(1, idxs=testwf.b0.w_a)

# ==============================================================================
## Same as above but with state in FreeMotion instead
# ==============================================================================
# This works with both quat and Euler
@component function TestWorldForce12(; name)
    pars = @parameters begin
        f[1:3]
    end

    systems = @named begin
        world = World()
        freemotion = FreeMotion(state=true, isroot=true, quat=false)
        forcea = WorldForce(resolve_frame=:frame_b)
        forceb = WorldForce(resolve_frame=:frame_b)
        b0 = Body(m=1)
        body = Body(m=1)
        tr = FixedTranslation(r=[1,0,0])
    end

    vars = @variables begin
    end

    f = collect(f)

    equations = Equation[
        connect(world.frame_b, freemotion.frame_a)
        connect(forcea.frame_b, b0.frame_a, tr.frame_a, freemotion.frame_b)
        connect(forceb.frame_b, body.frame_a, tr.frame_b)
        forcea.force.u ~ f
        forceb.force.u ~ -f
    ]

    return System(equations, t; name, systems)
end

@named testwf = TestWorldForce12()
testwf = complete(testwf)
ssys = multibody(testwf)
prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [1,0,0]], (0, 1))
sol = solve(prob, Tsit5())
# plot(sol)
@test sol(1, idxs=testwf.body.r_0) ≈ [0, 0.0, 0.0] atol=1e-3


prob = ODEProblem(ssys, [testwf.world.g => 0; collect(testwf.f) .=> [0,10,0]], (0, 1))
sol = solve(prob, Tsit5())
# plot(sol)
@test !iszero(sol(1, idxs=testwf.body.frame_a.r_0))
@test sol(1, idxs=testwf.body.frame_a.r_0) ≈ [-0.9300324062366484, 0.25508128572301375, 0.0] atol=1e-3

@test iszero(sol(0, idxs=testwf.forcea.frame_b.f + testwf.b0.frame_a.f + testwf.tr.frame_a.f))
@test iszero(sol(0, idxs=testwf.forceb.frame_b.f + testwf.body.frame_a.f + testwf.tr.frame_b.f))



##