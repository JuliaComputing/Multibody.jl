using Multibody, OrdinaryDiffEq, JuliaSimCompiler
using ModelingToolkit
using Multibody: connect_loop
using Test
t = Multibody.t

world = Multibody.world


W(args...; kwargs...) = Multibody.world

# @mtkmodel FourBar begin
#     @variables begin
#         j1_phi(t), [description = "Angle of revolute joint j1", output = true]
#         j2_s(t), [description = "Distance of prismatic joint j2", output = true]
#         j1_w(t), [description = "Axis speed of revolute joint j1", output = true]
#         j2_v(t), [description = "Axis velocity of prismatic joint j2", output = true]
#     end
#     @components begin
#         world = W()
#         j1 = Revolute(
#             n=[1,0,0],
#             # stateSelect=StateSelect.always,
#             w0=5.235987755982989,
#             # fixed=true
#             isroot = false,
#             )
#         j2 = Prismatic(n=[1,0,0], s0 = -0.2, isroot=false)
#         b1 = BodyShape(r=[0,0.5,0.1], isroot=false)
#         b2 = BodyShape(r=[0,0.2,0], isroot=false)
#         b3 = BodyShape(r=[-1,0.3,0.1])
#         rev = Revolute(n=[0,1,0], isroot=false)
#         rev1 = Revolute(isroot=true)
#         j3 = Revolute(n=[1,0,0])
#         j4 = Revolute(n=[0,1,0])
#         j5 = Revolute(n=[0,0,1])
#         b0 = FixedTranslation(r=[1.2,0,0])
#     end

#     @equations begin
#         connect(j2.frame_b, b2.frame_a)
#         connect(j1.frame_b, b1.frame_a)
#         # connect(rev.frame_a, b2.frame_b)
#         connect_loop(rev.frame_a, b2.frame_b)
#         connect(rev.frame_b, rev1.frame_a)
#         connect(rev1.frame_b, b3.frame_a) 
#         connect(world.frame_b, j1.frame_a)
#         connect(b1.frame_b, j3.frame_a)
#         connect(j3.frame_b, j4.frame_a)
#         connect(j4.frame_b, j5.frame_a)
#         connect(j5.frame_b, b3.frame_b)
#         connect(b0.frame_a, world.frame_b)
#         connect(b0.frame_b, j2.frame_a)
#         j1_phi ~ j1.phi
#         j2_s ~ j2.s
#         j1_w ~ j1.w
#         j2_v ~ j2.v
#     end
# end

# @named fourbar = FourBar()
@named begin
    j1 = Revolute(n = [1, 0, 0], w0 = 5.235987755982989, isroot = false)
    j2 = Prismatic(n = [1, 0, 0], s0 = -0.2, isroot=true)
    b1 = BodyShape(r = [0, 0.5, 0.1])
    b2 = BodyShape(r = [0, 0.2, 0])
    b3 = BodyShape(r = [-1, 0.3, 0.1])
    rev = Revolute(n = [0, 1, 0])
    rev1 = Revolute()
    j3 = Revolute(n = [1, 0, 0])
    j4 = Revolute(n = [0, 1, 0])
    j5 = Revolute(n = [0, 0, 1])
    b0 = FixedTranslation(r = [1.2, 0, 0])
end

connections = [connect(j2.frame_b, b2.frame_a)
               connect(j1.frame_b, b1.frame_a)
               Multibody.connect_loop(rev.frame_a, b2.frame_b)
            #    connect(rev.frame_a, b2.frame_b)
               connect(rev.frame_b, rev1.frame_a)
               connect(rev1.frame_b, b3.frame_a)
               connect(world.frame_b, j1.frame_a)
               connect(b1.frame_b, j3.frame_a)
               connect(j3.frame_b, j4.frame_a)
               connect(j4.frame_b, j5.frame_a)
               connect(j5.frame_b, b3.frame_b)
               connect(b0.frame_a, world.frame_b)
               connect(b0.frame_b, j2.frame_a)]
@named fourbar = ODESystem(connections, t,
                         systems = [
                             world,
                             j1,
                             j2,
                             b1,
                             b2,
                             b3,
                             rev,
                             rev1,
                             j3,
                             j4,
                             j5,
                             b0,
                         ])

# m = structural_simplify(fourbar)
m = structural_simplify(IRSystem(fourbar))

prob = ODEProblem(m, [], (0.0, 5.0))

sol = solve(prob, Rodas4(), u0 = prob.u0 .+ 0.01 .* randn.())



# ==============================================================================
## Trivial 4 bar
# ==============================================================================


# First test the structure without the loop closed, this makes a quadruple pendulum

systems = @named begin
    j1 = Revolute(isroot = true)
    j2 = Revolute(isroot = true)
    j3 = Revolute(isroot = true)
    j4 = Revolute(isroot = true)
    b1 = BodyShape(r = [1.0, 0, 0])
    b2 = BodyShape(r = [1.0, 0, 0])
    b3 = BodyShape(r = [1.0, 0, 0])
    b4 = BodyShape(r = [1.0, 0, 0])
end

connections = [
    connect(world.frame_b, j1.frame_a)
    connect(j1.frame_b, b1.frame_a)
    connect(b1.frame_b, j2.frame_a)
    connect(j2.frame_b, b2.frame_a)
    connect(b2.frame_b, j3.frame_a)
    connect(j3.frame_b, b3.frame_a)
    connect(b3.frame_b, j4.frame_a)
    connect(j4.frame_b, b4.frame_a)
]
@named fourbar = ODESystem(connections, t, systems = [world; systems])

# m = structural_simplify(fourbar)
m = structural_simplify(IRSystem(fourbar))
@test length(states(m)) == 8
prob = ODEProblem(m, [], (0.0, 10))


## Now close the loop
# When we do, we can only select one joint as root. The loop closure must not use a regular connect statement, instead, we use connect_loop
# For unknown reason, we must also change the connection to the world to use connect_loop for the system to balance
systems = @named begin
    j1 = Revolute(isroot = true)
    j2 = Revolute(isroot = false)
    j3 = Revolute(isroot = false)
    j4 = Revolute(isroot = false)
    b1 = BodyShape(r = [1.0, 0, 0])
    b2 = BodyShape(r = [1.0, 0, 0])
    b3 = BodyShape(r = [1.0, 0, 0])
    b4 = BodyShape(r = [1.0, 0, 0])
end

connections = [
    # connect(world.frame_b, j1.frame_a)
    connect(j1.frame_b, b1.frame_a)
    connect(b1.frame_b, j2.frame_a)
    connect(j2.frame_b, b2.frame_a)
    connect(b2.frame_b, j3.frame_a)
    connect(j3.frame_b, b3.frame_a)
    connect(b3.frame_b, j4.frame_a)
    connect(j4.frame_b, b4.frame_a)
    # connect(b4.frame_b, j1.frame_a)
    Multibody.connect_loop(b4.frame_b, j1.frame_a)
    Multibody.connect_loop(world.frame_b, j1.frame_a)
]
@named fourbar = ODESystem(connections, t, systems = [world; systems])

m = structural_simplify(IRSystem(fourbar))

@test_broken length(states(m)) == 2

prob = ODEProblem(m, [], (0.0, 5.0))

# Try the generated dynamics
du = zero(prob.u0)
prob.f.f(du, prob.u0, prob.p, 0)


sol = solve(prob, Rodas4())#, u0 = prob.u0 .+ 0.01 .* randn.())
isinteractive() && plot(sol, vars = [j1.phi, j2.phi, j3.phi, j4.phi])