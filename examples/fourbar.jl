using Multibody, OrdinaryDiffEq, JuliaSimCompiler
using ModelingToolkit
using Multibody: connect_loop
using Test
using Plots
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
systems = @named begin
    j1 = Revolute(n = [1, 0, 0], w0 = 5.235987755982989, isroot = true)
    j2 = Prismatic(n = [1, 0, 0], s0 = -0.2, isroot=false)
    b1 = BodyShape(r = [0, 0.5, 0.1])
    b2 = BodyShape(r = [0, 0.2, 0])
    b3 = BodyShape(r = [-1, 0.3, 0.1])
    rev = Revolute(n = [0, 1, 0], isroot = true, iscut=true)
    rev1 = Revolute()
    j3 = Revolute(n = [1, 0, 0], isroot = true)
    j4 = Revolute(n = [0, 1, 0], isroot = true)
    j5 = Revolute(n = [0, 0, 1], isroot = true)
    b0 = FixedTranslation(r = [1.2, 0, 0])
end

connections = [connect(j2.frame_b, b2.frame_a)

            #    Multibody.connect_loop2(j1.frame_b, b1.frame_a)
               connect(j1.frame_b, b1.frame_a)
            
            #    Multibody.connect_loop(rev.frame_a, b2.frame_b)
               connect(rev.frame_a, b2.frame_b)

               connect(rev.frame_b, rev1.frame_a)
               connect(rev1.frame_b, b3.frame_a)
               connect(world.frame_b, j1.frame_a)
               connect(b1.frame_b, j3.frame_a)m.mass
               connect(j3.frame_b, j4.frame_a)
               connect(j4.frame_b, j5.frame_a)
               connect(j5.frame_b, b3.frame_b)
               connect(b0.frame_a, world.frame_b)
               connect(b0.frame_b, j2.frame_a)
               ]
@named fourbar = ODESystem(connections, t, systems = [world; systems])

# m = structural_simplify(fourbar)
m = structural_simplify(IRSystem(fourbar))

prob = ODEProblem(m, [], (0.0, 5.0))

du = zero(prob.u0)
prob.f(du, prob.u0, prob.p, 0.0) 

sol = solve(prob, Rodas4())#, u0 = prob.u0 .+ 0.01 .* randn.())

# using SeeToDee, NonlinearSolve
# function dynamics(x,u,p,t)
#     dx = similar(x)
#     prob.f(dx,x,p,t)
# end

# Ts = 0.001
# nx = 2
# na = 7
# nu = 0
# x_inds = findall(!iszero, prob.f.mass_matrix.diag)
# a_inds = findall(iszero, prob.f.mass_matrix.diag)
# discrete_dynamics = SeeToDee.SimpleColloc2(dynamics, Ts, x_inds, a_inds, nu; n=3, solver=NonlinearSolve.NewtonRaphson())

# X = [prob.u0]
# i = 0
# for i in 1:1000
#     push!(X, discrete_dynamics(X[end], [], prob.p, i*Ts))
# end



# ==============================================================================
## Trivial 4 bar
# ==============================================================================


# First test the structure without the loop closed, this makes a quadruple pendulum


systems = @named begin
    j1 = Revolute(isroot = true)
    j2 = Revolute(isroot = true)
    j3 = Revolute(isroot = true)
    j4 = Revolute(isroot = true)
    b1 = BodyShape(m=1, r = [1.0, 0, 0], r_cm = 0.5*[1.0, 0, 0])
    b2 = BodyShape(m=3, r = [1.0, 0, 0], r_cm = 0.5*[1.0, 0, 0], radius=0.3)
    b3 = BodyShape(m=1, r = [-1.0, 0, 0], r_cm = 0.5*[-1.0, 0, 0])
    b4 = BodyShape(m=1, r = [1.0, 0, 0], r_cm = 0.5*[1.0, 0, 0])
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

@time "Quadpend" begin
    # m = structural_simplify(fourbar)
    m = structural_simplify(IRSystem(fourbar))
    @test length(states(m)) == 8
    prob = ODEProblem(m, [], (0.0, 10))

    @time sol = solve(prob, Rodas4())#, u0 = prob.u0 .+ 0.01 .* randn.())
end
isinteractive() && plot(sol)

# render(fourbar, sol; framerate=60)

## Now close the loop
# When we do, we can only select one joint as root. We must also replace one joint with a RevolutePlanarLoopConstraint 
systems = @named begin
    j1 = Revolute(isroot = true)
    j2 = Revolute(isroot = false)
    j3 = Revolute(isroot = false)
    j4 = RevolutePlanarLoopConstraint()
    b1 = BodyShape(r = [1.0, 0, 0])
    b2 = BodyShape(r = [1.0, 0, 0])
    b3 = BodyShape(r = [-1.0, 0, 0])
    b4 = BodyShape(r = [1.0, 0, 0])
end

connections = [
    connect(world.frame_b, j1.frame_a)
    # Multibody.connect_loop(world.frame_b, j1.frame_a)
    
    connect(j1.frame_b, b1.frame_a)
    connect(b1.frame_b, j2.frame_a)
    connect(j2.frame_b, b2.frame_a)
    connect(b2.frame_b, j3.frame_a)
    connect(j3.frame_b, b3.frame_a)
    connect(b3.frame_b, j4.frame_a)
    
    connect(j4.frame_b, b4.frame_a)
    # Multibody.connect_loop(j4.frame_b, b4.frame_a)
    
    connect(b4.frame_b, j1.frame_a)
    # Multibody.connect_loop(b4.frame_b, j1.frame_a)
    
]
@named fourbar = ODESystem(connections, t, systems = [world; systems])

##

# m = structural_simplify((fourbar), allow_parameter=false)
m = structural_simplify(IRSystem(fourbar)) # It does simplify

@test_broken length(states(m)) == 2

prob = ODEProblem(m, [], (0.0, 5.0))

# Try the generated dynamics
du = zero(prob.u0)
@time prob.f.f(du, prob.u0, prob.p, 0) # Yingbo: Singular


sol = solve(prob, Rodas4())#, u0 = prob.u0 .+ 0.01 .* randn.())
isinteractive() && plot(sol, vars = [j1.phi, j2.phi, j3.phi])