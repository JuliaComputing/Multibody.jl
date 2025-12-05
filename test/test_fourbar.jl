using Test
using Multibody
using ModelingToolkit
import ModelingToolkitStandardLibrary.Mechanical.Rotational
using OrdinaryDiffEq
using LinearAlgebra
# using JuliaSimCompiler

t = Multibody.t
D = Differential(t)
world = Multibody.world

systems = @named begin
    j1 = Revolute(n = [1, 0, 0], w0 = 5.235987755982989, state_priority=10.0, radius=0.1f0) # Increase state priority to ensure that this joint coordinate is chosen as state variable
    j2 = Prismatic(n = [1, 0, 0], s0 = -0.2)
    b1 = BodyShape(r = [0, 0.5, 0.1], radius=0.03)
    b2 = BodyShape(r = [0, 0.2, 0], radius=0.03)
    b3 = BodyShape(r = [-1, 0.3, 0.1], radius=0.03)
    rev = Revolute(n = [0, 1, 0])
    rev1 = Revolute()
    j3 = Revolute(n = [1, 0, 0])
    j4 = Revolute(n = [0, 1, 0])
    j5 = Revolute(n = [0, 0, 1], iscut=true)
    b0 = FixedTranslation(r = [1.2, 0, 0], radius=0)
end

connections = [connect(j2.frame_b, b2.frame_a)
               connect(j1.frame_b, b1.frame_a)
               connect(rev.frame_a, b2.frame_b)
               connect(rev.frame_b, rev1.frame_a)
               connect(rev1.frame_b, b3.frame_a)
               connect(world.frame_b, j1.frame_a)
               connect(b1.frame_b, j3.frame_a)
               connect(j3.frame_b, j4.frame_a)
               connect(j4.frame_b, j5.frame_a)
               connect(j5.frame_b, b3.frame_b)
               connect(b0.frame_a, world.frame_b)
               connect(b0.frame_b, j2.frame_a)
               ]
@named fourbar2 = System(connections, t, systems = [world; systems])
fourbar2 = complete(fourbar2)
ssys = structural_simplify(multibody(fourbar2))

prob = ODEProblem(ssys, [], (0.0, 1.4399)) # The end time is chosen to make the animation below appear to loop forever

sol = solve(prob, FBDF(autodiff=true));
@test SciMLBase.successful_retcode(sol)
# plot(sol, idxs=[j2.s]) # Plot the joint coordinate of the prismatic joint (green in the animation below)


systems = @named begin
    j1 = Revolute(n = [1, 0, 0], w0 = 5.235987755983, state_priority=12.0, radius=0.1f0) # Increase state priority to ensure that this joint coordinate is chosen as state variable
    j2 = Prismatic(n = [1, 0, 0], s0 = -0.2)
    b1 = BodyShape(r = [0, 0.5, 0.1], radius=0.03)
    b2 = BodyShape(r = [0, 0.2, 0], radius=0.03)
    b3 = FixedTranslation(r = [1.2, 0, 0], radius=0)
    universalSpherical = UniversalSpherical(n1_a = [0, 1, 0], rRod_ia = [-1, 0.3, 0.1])
end

connections = [connect(j2.frame_b, b2.frame_a)
               connect(j1.frame_b, b1.frame_a)
               connect(j1.frame_a, world.frame_b)
               connect(b1.frame_b, universalSpherical.frame_b)
               connect(universalSpherical.frame_a, b2.frame_b)
               connect(b3.frame_a, world.frame_b)
               connect(b3.frame_b, j2.frame_a)
]

@named model = System(connections, t, systems = [world; systems])
model = complete(model)
ssys = structural_simplify(multibody(model))
prob = ODEProblem(ssys, [], (0.0, 1.4399)) # The end time is chosen to make the animation below appear to loop forever
sol2 = solve(prob, FBDF(autodiff=true)) # 3.9x faster than above
# plot(sol2, idxs=[j2.s]) # Plot the joint coordinate of the prismatic joint (green in the animation below)


using LinearAlgebra
@test norm(sol(0:0.1:1, idxs=j2.s) - sol2(0:0.1:1, idxs=j2.s)) < 0.2