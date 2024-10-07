using Multibody, ModelingToolkit, JuliaSimCompiler, Test
@mtkmodel FallingBody begin
    @components begin
        my_world = World(g = 1, n = [0, 1, 0])
        body = Body(isroot=true)
    end
end

@named model = FallingBody()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, [], (0, 1))
sol = solve(prob, Tsit5())

tv = 0:0.1:1
@test iszero(sol(tv, idxs=model.body.r_0[1]))
@test sol(tv, idxs=model.body.r_0[2]) ≈ tv.^2 ./ 2 atol=1e-6
@test iszero(sol(tv, idxs=model.body.r_0[3]))


# Change g
prob = ODEProblem(ssys, [model.my_world.g .=> 2], (0, 1))
sol = solve(prob, Tsit5())
@test iszero(sol(tv, idxs=model.body.r_0[1]))
@test sol(tv, idxs=model.body.r_0[2]) ≈ 2*tv.^2 ./ 2 atol=1e-6
@test iszero(sol(tv, idxs=model.body.r_0[3]))

# Change n
prob = ODEProblem(ssys, collect(model.my_world.n) .=> [1, 0, 0], (0, 1))
sol = solve(prob, Tsit5())
@test sol(tv, idxs=model.body.r_0[1]) ≈ tv.^2 ./ 2 atol=1e-6
@test iszero(sol(tv, idxs=model.body.r_0[2]))
@test iszero(sol(tv, idxs=model.body.r_0[3]))


## World in more than one place
# This should result in too many equations
@mtkmodel FallingBodyOuter begin
    @components begin
        inner_model = FallingBody()
        my_world_outer = World(g = 2, n = [0, 1, 0])
        body = Body(isroot=true)
    end
end

@named model = FallingBodyOuter()
model = complete(model)
@test_throws ModelingToolkit.ExtraEquationsSystemException structural_simplify(IRSystem(model))

