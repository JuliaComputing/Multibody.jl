using Multibody, ModelingToolkit, Test, OrdinaryDiffEq
@component function FallingBody(; name)
    systems = @named begin
        my_world = World(g = 1, n = [0, 1, 0])
        body = Body(isroot=true)
    end

    pars = @parameters begin
    end

    vars = @variables begin
    end

    equations = [
    ]

    return System(equations, t; name, systems)
end

@named model = FallingBody()
model = complete(model)
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, [], (0, 1))
sol = solve(prob, Rodas5P())

tv = 0:0.1:1
@test iszero(sol(tv, idxs=model.body.r_0[1]))
@test sol(tv, idxs=model.body.r_0[2]) ≈ tv.^2 ./ 2 atol=1e-6
@test iszero(sol(tv, idxs=model.body.r_0[3]))


# Change g
prob = ODEProblem(ssys, [model.my_world.g .=> 2], (0, 1))
sol = solve(prob, Rodas5P())
@test iszero(sol(tv, idxs=model.body.r_0[1]))
@test sol(tv, idxs=model.body.r_0[2]) ≈ 2*tv.^2 ./ 2 atol=1e-6
@test iszero(sol(tv, idxs=model.body.r_0[3]))

# Change n
prob = ODEProblem(ssys, collect(model.my_world.n) .=> [1, 0, 0], (0, 1))
sol = solve(prob, Rodas5P())
@test sol(tv, idxs=model.body.r_0[1]) ≈ tv.^2 ./ 2 atol=1e-6
@test iszero(sol(tv, idxs=model.body.r_0[2]))
@test iszero(sol(tv, idxs=model.body.r_0[3]))


## World in more than one place
# This should result in too many equations
@component function FallingBodyOuter(; name)
    systems = @named begin
        inner_model = FallingBody()
        my_world_outer = World(g = 2, n = [0, 1, 0])
        body = Body(isroot=true)
    end

    pars = @parameters begin
    end

    vars = @variables begin
    end

    equations = [
    ]

    return System(equations, t; name, systems)
end

@named model = FallingBodyOuter()
model = complete(model)
@test_throws ModelingToolkit.ExtraEquationsSystemException structural_simplify(IRSystem(model))

