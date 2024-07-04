function build_fun(prob::ODEProblem, inputs, outputs)
    ssys = prob.f.sys
    fkine = JuliaSimCompiler.build_explicit_observed_function(ssys, outputs)
    x = nameof.(unknowns(ssys))
    nx = length(x)

    # ins = [JuliaSimCompiler.ADT.IRElement((i.val.f)) for i in inputs]
    ins = [(nameof(i.val.f)) for i in inputs]
    
    inds = [findfirst(==(i), x) for i in ins]
    any(isnothing, inds) && error("Could not locate all inputs in the state vector, failed to find $(inputs[inds .== nothing])")

    let fkine = fkine, bufferf64 = zeros(nx)
        function kine(q::AbstractVector{Float64}, p=prob.p, t=0)
            bufferf64[inds] .= q
            fkine(bufferf64, p, t)
        end
        function kine(q, p=prob.p, t=0)
            T = eltype(q)
            buffer = zeros(T, nx)
            buffer[inds] .= q
            fkine(buffer, p, t)
        end
    end
end