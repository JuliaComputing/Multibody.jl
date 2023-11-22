
"""
    CraigBampton(; M, D = nothing, K, T)

A Craig-Bampton approximation of a flexible body.

# Arguments
- `M`: Mass matrix, `n x n`
- `D`: Damping matrix (optional), `n x n`
- `K`: Stiffness matrix, `n x n`
- `T`: A vector of transformations (SE(3)) from common frame of reference to each interface node. 

# Internals
The transformations from the common reference frame to the interface nodes are modeled using [`FixedRotation`](@ref). Each interface node is located at `frame_b` of the transformation model, e.g., the first interface node frame is `cb.interface1.frame_b`.
"""
function CraigBampton(; M, D=nothing, K, T, name)

    ntot = LinearAlgebra.checksquare(M)
    ntot % 6 == 0 || error("M must be square and have a side length that is a multiple of 6")
    n = ntot ÷ 6
    size(K) == (ntot, ntot) || error("K must be $ntot x $ntot")
    length(T) == n || error("T must be a vector of length $n")
    @parameters begin
        M[1:ntot, 1:ntot] = M, [description = "CB Mass matrix"]
        K[1:ntot, 1:ntot] = K, [description = "CB Stiffness matrix"]
    end
    M,K = collect.((M, K))

    if D !== nothing
        size(D) == (ntot, ntot) || error("D must be $ntot x $ntot")
        @parameters begin
            D[1:ntot, 1:ntot] = D, [description = "CB Damping matrix"]
        end
        D = collect(D)
    end

    translations = map(T2t, T)
    rot_axes = map(rotation_axis ∘ T2R, T)
    rot_angles = map(rotation_angle ∘ T2R, T)

    interface_nodes = map(1:n) do i
        FreeMotion(
            # r_rel_a = translations[i],
            name = Symbol("interface$i"),
            # isroot = false,
        )
    end

    node_transforms = map(1:n) do i
        FixedRotation( # QUESTION: This prevents the nodes from moving w.r.t. each other, is this desired?
            r = translations[i],
            n = rot_axes[i],
            angle = rot_angles[i];
            name = Symbol("transform$i")
        )
    end

    @named common_reference_frame = Frame()

    # Generalized reaction forces at interface nodes
    forces = [[collect(node.frame_b.f); collect(node.frame_b.tau)] for node in interface_nodes]
    forces = reduce(vcat, forces)

    @variables begin
        x(t)[1:ntot]=0, [description = "CB Displacement"]
        v(t)[1:ntot]=0, [description = "CB Velocity"]
        a(t)[1:ntot]=0, [description = "CB Acceleration"]
        f(t)[1:ntot]=0, [description = "CB Reaction forces"]
    end
    x, v, a, f = collect.((x, v, a, f))

    if D === nothing
        eqs = M*a       + K*x .~ f
    else
        eqs = M*a + D*v + K*x .~ f
    end
    append!(eqs, Multibody.D.(x) .~ v)
    append!(eqs, Multibody.D.(v) .~ a)
    append!(eqs, forces .~ f)
    for i = 1:n
        push!(eqs, connect(common_reference_frame, node_transforms[i].frame_a))
        push!(eqs, connect(node_transforms[i].frame_b, interface_nodes[i].frame_a))
    end


    inds = 1:3
    for i = 1:n
        # assert equality between deflections and free motion coordinates
        append!(eqs, collect(x[inds] .~ interface_nodes[i].r_rel_a))
        inds = inds .+ 3
        append!(eqs, collect(x[inds] .~ interface_nodes[i].phi))
        inds = inds .+ 3
    end

    ODESystem(eqs, t; systems = [common_reference_frame; interface_nodes; node_transforms], name)

end


#=
TODO
- Handle the rigid body mode
- Add support for internal dynamic modes
=#