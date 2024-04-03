
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
        FixedRotation( # NOTE: This prevents the nodes from moving w.r.t. each other, the movement is handled by the interface nodes
            r = translations[i],
            n = rot_axes[i],
            angle = rot_angles[i];
            name = Symbol("transform$i")
        )
    end

    @named common_reference_frame = Frame()

    # QUESTION: Should each mode also be associated with a Body attached to the interface node, + one body for the common frame? Can M be divided up into pieces for each body? If M is block-diagonal with at most 3x3 blocks, this would be easy.

    # Generalized reaction forces at interface nodes
    forces_b = [[collect(node.frame_b.f); collect(node.frame_b.tau)] for node in interface_nodes]
    forces_a = [[collect(node.frame_a.f); collect(node.frame_a.tau)] for node in interface_nodes]
    forces_b = reduce(vcat, forces_b)
    forces_a = reduce(vcat, forces_a)

    @variables begin
        x(t)[1:ntot]=0, [description = "CB Displacement"]
        v(t)[1:ntot]=0, [description = "CB Velocity"]
        a(t)[1:ntot]=0, [description = "CB Acceleration"]
        f(t)[1:ntot]=0, [description = "CB Reaction forces"]
    end
    x, v, a, f = collect.((x, v, a, f))
    
    
    
    #=
    u = Φ*q # + H.O.T.
    r = r_R + c + u

    x = q in the modelica FlexibleBodies paper
    there is a transformation Φ that relates elastic displacement u to the modal coordinates q according to u = Φq + H.O.T. This transformation depends on the point on the flexible object one is considering, i.e., Φ = Φ(c) ∈ R(3 × N)

    Currently, we do not distinguish between u and q and have more or less u = q[1:3] = x[1:3]
    =#

    if D === nothing
        eqs = M*a       + K*x .~ f # + gravity forces + quadratic velocity forces
    else
        eqs = M*a + D*v + K*x .~ f # NOTE: see eq (7) in https:#www.bausch-gall.de/Session1d1.pdf This models the scenario when the motion of the reference frame is constrained to be zero
    end
    append!(eqs, Multibody.D.(x) .~ v)
    append!(eqs, Multibody.D.(v) .~ a)
    append!(eqs, forces_b .~ f)
    # append!(eqs, forces_a .~ f) # The free motion is not completely free, it carries force through the interface
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
- Handle the rigid body mode and gravity
- Add support for internal dynamic modes
=#

function FlexibleCoupling(; D=nothing, K, T, name)

    ntot = LinearAlgebra.checksquare(K)
    ntot == 6 || error("K must be a 6×6 matrix")
    @parameters begin
        K[1:ntot, 1:ntot] = K, [description = "CB Stiffness matrix"]
    end
    K = collect(K)

    if D !== nothing
        size(D) == (ntot, ntot) || error("D must be $ntot x $ntot")
        @parameters begin
            D[1:ntot, 1:ntot] = D, [description = "CB Damping matrix"]
        end
        D = collect(D)
    end

    translation = T2t(T)
    rot_angles = rotation_angle(T2R(T))
    @named fm = FreeMotion(r_rel_a = translation, phi = rot_angles, isroot=true)
    @unpack frame_a, frame_b = fm
    r_rel_a, phi = fm.r_rel_a, fm.phi
    @variables begin
        x(t)[1:ntot]=0, [state_priority = 10, description = "CB Displacement"]
        v(t)[1:ntot]=0, [state_priority = 10, description = "CB Velocity"]
        a(t)[1:ntot]=0, [state_priority = 10, description = "CB Acceleration"]
        f(t)[1:ntot]=0, [state_priority = 10, description = "CB Reaction forces"]
    end

    @named begin
        frame_a = Frame()
        frame_b = Frame()
    end

    @named force = BasicForce(resolveInFrame = :frame_a)
    @named torque = BasicTorque(resolveInFrame = :frame_a)

    # Ra = ori(frame_a)

    # fa = frame_a.f |> collect
    # fb = frame_b.f |> collect
    # taua = frame_a.tau |> collect
    # taub = frame_b.tau |> collect
    # r = translation + x[1:3]

    # Rdeflect = axesRotations([1,2,3], x[4:6], D.(x[4:x]), name = :R_ar)

    # eqs = [
    #     D === nothing ? (K*x .~ f) : (D*v + K*x .~ f)
    #     Multibody.D.(x) .~ v
    #     Multibody.D.(v) .~ a
    #     fa .+ fb .~ f[1:3]
    #     taua .+ taub + cross(r, fb) .~ f[4:6]
    #     frame_b.r_0 .~ frame_a.r_0 + resolve1(ori(frame_a), r)
    #     ori(frame_b) ~ Rdeflect*ori(frame_a)
    # ]

    
    eqs = [
        (D === nothing ? (K*x) : (D*v + K*x)) .~ f
        Multibody.D.(x) .~ v
        Multibody.D.(v) .~ a

        force.force.u .~ f[1:3]
        torque.torque.u .~ f[4:6]
        connect(frame_a, fm.frame_a, force.frame_a, torque.frame_a)
        connect(frame_b, fm.frame_b, force.frame_b, torque.frame_b)

        # frame_a.f + frame_b.f .~ force.force.u
        # frame_a.tau + frame_b.tau .~ torque.torque.u
        force.r_0 .~ r_rel_a
        # torque.r_0 .~ r_rel_a

        x[1:3] .~ r_rel_a .- translation
        x[4:6] .~ phi .- rot_angles

    ]

    ODESystem(eqs, t; name, systems=[force, torque, frame_a, frame_b, fm])

end


function FlexibleLink(; b1, b2, D=nothing, K, T, name)

    ntot = LinearAlgebra.checksquare(K)
    ntot == 6 || error("K must be a 6×6 matrix")
    @parameters begin
        K[1:ntot, 1:ntot] = K, [description = "CB Stiffness matrix"]
    end
    K = collect(K)

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
        FixedRotation( # NOTE: This prevents the nodes from moving w.r.t. each other, the movement is handled by the interface nodes
            r = translations[i],
            n = rot_axes[i],
            angle = rot_angles[i];
            name = Symbol("transform$i")
        )
    end

    @named common_reference_frame = Frame()

    # QUESTION: Should each mode also be associated with a Body attached to the interface node, + one body for the common frame? Can M be divided up into pieces for each body? If M is block-diagonal with at most 3x3 blocks, this would be easy.

    # Generalized reaction forces at interface nodes
    forces_b = [[collect(node.frame_b.f); collect(node.frame_b.tau)] for node in interface_nodes]
    forces_a = [[collect(node.frame_a.f); collect(node.frame_a.tau)] for node in interface_nodes]
    forces_b = reduce(vcat, forces_b)
    forces_a = reduce(vcat, forces_a)

    @variables begin
        x(t)[1:ntot]=0, [description = "CB Displacement"]
        v(t)[1:ntot]=0, [description = "CB Velocity"]
        a(t)[1:ntot]=0, [description = "CB Acceleration"]
        f(t)[1:ntot]=0, [description = "CB Reaction forces"]
    end
    x, v, a, f = collect.((x, v, a, f))
    
    
    
    #=
    u = Φ*q # + H.O.T.
    r = r_R + c + u
    
    x = q in the modelica FlexibleBodies paper
    there is a transformation Φ that relates elastic displacement u to the modal coordinates q according to u = Φq + H.O.T. This transformation depends on the point on the flexible object one is considering, i.e., Φ = Φ(c) ∈ R(3 × N)

    Currently, we do not distinguish between u and q and have more or less u = q[1:3] = x[1:3]
    =#

    if D === nothing
        eqs =       K*x .~ f # + gravity forces + quadratic velocity forces
    else
        eqs = D*v + K*x .~ f # NOTE: see eq (7) in https:#www.bausch-gall.de/Session1d1.pdf This models the scenario when the motion of the reference frame is constrained to be zero
    end
    append!(eqs, Multibody.D.(x) .~ v)
    append!(eqs, Multibody.D.(v) .~ a)
    append!(eqs, forces_b .~ f)
    # append!(eqs, forces_a .~ f) # The free motion is not completely free, it carries force through the interface
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











# # From https://www.tandfonline.com/doi/epdf/10.1080/13873954.2013.807433?needAccess=false
# # Object-oriented modelling of general flexible multibody systemsGianni Ferrettia*, Alberto Levaaand Bruno Scaglioni

# dq = der(q) # elastic coordinates velocities
# ddq = der(dq) # elastic coordinates accelerations
# [data.inv1*I(3),mD̄̃c',transpose(data.Ct)]*[aa_a - g_0;za_a;ddq] = h_w_r+h_e_r # (2)
# [mD̄̃c,Jbar,Cr']*[aa_a - g_0;za_a;ddq]  =  h_w_θ+h_e_θ # (2)
# [data.Ct,Cr,data.Me]*[aa_a - g_0;za_a;ddq] = h_w_f+h_e_f -(data.Ke*q)−((d/100)*(alpha*data.Me+beta*data.Ke)*dq) # (2)
# mD̄c = data.inv2+sum(data.inv3[:,i]*q[i] for i in 1:M) # (4)
# mD̄̃c = skew(mD̄c) # (4)
# Cr' = data.inv4+sum(data.inv5[i,:,:]*q[i] for i in 1:M) # (6)
# Jbar = data.inv7−sum((transpose(data.inv8[i,:,:])+data.inv8[i,:,:])*q[i] for i in 1:M)−sum(sum(data.inv9[i,j,:,:]*q[j] for j in i:M)*q[i] for i in 1:M) # (8)
# h_w_r = (-cross(wa_a,cross(wa_a, mD̄c))−2*cross(wa_a, (transpose(data.Ct)*dq))) # (9)
# h_w_θ = (-cross(wa_a, Jbar*wa_a) - (der(Jbar)*wa_a)−cross(wa_a, (Cr'*dq))) # (10)
# h_w_f = (hwfint'*wa_a−2*[wa_a[1]*(data.inv11[3,2,:,:]-data.inv11[2,3,:,:])+wa_a[2]*(data.inv11[1,3,:,:]-data.inv11[3,1,:,:])+wa_a[3]*(data.inv11[2,1,:,:]-data.inv11[1,2,:,:])]*dq) # (12)
# for i in 1:3
#     for j in 1:3
#         D[i,j,:] = data.inv10[i,j,:]+data.inv11[i,j,:,:]*q # (13)
#     end 
# end 
# sumd = sum(D[i,i,:] for i in 1:3) # (13)
# for i in 1:3
#     hwfint[i,:] = wa_a[i]*sumd-wa_a*D[i,:,:] # (13)
# end
# h_e_r = (sum(fb_a[i,:] for i in 1:data.Nc)) # (14)
# h_e_θ = (sum((taub_a[i,:]+cross((data.r0b[i,:]+data.Sb[i,:,:]*q),fb_a[i,:])) for i in 1:data.Nc)) # (15)
# h_e_f = sum((transpose(data.Sb[i,:,:])*fb_a[i,:]+transpose(data.SbCap[i,:,:])*taub_a[i,:]) for i in 1:data.Nc) # (16)

