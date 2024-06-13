function centraldiff(v::AbstractMatrix)
    dv = Base.diff(v, dims=1)/2
    a1 = [dv[[1],:];dv]
    a2 = [dv;dv[[end],:]]
    a = a1+a2
end

function centraldiff(v::AbstractVector)
    dv = Base.diff(v)/2
    a1 = [dv[1];dv]
    a2 = [dv;dv[end]]
    a = a1+a2
end



"""
    q,qd,qdd,t_end = point_to_point(time; q0 = 0.0, q1 = 1.0, t0 = 0, qd_max = 1, qdd_max = 1)

Generate a minimum-time point-to-point trajectory with specified start and endpoints, not exceeding specified speed and acceleration limits.

The trajectory produced by this function will typically exhibit piecewise constant accleration, piecewise linear velocity and piecewise quadratic position curves.

If a vector of `time` points is provided, the function returns matrices `q,qd,qdd` of size `(length(time), n_dims)`. If a scalar `time` point is provided, the function returns `q,qd,qdd` as vectors with the specified dimension (same dimension as `q0`).
`t_end` is the time at which the trajectory will reach the specified end position.

# Arguments:
- `time`: A scalar or a vector of time points.
- `q0`: Initial coordinate, may be a scalar or a vector.
- `q1`: End coordinate
- `t0`: Tiem at which the motion starts. If `time` contains time points before `t0`, the trajectory will stand still at `q0` until `time` reaches `t0`.
- `qd_max`: Maximum allowed speed.
- `qdd_max`: Maximum allowed acceleration.

See also [`KinematicPTP`](@ref) and [`traj5`](@ref).
"""
function point_to_point(time; q0 = 0.0, q1 = 1.0, t0 = 0, qd_max = 1, qdd_max = 1)
    nout = max(length(q0), length(q1), length(qd_max), length(qdd_max))
    p_q0 = q0
    p_q1 = q1
    p_qd_max = qd_max
    p_qdd_max = qdd_max
    p_deltaq = p_q1 .- p_q0

    T = Base.promote_eltype(q0, q1, t0, qd_max, qdd_max)

    aux1 = p_deltaq ./ p_qd_max
    aux2 = p_deltaq ./ p_qdd_max

    sd_max_inv = maximum(abs, aux1)
    sdd_max_inv = maximum(abs, aux2)

    if sd_max_inv <= eps() || sdd_max_inv <= eps()
        if time isa Number
            return q0, zero(q0), zero(q0), float(time)
        else
            return zeros(T, nout) .+ q0', zeros(T, length(time), nout), zeros(T, length(time), nout), float(first(time))
        end
    end

    sd_max = 1 / sd_max_inv
    sdd_max = 1 / sdd_max_inv
    Ta1 = sqrt(1 / sdd_max)
    Ta2 = sd_max / sdd_max
    noWphase = Ta2 >= Ta1
    Tv = if noWphase
        Ta1
    else
        1 / sd_max
    end
    Te = if noWphase
        Ta1 + Ta1
    else
        Tv + Ta2
    end
    Ta1s = Ta1 + t0
    Ta2s = Ta2 + t0
    Tvs = Tv + t0
    Tes = Te + t0
    sd_max2 = sdd_max * Ta1
    s1 = sdd_max * (noWphase ?
                    Ta1 * Ta1 : Ta2 * Ta2) / 2
    s2 = s1 + (noWphase ? sd_max2 * (Te - Ta1) - (sdd_max / 2) * (Te - Ta1)^2 :
            sd_max * (Tv - Ta2))

    s3 = s2 + sd_max * (Te - Tv) - (sdd_max / 2) * (Te - Tv) * (Te - Tv)

    if !(time isa Number)
        N = length(time)
        Q = zeros(T, N, nout)
        Qd = zeros(T, N, nout)
        Qdd = zeros(T, N, nout)
    end

    for (i, t) in enumerate(time)
        sd = sdd = zero(T)

        if t < t0
            s = zero(T)
        elseif noWphase
            if t < Ta1s
                s = (sdd_max / 2) * (t - t0)^2
                sd = sdd_max*(t - t0)
                sdd = sdd_max
            elseif t < Tes
                s = s1 + sd_max2 * (t - Ta1s) -
                    (sdd_max / 2) * (t - Ta1s)^2
                sd = sd_max2 + (Ta1s - t)*sdd_max
                sdd = -sdd_max
            else
                s = s2
            end
        elseif t < Ta2s
            s = (sdd_max / 2) * (t - t0)^2
            sd = sdd_max*(t - t0)
            sdd = sdd_max
        elseif t < Tvs
            s = s1 + sd_max * (t - Ta2s)
            sd = sd_max
        elseif t < Tes
            s = s2 + sd_max * (t - Tvs) - (sdd_max / 2) * (t - Tvs)^2
            sd = sd_max + (Tvs - t)*sdd_max
            sdd = -sdd_max
        else
            s = s3
        end
        
        if time isa Number
            qdd = p_deltaq * sdd
            qd = p_deltaq * sd
            t1 = Tes
            q = @. p_q0 + p_deltaq * s
            return q, qd, qdd, t1
        else
            @. Qdd[i, :] = p_deltaq * sdd
            @. Qd[i, :] = p_deltaq * sd
            @. Q[i, :] = p_q0 + p_deltaq * s
        end
    end
    @assert !(time isa Number)
    return Q, Qd, Qdd, t1

end

