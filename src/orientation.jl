using Rotations

_norm(x) = sqrt(sum(abs2(x) for x in x)) # Workaround for buggy symbolic arrays
_normalize(x) = x ./ _norm(x)

const R3{T} = RotMatrix{3, T}

abstract type Orientation end

"""
    RotationMatrix

A struct representing a 3D orientation as a rotation matrix.

If `ODESystem` is called on a `RotationMatrix` object `o`, symbolic variables for `o.R` and `o.w` are created and the value of `o.R` is used as the default value for the symbolic `R`.

# Fields:
- `R::R3`: The rotation 3×3 matrix ∈ SO(3)
- `w`: The angular velocity vector
"""
struct RotationMatrix <: Orientation
    R::R3
    w::Any
end

RotationMatrix(R::AbstractMatrix, w) = RotationMatrix(R3(R), w)

RotationMatrix() = RotationMatrix(R3(1.0I(3)), zeros(3))


"""
    NumRotationMatrix(; R = collect(1.0 * I(3)), w = zeros(3), name, varw = false)

Create a new [`RotationMatrix`](@ref) struct with symbolic elements. `R,w` determine default values.

The primary difference between `NumRotationMatrix` and `RotationMatrix` is that the `NumRotationMatrix` constructor is used in the constructor of a [`Frame`](@ref) in order to introduce the frame variables, whereas `RorationMatrix` (the struct) only wraps existing variables.

- `varw`: If true, `w` is a variable, otherwise it is derived from the derivative of `R` as `w = get_w(R)`.

Never call this function directly from a component constructor, instead call `f = Frame(); R = ori(f)` and add `f` to the subsystems.
"""
function NumRotationMatrix(; R = collect(1.0I(3)), w = zeros(3), name=:R, varw = false, state_priority=nothing)
    # The reason for not calling this directly is that all R vaiables have to have the same name since they are treated as connector variables (otherwise a connection error is thrown). A component with more than one rotation matrix will thus have two different R variables that overwrite each other
    R = at_variables_t(:R, 1:3, 1:3; default = R, state_priority) #[description="Orientation rotation matrix ∈ SO(3)"]
    # @variables w(t)[1:3]=w [description="angular velocity"]
    # R = collect(R)
    # R = ModelingToolkit.renamespace.(name, R) .|> Num
    if varw
        w = at_variables_t(:w, 1:3, default = w)
    else
        w = get_w(R)
    end
    R, w = collect.((R, w))
    RotationMatrix(R, w)
end

nullrotation() = RotationMatrix()

function ModelingToolkit.ODESystem(RM::RotationMatrix; name)
    # @variables R(t)[1:3, 1:3]=Matrix(RM) [description="Orientation rotation matrix ∈ SO(3)"]
    # @variables w(t)[1:3]=w [description="angular velocity"]
    # R,w = collect.((R,w))
    R = at_variables_t(:R, 1:3, 1:3)
    w = at_variables_t(:w, 1:3)

    defaults = Dict(R .=> RM)
    ODESystem(Equation[], t, [vec(R); w], []; name, defaults)
end

Base.:*(R1::RotationMatrix, x::AbstractArray) = R1.R * x
Base.:*(x::AbstractArray, R2::RotationMatrix) = x * R2.R
function Base.:*(R1::RotationMatrix, R2::RotationMatrix)
    RotationMatrix(R1.R.mat * R2.R.mat, R1 * R2.w + collect(R1.w))
end
LinearAlgebra.adjoint(R::RotationMatrix) = RotationMatrix(R.R', -R.w)

function (D::Differential)(RM::RotationMatrix)
    # https://build.openmodelica.org/Documentation/Modelica.Mechanics.MultiBody.Frames.Orientation.html
    R = RM.R
    DR = D.(RM.R)
    Dw = [R[3, :]'DR[2, :], -R[3, :]'DR[1, :], R[2, :]'DR[1, :]]
    RotationMatrix(DR, Dw)
end

"""
    get_w(R)

Compute the angular velocity `w` from the rotation matrix `R` and its derivative `DR = D.(R)`.
"""
function get_w(R::AbstractMatrix)
    R = collect(R)
    DR = collect(D.(R))
    [R[3, :]'DR[2, :], -R[3, :]'DR[1, :], R[2, :]'DR[1, :]] |> collect
end

function get_w(RM)
    R = RM.R
    DR = D.(RM.R)
    [R[3, :]'DR[2, :], -R[3, :]'DR[1, :], R[2, :]'DR[1, :]]
end

function get_w(Q::AbstractVector)
    Q = collect(Q)
    DQ = collect(D.(Q))
    angular_velocity2(Q, DQ)
end

"""
    h2 = resolve2(R21, h1)

`R21` is a 3x3 matrix that transforms a vector from frame 1 to frame 2. `h1` is a
vector resolved in frame 1. `h2` is the same vector in frame 2.

Typical usage:
```julia
resolve2(ori(frame_a), a_0 - g_0)
```
"""
resolve2(R21::RotationMatrix, v1) = R21 * collect(v1)

"""
    h1 = resolve1(R21, h2)

`R12` is a 3x3 matrix that transforms a vector from frame 1 to frame 2. `h2` is a
vector resolved in frame 2. `h1` is the same vector in frame 1.

Typical usage:
```julia
resolve1(ori(frame_a), r_ab)
```
"""
resolve1(R21::RotationMatrix, v2) = R21'collect(v2)

resolve1(sys::ODESystem, v) = resolve1(ori(sys), v)
resolve2(sys::ODESystem, v) = resolve2(ori(sys), v)

function resolve_relative(v1, R1, R2)
    R1 isa ODESystem && (R1 = ori(R1))
    R2 isa ODESystem && (R2 = ori(R2))
    v2 = resolve2(R2, resolve1(R1, v1))
end

skew(s) = [0 -s[3] s[2]; s[3] 0 -s[1]; -s[2] s[1] 0]
skewcoords(R::AbstractMatrix) = [R[3, 2]; R[1, 3]; R[2, 1]]

function planar_rotation(axis, phi, der_angle)
    length(axis) == 3 || error("axis must be a 3-vector")
    axis = collect(axis)
    ee = collect(axis * axis')
    R = ee + (I(3) - ee) * cos(phi) - skew(axis) * sin(phi)
    w = axis * phi
    RotationMatrix(R, w)
end

"""
    R2 = absolute_rotation(R1, R_rel)

- `R1`: `Orientation` object to rotate frame 0 into frame 1
- `R_rel`: `Orientation` object to rotate frame 1 into frame 2
- `R2`: `Orientation` object to rotate frame 0 into frame 2
"""
function absolute_rotation(R1, R_rel)
    # R2 = R_rel.R*R1.R
    # w = resolve2(R_rel, R1.w) + R_rel.w
    # RotationMatrix(R2, w)
    R1 isa ODESystem && (R1 = ori(R1))
    R_rel isa ODESystem && (R_rel = ori(R_rel))
    R_rel * R1
end

function relative_rotation(R1, R2)
    R1 isa ODESystem && (R1 = ori(R1))
    R2 isa ODESystem && (R2 = ori(R2))
    R = R2'R1
    w = R2.w - resolve2(R2, resolve1(R1, R1.w))
    RotationMatrix(R.R, w)
end

function inverse_rotation(R)
    R isa ODESystem && (R = ori(R))
    Ri = R.R'
    wi = -resolve1(R, R.w)
    RotationMatrix(Ri, wi)
end

function Base.:~(R1::RotationMatrix, R2::RotationMatrix)
    # [vec(R1.R.mat .~ R2.R.mat);
    #     R1.w .~ R2.w]
    vec(R1.R.mat .~ R2.R.mat)
end

"""
    connect_orientation(R1,R2; iscut=false)

Connect two rotation matrices together, optionally introducing a cut joint. A normal connection of two rotation matrices introduces 9 constraints, while a cut connection introduces 3 constraints only. This is useful to open kinematic loops, see [Using cut joints](@ref) (docs page) for an example where this is used.
"""
function connect_orientation(R1,R2; iscut=false)
    if iscut
        residue(R1, R2) .~ 0
    else
        R1 ~ R2
    end
end

function angular_velocity2(R::RotationMatrix, w=R.w)
    w
end

function angular_velocity1(R::RotationMatrix, w=R.w)
    resolve1(R, w)
end

function orientation_constraint(R::RotationMatrix)
    T = R.R
    [T[:, 1]'T[:, 1] - 1
     T[:, 2]'T[:, 2] - 1
     T[:, 3]'T[:, 3] - 1
     T[:, 1]'T[:, 2]
     T[:, 1]'T[:, 3]
     T[:, 2]'T[:, 3]]
end

orientation_constraint(R1, R2) = orientation_constraint(R1'R2)

function residue(R1, R2)
    # https://github.com/modelica/ModelicaStandardLibrary/blob/master/Modelica/Mechanics/MultiBody/Frames/Orientation.mo
    R1 isa ODESystem && (R1 = ori(R1))
    R2 isa ODESystem && (R2 = ori(R2))
    R1 = R1.R
    R2 = R2.R
    [atan(cross(R1[1, :], R1[2, :]) ⋅ R2[2, :], R1[1, :] ⋅ R2[1, :])
     atan(-cross(R1[1, :], R1[2, :]) ⋅ R2[1, :], R1[2, :] ⋅ R2[2, :])
     atan(R1[2, :] ⋅ R2[1, :], R1[3, :] ⋅ R2[3, :])]

    # [
    # cross(R1[1, :], R1[2, :])'R2[2, :]
    # -cross(R1[1, :], R1[2, :])'R2[1, :]
    # R1[2, :]'R2[1, :]
    # ]
end

function connect_loop(F1, F2)
    F1.metadata[:loop_opening] = true
    # connect(F1, F2)
    # orientation_constraint(ori(F1)'ori(F2)) .~ 0
    residue(F1, F2) .~ 0
end

## Quaternions
struct Quaternion <: Orientation
    Q
    w::Any
end

Base.getindex(Q::Quaternion, i) = Q.Q[i]

"""
Never call this function directly from a component constructor, instead call `f = Frame(); R = ori(f)` and add `f` to the subsystems.
"""
function NumQuaternion(; Q = [1.0, 0, 0, 0.0], w = zeros(3), name, varw = false)
    # The reason for not calling this directly is that all R vaiables have to have the same name since they are treated as connector variables (otherwise a connection error is thrown). A component with more than one rotation matrix will thus have two different R variables that overwrite each other
    # Q = at_variables_t(:Q, 1:4, default = Q) #[description="Orientation rotation matrix ∈ SO(3)"]
    @variables Q(t)[1:4] = [1.0,0,0,0]
    if varw
        @variables w(t)[1:3]=w [description="angular velocity"]
        # w = at_variables_t(:w, 1:3, default = w)
    else
        w = get_w(Q)
    end
    Q, w = collect.((Q, w))
    Quaternion(Q, w)
end


orientation_constraint(q::AbstractVector) = q'q - 1
orientation_constraint(q::Quaternion) = orientation_constraint(q.Q)

# function angular_velocity2(q::AbstractVector, q̇)
#     Q = [q[4] q[3] -q[2] -q[1]; -q[3] q[4] q[1] -q[2]; q[2] -q[1] q[4] -q[3]]
#     2 * Q * q̇
# end

Base.:/(q::Rotations.Quaternions.Quaternion, x::Num) = Rotations.Quaternions.Quaternion(q.s / x, q.v1 / x, q.v2 / x, q.v3 / x)
function from_Q(Q2, w)
    # Q2 = to_q(Q) # Due to different conventions
    q = Rotations.QuatRotation(Q2)
    R = RotMatrix(q)
    RotationMatrix(R, w)
end

to_q(Q::AbstractVector) = SA[Q[4], Q[1], Q[2], Q[3]]
to_q(Q::Rotations.QuatRotation) = to_q(vec(Q))
to_mb(Q::AbstractVector) = SA[Q[2], Q[3], Q[4], Q[1]]
to_mb(Q::Rotations.QuatRotation) = to_mb(vec(Q))
Base.vec(Q::Rotations.QuatRotation) = SA[Q.q.s, Q.q.v1, Q.q.v2, Q.q.v3]

# function angular_velocity1(Q, der_Q)
#     2*([Q[4] -Q[3] Q[2] -Q[1]; Q[3] Q[4] -Q[1] -Q[2]; -Q[2] Q[1] Q[4] -Q[3]]*der_Q)
# end

# function angular_velocity2(Q, der_Q)
#     2*([Q[4]  Q[3] -Q[2] -Q[1]; -Q[3] Q[4] Q[1] -Q[2]; Q[2] -Q[1] Q[4] -Q[3]]*der_Q)
# end


## Euler

function axes_rotations(sequence, angles, der_angles, name = :R_ar)
    R = axis_rotation(sequence[3], angles[3]) *
        axis_rotation(sequence[2], angles[2]) *
        axis_rotation(sequence[1], angles[1])

    w = axis(sequence[3]) * der_angles[3] +
        resolve2(axis_rotation(sequence[3], angles[3]), axis(sequence[2]) * der_angles[2]) +
        resolve2(axis_rotation(sequence[3], angles[3]) *
                 axis_rotation(sequence[2], angles[2]),
                 axis(sequence[1]) * der_angles[1])
    RotationMatrix(R.R, w)
end

axis(s) = float.(s .== (1:3))

"""
    axis_rotation(sequence, angle; name = :R)

Generate a rotation matrix for a rotation around the specified axis.

- `sequence`: The axis to rotate around (1: x-axis, 2: y-axis, 3: z-axis)
- `angle`: The angle of rotation (in radians)

Returns a `RotationMatrix` object.
"""
function axis_rotation(sequence, angle; name = :R)
    if sequence == 1
        return RotationMatrix(rotx(angle), zeros(3))
    elseif sequence == 2
        return RotationMatrix(roty(angle), zeros(3))
    elseif sequence == 3
        return RotationMatrix(rotz(angle), zeros(3))
    else
        error("Invalid sequence $sequence")
    end
end

"""
    rotx(t, deg = false)

Generate a rotation matrix for a rotation around the x-axis.

- `t`: The angle of rotation (in radians, unless `deg` is set to true)
- `deg`: (Optional) If true, the angle is in degrees

Returns a 3x3 rotation matrix.
"""
function rotx(t, deg = false)
    if deg
        t *= pi / 180
    end
    ct = cos(t)
    st = sin(t)
    R = [1 0 0
         0 ct -st
         0 st ct]
end

"""
    roty(t, deg = false)

Generate a rotation matrix for a rotation around the y-axis.

- `t`: The angle of rotation (in radians, unless `deg` is set to true)
- `deg`: (Optional) If true, the angle is in degrees

Returns a 3x3 rotation matrix.
"""
function roty(t, deg = false)
    if deg
        t *= pi / 180
    end
    ct = cos(t)
    st = sin(t)
    R = [ct 0 st
         0 1 0
         -st 0 ct]
end

"""
    rotz(t, deg = false)

Generate a rotation matrix for a rotation around the z-axis.

- `t`: The angle of rotation (in radians, unless `deg` is set to true)
- `deg`: (Optional) If true, the angle is in degrees

Returns a 3x3 rotation matrix.
"""
function rotz(t, deg = false)
    if deg
        t *= pi / 180
    end
    ct = cos(t)
    st = sin(t)
    R = [ct -st 0
         st ct 0
         0 0 1]
end

function from_nxy(n_x, n_y)
    e_x = norm(n_x) < 1e-10 ? [1.0, 0, 0] : _normalize(n_x)
    e_y = norm(n_y) < 1e-10 ? [0, 1.0, 0] : _normalize(n_y)
    n_z_aux = cross(e_x, e_y)
    n_y_aux = #if n_z_aux' * n_z_aux > 1.0e-6 # TODO: MTK too buggy with ifelse to handle this logic
        e_y
    # elseif abs(e_x[1]) > 1.0e-6
    #     [0, 1.0, 0]
    # else
    #     [1.0, 0, 0]
    # end
    e_z_aux = cross(e_x, n_y_aux)
    e_z = _normalize(e_z_aux)
    RotationMatrix([e_x cross(e_z, e_x) e_z]', zeros(3))
end

# function from_nxy(n_x, n_y)
#     e_x = norm(n_x) < 1e-10 ? [1.0, 0, 0] : _normalize(n_x)
#     e_y = norm(n_y) < 1e-10 ? [0, 1.0, 0] : _normalize(n_y)
#     n_z_aux = cross(e_x, e_y)
#     n_y_aux = ifelse(
#         n_z_aux' * n_z_aux > 1.0e-6,
#         e_y,
#         ifelse(
#             abs(e_x[1]) > 1.0e-6,
#             [0, 1.0, 0],
#             [1.0, 0, 0])
#     )
#     e_z_aux = cross(e_x, n_y_aux)
#     e_z = _normalize(e_z_aux)
#     RotationMatrix([e_x e_y e_z], zeros(3))
# end

function resolve_dyade1(R, D2)
    R'D2*R
end

function resolve_dyade2(R, D1)
    R*D1*R'
end

"""
    R_W_F = get_rot(sol, frame, t)

Extract a 3×3 rotation matrix ∈ SO(3) from a solution at time `t`.

The rotation matrix returned, ``R_W^F``, is such that when a vector ``p_F`` expressed in the local `frame` of reference ``F`` is multiplied by ``R_W^F`` as ``Rp``, the resulting vector is ``p_W`` expressed in the world frame:
```math
p_W = R_W^F  p_F
```

The columns of ``R_W_F`` indicate are the basis vectors of the frame ``F`` expressed in the world coordinate frame.

See also [`get_trans`](@ref), [`get_frame`](@ref), [Orientations and directions](@ref) (docs section).
"""
function get_rot(sol, frame, t)
    Rotations.RotMatrix3(reshape(sol(t, idxs = vec(ori(frame).R.mat')), 3, 3))
end

"""
    get_trans(sol, frame, t)

Extract the translational part of a frame from a solution at time `t`.
See also [`get_rot`](@ref), [`get_frame`](@ref), [Orientations and directions](@ref) (docs section).
"""
function get_trans(sol, frame, t)
    SVector{3}(sol(t, idxs = collect(frame.r_0)))
end

"""
    T_W_F = get_frame(sol, frame, t)

Extract a 4×4 transformation matrix ∈ SE(3) from a solution at time `t`.

The transformation matrix returned, ``T_W^F``, is such that when a homogenous-coordinate vector ``p_F``, expressed in the local `frame` of reference ``F`` is multiplied by ``T_W^F`` as ``Tp``, the resulting vector is ``p_W`` expressed in the world frame:
```math
p_W = T_W^F  p_F
```

See also [`get_trans`](@ref) and [`get_rot`](@ref), [Orientations and directions](@ref) (docs section).
"""
function get_frame(sol, frame, t)
    R = get_rot(sol, frame, t)
    tr = get_trans(sol, frame, t)
    [R tr; 0 0 0 1]
end

function nonunit_quaternion_equations(R, w)
    @variables Q(t)[1:4]=[1,0,0,0], [description="Unit quaternion with [w,i,j,k]"] # normalized
    @variables Q̂(t)[1:4]=[1,0,0,0], [description="Non-unit quaternion with [w,i,j,k]"] # Non-normalized
    @variables n(t)=1 c(t)=0
    @parameters k = 0.1
    Q̂ = collect(Q̂)
    Q = collect(Q)
    # w is used in Ω, and Ω determines D(Q̂)
    # This corresponds to modelica's 
    # frame_a.R = from_Q(Q, angularVelocity2(Q, der(Q)));
    # where angularVelocity2(Q, der(Q)) = 2*([Q[4]  Q[3] -Q[2] -Q[1]; -Q[3] Q[4] Q[1] -Q[2]; Q[2] -Q[1] Q[4] -Q[3]]*der_Q)
    # They also have w_a = angularVelocity2(frame_a.R) even for quaternions, so w_a = angularVelocity2(Q, der(Q)), this is their link between w_a and D(Q), while ours is D(Q̂) .~ (Ω * Q̂)
    Ω = [0 -w[1] -w[2] -w[3]; w[1] 0 w[3] -w[2]; w[2] -w[3] 0 w[1]; w[3] w[2] -w[1] 0]
    # QR = from_Q(Q, angular_velocity2(Q, D.(Q)))
    QR = from_Q(Q, w)
    [
        n ~ Q̂'Q̂
        c ~ k * (1 - n)
        D.(Q̂) .~ (Ω' * Q̂) ./ 2 + c * Q̂ # We use Ω' which is the same as using -w to handle the fact that w coming in here is typically described frame_a rather than in frame_b, the paper is formulated with w being expressed in the rotating body frame (frame_b)
        Q .~ Q̂ ./ sqrt(n)
        R ~ QR
    ]
end