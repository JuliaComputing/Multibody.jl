using Rotations

const R3{T} = RotMatrix{3, T}

abstract type Orientation end

struct RotationMatrix <: Orientation
    R::R3
    w::Any
end

RotationMatrix(R::AbstractMatrix, w) = RotationMatrix(R3(R), w)

RotationMatrix() = RotationMatrix(R3(1.0I(3)), zeros(3))

function NumRotationMatrix(; R = collect(1.0I(3)), w = zeros(3), name, varw = false)
    R = at_variables_t(:R, 1:3, 1:3, default = R) #[description="Orientation rotation matrix ∈ SO(3)"]
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

nullRotation() = RotationMatrix()

function ModelingToolkit.ODESystem(RM::RotationMatrix; name)
    # @variables R(t)[1:3, 1:3]=Matrix(RM) [description="Orientation rotation matrix ∈ SO(3)"]
    # @variables w(t)[1:3]=w [description="angular velocity"]
    # R,w = collect.((R,w))
    R = at_variables_t(:R, 1:3, 1:3)
    w = at_variables_t(:w, 1:3)

    defaults = Dict(R .=> RM)
    ODESystem(Equation[], t, [vec(R); w], []; name, defaults)
end

Base.:*(R1::RotationMatrix, x::AbstractVector) = R1.R * x
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

function resolveRelative(v1, R1, R2)
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
    R2 = absoluteRotation(R1, R_rel)

- `R1`: `Orientation` object to rotate frame 0 into frame 1
- `R_rel`: `Orientation` object to rotate frame 1 into frame 2
- `R2`: `Orientation` object to rotate frame 0 into frame 2
"""
function absoluteRotation(R1, R_rel)
    # R2 = R_rel.R*R1.R
    # w = resolve2(R_rel, R1.w) + R_rel.w
    # RotationMatrix(R2, w)
    R1 isa ODESystem && (R1 = ori(R1))
    R_rel isa ODESystem && (R_rel = ori(R_rel))
    R_rel * R1
end

function relativeRotation(R1, R2)
    R1 isa ODESystem && (R1 = ori(R1))
    R2 isa ODESystem && (R2 = ori(R2))
    R = R2'R1
    w = R2.w - resolve2(R2, resolve1(R1, R1.w))
    RotationMatrix(R.R, w)
end

function inverseRotation(R)
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

function angularVelocity2(R::RotationMatrix)
    R.w
end

function angular_velocity1(R::RotationMatrix)
    resolve1(R, R.w)
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

function residue(O1, O2)
    # https://github.com/modelica/ModelicaStandardLibrary/blob/master/Modelica/Mechanics/MultiBody/Frames/Orientation.mo
    R1 = O1.R
    R2 = O2.R
    [atan(cross(R1[1, :], R1[2, :]) ⋅ R2[2, :], R1[1, :] ⋅ R2[1, :])
     atan(-cross(R1[1, :], R1[2, :]) ⋅ R2[1, :], R1[2, :] ⋅ R2[2, :])
     atan(R1[2, :] ⋅ R2[1, :], R1[3, :] ⋅ R2[3, :])]
end

## Quaternions
orientation_constraint(q::AbstractVector) = q'q - 1

function angularVelocity2(q::AbstractVector, q̇)
    Q = [q[4] q[3] -q[2] -q[1]; -q[3] q[4] q[1] -q[2]; q[2] -q[1] q[4] -q[3]]
    2 * Q * q̇
end

function axesRotations(sequence, angles, der_angles, name = :R_ar)
    R = axisRotation(sequence[3], angles[3]) *
        axisRotation(sequence[2], angles[2]) *
        axisRotation(sequence[1], angles[1])

    w = axis(sequence[3]) * der_angles[3] +
        resolve2(axisRotation(sequence[3], angles[3]), axis(sequence[2]) * der_angles[2]) +
        resolve2(axisRotation(sequence[3], angles[3]) *
                 axisRotation(sequence[2], angles[2]),
                 axis(sequence[1]) * der_angles[1])
    RotationMatrix(R.R, w)
end

axis(s) = float.(s .== (1:3))

"""
    axisRotation(sequence, angle; name = :R)

Generate a rotation matrix for a rotation around the specified axis.

- `sequence`: The axis to rotate around (1: x-axis, 2: y-axis, 3: z-axis)
- `angle`: The angle of rotation (in radians)

Returns a `RotationMatrix` object.
"""
function axisRotation(sequence, angle; name = :R)
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
