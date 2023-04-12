using Rotations

const R3{T} = RotMatrix{3, T}

abstract type Orientation end

struct RotationMatrix <: Orientation
    R::R3
    w
end

RotationMatrix(R::AbstractMatrix, w) = RotationMatrix(R3(R), w)

RotationMatrix() = RotationMatrix(R3(1.0I(3)), zeros(3))

function NumRotationMatrix(; R = collect(1.0I(3)), w = zeros(3), name, derived_w = true)
    R = at_variables_t(:R, 1:3, 1:3) #[description="Orientation rotation matrix ∈ SO(3)"]
    # @variables w(t)[1:3]=w [description="angular velocity"]
    # R = collect(R)
    # R = ModelingToolkit.renamespace.(name, R) .|> Num
    if derived_w
        w = get_w(R)
    else
        w = at_variables_t(:w, 1:3)
    end
    R,w = collect.((R,w))
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

Base.:*(R1::RotationMatrix, x::AbstractVector) = R1.R*x
Base.:*(R1::RotationMatrix, R2::RotationMatrix) = RotationMatrix(R1.R.mat*R2.R.mat, R1*R2.w + collect(R1.w))
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
"""
resolve2(R21::RotationMatrix, v1) = R21 * collect(v1)

"""
    h1 = resolve1(R21, h2)

`R12` is a 3x3 matrix that transforms a vector from frame 1 to frame 2. `h2` is a
vector resolved in frame 2. `h1` is the same vector in frame 1.
"""
resolve1(R21::RotationMatrix, v2) = R21'collect(v2)



skew(s) = [0 -s[3] s[2];s[3] 0 -s[1]; -s[2] s[1] 0]
skewcoords(R::AbstractMatrix) = [R[3,2];R[1,3];R[2,1]]

function planar_rotation(axis, ϕ, ϕ̇)
    length(axis) == 3 || error("axis must be a 3-vector")
    axis = collect(axis)
    ee = collect(axis*axis')
    R = ee + (I(3) - ee)*cos(ϕ) - skew(axis)*sin(ϕ)
    w = axis*ϕ̇
    RotationMatrix(R, w)
end

"""
    R2 = abs_rotation(R1, R_rel)

- `R1`: `Orientation` object to rotate frame 0 into frame 1
- `R_rel`: `Orientation` object to rotate frame 1 into frame 2
- `R2`: `Orientation` object to rotate frame 0 into frame 2
"""
function abs_rotation(R1, R_rel)
    # R2 = R_rel.R*R1.R
    # w = resolve2(R_rel, R1.w) + R_rel.w
    # RotationMatrix(R2, w)
    R_rel*R1
end

function Base.:~(R1::RotationMatrix, R2::RotationMatrix)
    # [vec(R1.R.mat .~ R2.R.mat);
    #     R1.w .~ R2.w]
    vec(R1.R.mat .~ R2.R.mat)
end

function angular_velocity2(R::RotationMatrix)
    R.w
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
    [
        atan(cross(R1[1, :], R1[2, :])⋅R2[2, :],R1[1,:]⋅R2[1,:])
        atan(-cross(R1[1, :],R1[2, :])⋅R2[1, :],R1[2,:]⋅R2[2,:])
        atan(R1[2, :]⋅R2[1, :],R1[3,:]⋅R2[3,:])
    ]
end

## Quaternions
orientation_constraint(q::AbstractVector) = q'q - 1

function angular_velocity2(q::AbstractVector, q̇) 
    Q = [q[4] q[3] -q[2] -q[1]; -q[3] q[4] q[1] -q[2]; q[2] -q[1] q[4] -q[3]]
    2*Q*q̇
end