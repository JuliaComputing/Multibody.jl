using Test

Ts = 0.001
t = -1:Ts:3
q, qd, qdd = point_to_point(t)

q1 = [1, 1.2]
qd_max = [0.7, 1.2]
qdd_max = [0.9, 1.1]
q, qd, qdd = point_to_point(t; q1, qd_max)
@test all(q .<= q1')
@test all(qd .<= qd_max')
@test all(qdd.<= qdd_max')

_, _, _, t1 = point_to_point(0; q1, qd_max)
ci = findlast(qd .!= 0)
t02 = t[ci[1]]
@test t02 ≈ t1 rtol=1e-3

_qd = Multibody.centraldiff(q) ./ Ts
_qdd = Multibody.centraldiff(qd)./ Ts

# Plots.plot(t, [q qd qdd], layout=(2,1))
# Plots.plot!(t, [_qd _qdd], sp=[1 2])


N = 200
t0 = randn()
q0 = randn(N)
q1 = randn(N)
qd_max = rand(N)
qdd_max = rand(N)
O = zeros(N)
q, qd, qdd, t1 = point_to_point(0; q0, q1, qd_max, qdd_max, t0)
Ts = 0.01
Tf = t1 + 0.1
t = (t0 - 1):Ts:Tf
q, qd, qdd = point_to_point(t; q0, q1, qd_max, qdd_max, t0)
@test q[1,:] ≈ q0 norm=maximum atol=1e-3
@test q[end,:] ≈ q1 norm=maximum atol=1e-3
@test qd[1,:] == O
@test qd[end,:] == O
@test qdd[1,:] == O
@test qdd[end,:] == O
@test all(min.(q0, q1)' .- 1e-4 .<= q .<= max.(q0, q1)' .+ 1e-4)
@test all(abs.(qd) .<= 1.001.*qd_max')
@test all(abs.(qdd) .<= 1.001.*qdd_max')


