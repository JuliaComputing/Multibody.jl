# Trajectory planning

Two methods of planning trajectories are available
- [`point_to_point`](@ref): Generate a minimum-time point-to-point trajectory with specified start and endpoints, not exceeding specified speed and acceleration limits.
- [`traj5`](@ref): Generate a 5:th order polynomial trajectory with specified start and end points. Additionally allows specification of start and end values for velocity and acceleration.

Components that make use of these trajectory generators is provided:
- [`KinematicPTP`](@ref)
- [`Kinematic5`](@ref)

These both have output connectors of type `RealOutput` called `q, qd, qdd` for positions, velocities and accelerations.

See [Industrial robot](@ref) for an example making use of the [`point_to_point`](@ref) planner.


## Example

### Point-to-point trajectory
```@example TRAJ
using Multibody, Plots
Ts = 0.001
t = -1:Ts:3

q1 = [1, 1.2]           # Final point (2 DOF)
qd_max = [0.7, 1.2]     # Max velocity (2 DOF)
qdd_max = [0.9, 1.1]    # Max acceleration (2 DOF)
q, qd, qdd = point_to_point(t; q1, qd_max, qdd_max)

plot(t, [q qd qdd], ylabel=["\$q\$" "\$\\dot{q}\$" "\$\\ddot{q}\$"], layout=(3,1), l=2, sp=[1 1 2 2 3 3], legend=false)
hline!([qd_max' qdd_max'], l=(2, :dash), sp=[2 2 3 3], c=[1 2 1 2], legend=false)
```

### 5:th order polynomial trajectory
```@example TRAJ
t = 0:Ts:3
q1 = 1
q, qd, qdd = traj5(t; q1)

plot(t, [q qd qdd], ylabel=["\$q\$" "\$\\dot{q}\$" "\$\\ddot{q}\$"], layout=(3,1), l=2, legend=false)
```



## Docstrings

```@index
Pages   = ["trajectory_planning.md"]
```


```@autodocs
Modules = [Multibody]
Pages   = ["path_planning.jl", "ptp.jl"]
```
