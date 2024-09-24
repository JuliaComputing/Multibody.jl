# Trajectory_planning

Two methods of planning trajectories are available
- [`point_to_point`](@ref): Generate a minimum-time point-to-point trajectory with specified start and endpoints, not exceeding specified speed and acceleration limits.
- [`traj5`](@ref): Generate a 5:th order polynomial trajectory with specified start and end points. Additionally allows specification of start and end values for velocity and acceleration.

Components that make use of these trajectory generators is provided:
- [`KinematicPTP`](@ref)
- [`Kinematic5`](@ref)

These both have output connectors of type `RealOutput` called `q, qd, qdd` for positions, velocities and accelerations.

See [Industrial robot](@ref) for an example making use of the [`point_to_point`](@ref) planner.

## Docstrings

```@index
```


```@autodocs
Modules = [Multibody]
Pages   = ["path_planning.jl", "ptp.jl"]
```