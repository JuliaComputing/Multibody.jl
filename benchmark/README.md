# JSC-era benchmark of Multibody.jl

This branch (`jsc-benchmark`) reproduces timings of Multibody.jl when it still
used JuliaSimCompiler (JSC) as the simplification + codegen backend, after JSC
was deprecated. The branch is based on commit `98c867f` (2024-07-17,
Multibody.jl `v0.1.1`), which is squarely in the JSC 0.1.15 era.

## What you need

- `git` with access to the private repo
  `JuliaComputing/Multibody-JSC-Benchmark-Vendor` (vendored as the `vendor/`
  submodule).
- `juliaup` with Julia `1.10` installed (`juliaup add 1.10`).
- Network access to the General registry — every non-JSC dep is pinned to a
  version that is still in General.

You do **not** need JuliaHubRegistry access. The vendored submodule supplies
the patched JuliaSimCompiler / JuliaSimBase / JuliaSimCompilerRuntime trio
directly.

## Setup (one-time)

```bash
git clone --branch jsc-benchmark --recurse-submodules \
    git@github.com:JuliaComputing/Multibody.jl.git Multibody.jl-jsc
cd Multibody.jl-jsc
julia +1.10 --project=benchmark -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

(If you already cloned without `--recurse-submodules`, run
`git submodule update --init --recursive` from inside the clone.)

The `benchmark/Manifest.toml` pins:

| Package | Version |
|---|---|
| Julia | 1.10 |
| JuliaSimCompiler | 0.1.15 (vendored, patched) |
| JuliaSimBase | 1.1.3 (vendored) |
| JuliaSimCompilerRuntime | 1.0.2 (vendored) |
| ModelingToolkit | 9.26.0 |
| ModelingToolkitStandardLibrary | 2.10.0 |
| OrdinaryDiffEq | 6.87.0 |
| SciMLBase | 2.46.0 |
| Symbolics / SymbolicUtils | 5.34.0 / 2.1.2 |
| NonlinearSolve / DiffEqBase | 3.13.1 / 6.151.5 |
| Multibody | dev-tracked at this branch's root (`v0.1.1`) |

## Run the robot benchmark

```bash
julia +1.10 --project=benchmark benchmark/robot_bench.jl
```

The script runs the `Robot6DOF` example twice and prints `@time` for
`structural_simplify(IRSystem(robot))`, `ODEProblem` construction, and `solve`,
then a `@btime` of `solve` on the already-compiled problem.

Reference numbers from the original benchmarking machine (Julia 1.10.11,
x86_64 Linux):

```
========== WARM-UP RUN (includes compilation) ==========
structural_simplify (warm-up):  27.5 s  (45.4 M alloc, 2.98 GiB)
ODEProblem creation (warm-up):   1.7 s  ( 4.9 M alloc,  333 MiB)
solve (warm-up):                15.1 s  (18.7 M alloc, 1.06 GiB)

========== TIMED RUN (cold structural_simplify/ODEProblem) ==========
structural_simplify:             3.7 s  (28.7 M alloc, 1.88 GiB)
ODEProblem creation:             0.32 s ( 1.7 M alloc,  103 MiB)
first solve:                     8.9 s  (10.3 M alloc,  521 MiB) [98% compile]

========== BenchmarkTools @btime (solve only, already-compiled) ==========
                                 144 ms (2,274,911 alloc, 39.68 MiB)
```

The 144 ms / 39.7 MiB on the warm `solve` matches the historical in-source
comment in `test/test_robot.jl` (`152.225 ms (2272926 allocations: 40.08 MiB)`)
on this hardware, confirming the recovered environment reproduces the original
performance.

## Patches applied to the vendored JSC

See `vendor/README.md`. In short:

1. `scheduled_system.jl` `id_to_clock::Union{Nothing, Vector{SciMLBase.TimeDomain}}`
   → `Union{Nothing, Vector}` (MTK 9.26 has its own `TimeDomain` abstract type
   that does not unify with `SciMLBase.TimeDomain`).
2. `JuliaSimCompiler.jl` SnoopPrecompile workload disabled (hits the same
   constructor mismatch during package precompile and is not needed for the
   benchmark).

Neither patch touches the structural-simplify or solve hot paths.

## Notes / pitfalls

- The branch is based on `98c867f`, *not* the last JSC-era commit `8e3f8c1`,
  because Multibody at `8e3f8c1` uses MTK 9.42-only syntax in PlanarMechanics
  (`@mtkmodel` defaults like `[0, 0]`) that MTK 9.26 cannot parse. JSC 0.1.15
  itself is internally incompatible with MTK 9.32+ (the `ScheduledSystem`
  constructor changed), so we sit at MTK 9.26 and Multibody `v0.1.1`.
- The vendored JSC trio is proprietary. The vendor repo
  `JuliaComputing/Multibody-JSC-Benchmark-Vendor` is **private**.
