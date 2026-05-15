using Multibody
using ModelingToolkit
using JuliaSimCompiler
t = Multibody.t

# WheelInWorld adapted from test_wheels.jl @ 38abc92 to the 98c867f API
# (RollingWheel here takes x0/y0, not x0/z0 — the y/z convention changed later).
@mtkmodel WheelInWorld begin
    @components begin
        world = World()
        wheel = RollingWheel(
            radius     = 0.3,
            m          = 2,
            I_axis     = 0.06,
            I_long     = 0.12,
            x0         = 0.2,
            y0         = 0.2,
            der_angles = [0, 5, 1],
        )
    end
end

@named worldwheel = WheelInWorld()
worldwheel = complete(worldwheel)

@info "structural_simplify(IRSystem(worldwheel))"
@time "structural_simplify" ssys = structural_simplify(IRSystem(worldwheel))

println("\n=== unknowns(ssys) — JSC's state realization ===")
us = unknowns(ssys)
for (i, u) in enumerate(us)
    println(lpad(i, 3), ": ", u)
end
println("\nTotal: ", length(us), " state variables")
