#=
# Sloshing Water Tank on a Cart

This example demonstrates the CraigBampton component by modeling a cart on a
horizontal prismatic joint with a water tank on top. The water sloshing dynamics
are modeled using Craig-Bampton modal reduction.

## Physical Setup

```
    ┌─────────────────┐
    │   Water Tank    │  ← Sloshing water (CB model)
    │   ~~~~~~~~~~~~  │
    │                 │
    └────────┬────────┘
             │
    ┌────────┴────────┐
    │      Cart       │  ← Rigid body
    └────────┬────────┘
             │
    ═══════════════════  ← Prismatic joint (horizontal)
         World
```

## Sloshing Physics

For a rectangular tank with length L and water depth h, the sloshing frequency is:
    ω_n = sqrt(g * π * n / L * tanh(π * n * h / L))

The participating (effective sloshing) mass for mode n is:
    m_n = 8 * ρ * L * W * h / (π³ * n³) * tanh(π * n * h / L)
=#

using Multibody
using ModelingToolkit
using OrdinaryDiffEq
using LinearAlgebra
using Plots

t = Multibody.t
D = Differential(t)
world = Multibody.world

# =============================================================================
# Tank and fluid parameters
# =============================================================================

# Tank geometry
L_tank = 1.0        # Tank length in x-direction (m)
W_tank = 0.5        # Tank width in y-direction (m)
h_water = 0.4       # Water depth (m)
h_tank = 0.5        # Tank height (m)
t_wall = 0.01       # Wall thickness (m)

# Material properties
ρ_water = 1000.0    # Water density (kg/m³)
ρ_tank = 2700.0     # Tank material density (aluminum, kg/m³)
g = 9.81            # Gravitational acceleration (m/s²)

# Cart parameters
m_cart = 50.0       # Cart mass (kg)
cart_length = 0.4   # Cart length (m)

# Damping
ζ_slosh = 0.02      # Sloshing damping ratio (2%)

# =============================================================================
# Compute sloshing model parameters
# =============================================================================

# Total water mass
m_water = ρ_water * L_tank * W_tank * h_water
println("Water mass: $(m_water) kg")

# Approximate tank structure mass (simplified: just the walls)
m_tank_structure = ρ_tank * t_wall * (
    2 * L_tank * h_tank +           # Front and back walls
    2 * W_tank * h_tank +           # Side walls
    L_tank * W_tank                  # Bottom
) * W_tank  # This is approximate
m_tank_structure = 20.0  # Simplified: assume 20 kg tank structure

# Sloshing frequencies for first two modes
function sloshing_frequency(n, L, h, g)
    sqrt(g * π * n / L * tanh(π * n * h / L))
end

ω_1 = sloshing_frequency(1, L_tank, h_water, g)
ω_2 = sloshing_frequency(2, L_tank, h_water, g)
println("Sloshing frequencies: ω₁ = $(round(ω_1, digits=2)) rad/s, ω₂ = $(round(ω_2, digits=2)) rad/s")
println("Sloshing periods: T₁ = $(round(2π/ω_1, digits=2)) s, T₂ = $(round(2π/ω_2, digits=2)) s")

# Participating mass for each mode (mass that effectively sloshes)
function participating_mass(n, ρ, L, W, h)
    8 * ρ * L * W * h / (π^3 * n^3) * tanh(π * n * h / L)
end

m_1 = participating_mass(1, ρ_water, L_tank, W_tank, h_water)
m_2 = participating_mass(2, ρ_water, L_tank, W_tank, h_water)
println("Participating masses: m₁ = $(round(m_1, digits=2)) kg, m₂ = $(round(m_2, digits=2)) kg")
println("Total participating: $(round(m_1 + m_2, digits=2)) kg ($(round(100*(m_1+m_2)/m_water, digits=1))% of water)")

# =============================================================================
# Build Craig-Bampton matrices for sloshing model
# =============================================================================

n_modes = 2
n_boundaries = 2  # Tank bottom (attached) and top (free)
n_bdof = 6 * n_boundaries

# Total mass at boundary (tank + water)
m_total = m_tank_structure + m_water

# Approximate moments of inertia for the tank+water assembly
I_xx = m_total * (W_tank^2 + h_tank^2) / 12
I_yy = m_total * (L_tank^2 + h_tank^2) / 12
I_zz = m_total * (L_tank^2 + W_tank^2) / 12

# M_BB: Boundary mass matrix
# For boundary 1 (bottom, attached to cart): full mass
# For boundary 2 (top, free): small mass (just for numerical stability)
M_BB = zeros(n_bdof, n_bdof)
# Boundary 1 (indices 1-6): main mass
M_BB[1,1] = m_total      # x translation
M_BB[2,2] = m_total      # y translation
M_BB[3,3] = m_total      # z translation
M_BB[4,4] = I_xx         # rotation about x
M_BB[5,5] = I_yy         # rotation about y
M_BB[6,6] = I_zz         # rotation about z
# Boundary 2 (indices 7-12): small mass for free end
M_BB[7,7] = 0.1
M_BB[8,8] = 0.1
M_BB[9,9] = 0.1
M_BB[10,10] = 0.01
M_BB[11,11] = 0.01
M_BB[12,12] = 0.01

# M_BI: Boundary-modal coupling
# Horizontal (x) acceleration at boundary 1 excites sloshing modes
M_BI = zeros(n_bdof, n_modes)
M_BI[1, 1] = m_1    # x-acceleration couples to mode 1
M_BI[1, 2] = m_2    # x-acceleration couples to mode 2

# M_II: Modal mass matrix (identity for mass-normalized modes)
M_II = Matrix(1.0I, n_modes, n_modes)

# K_BB: Boundary stiffness (zero for rigid connection)
K_BB = zeros(n_bdof, n_bdof)
# Add small stiffness for the free boundary to prevent drift
K_BB[7,7] = 100.0
K_BB[8,8] = 100.0
K_BB[9,9] = 100.0
K_BB[10,10] = 10.0
K_BB[11,11] = 10.0
K_BB[12,12] = 10.0

# K_II: Modal stiffness matrix (diagonal with ω²)
ω_modes = [ω_1, ω_2]
K_II = Diagonal(ω_modes.^2)

# Damping ratios
ζ_modes = [ζ_slosh, ζ_slosh]

# =============================================================================
# Create the multibody model
# =============================================================================

systems = @named begin
    # Prismatic joint for cart horizontal motion
    prismatic = Prismatic(n=[1, 0, 0], axisflange=true)

    # Cart body
    cart = BodyShape(m=m_cart, r=[cart_length, 0, 0], radius=0.05, color=[0.3, 0.3, 0.3, 1])

    # Mounting point on top of cart for tank
    mount = FixedTranslation(r=[cart_length/2, h_tank/2, 0])

    # Sloshing water tank (Craig-Bampton model)
    tank = CraigBampton(
        boundary_positions = [[0.0, 0.0, 0.0], [0.0, h_tank, 0.0]],  # Bottom and top
        M_BB = M_BB,
        M_BI = M_BI,
        M_II = M_II,
        K_BB = K_BB,
        K_II = K_II,
        ζ = ζ_modes,
        ω = ω_modes,
    )

    # Optional: damper on prismatic joint for cart friction
    cart_damper = Multibody.Translational.Damper(d=5.0)
end

# Connect the system
eqs = [
    connect(world.frame_b, prismatic.frame_a)
    connect(prismatic.frame_b, cart.frame_a)
    connect(cart.frame_b, mount.frame_a)
    connect(mount.frame_b, tank.frame_a)
    connect(prismatic.axis, cart_damper.flange_a)
    connect(prismatic.support, cart_damper.flange_b)
]

@named model = System(eqs, t; systems=[world; systems])

# =============================================================================
# Simulate: Free response - initial cart velocity
# =============================================================================

println("\n=== Building and simulating model ===")
ssys = multibody(model)

# Initial conditions: cart has initial velocity
v0_cart = 1.0  # m/s initial velocity
prob = ODEProblem(ssys, [
    prismatic.v => v0_cart,
    prismatic.s => 0.0,
], (0.0, 20.0))

sol = solve(prob, Tsit5())
println("Simulation complete: $(sol.retcode)")

# =============================================================================
# Plot results
# =============================================================================

# Cart position and velocity
p1 = plot(sol, idxs=[prismatic.s], label="Cart position (m)", ylabel="Position")
p2 = plot(sol, idxs=[prismatic.v], label="Cart velocity (m/s)", ylabel="Velocity")

# Sloshing mode amplitudes
p3 = plot(sol, idxs=[tank.η[1]], label="Mode 1 (ω=$(round(ω_1,digits=1)) rad/s)", ylabel="η")
plot!(p3, sol, idxs=[tank.η[2]], label="Mode 2 (ω=$(round(ω_2,digits=1)) rad/s)")

# Combined plot
p = plot(p1, p2, p3, layout=(3,1), size=(800, 600),
         title=["Cart with Sloshing Tank" "" ""])
display(p)

# =============================================================================
# Analysis: Energy transfer
# =============================================================================

println("\n=== Energy Analysis ===")
# Initial kinetic energy of cart
KE_initial = 0.5 * m_cart * v0_cart^2
println("Initial cart KE: $(round(KE_initial, digits=2)) J")

# Final cart velocity (should be reduced due to sloshing)
v_final = sol[prismatic.v][end]
KE_final = 0.5 * m_cart * v_final^2
println("Final cart KE: $(round(KE_final, digits=2)) J")
println("Energy dissipated by sloshing damping: $(round(KE_initial - KE_final, digits=2)) J")
