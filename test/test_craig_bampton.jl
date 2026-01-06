using Multibody
using ModelingToolkit
using OrdinaryDiffEq
using LinearAlgebra
using Test

t = Multibody.t
D = Differential(t)
world = Multibody.world


@testset "Component construction" begin
    # Test basic construction with 2 boundaries and 2 modes
    n_boundaries = 2
    n_modes = 2
    n_bdof = 6 * n_boundaries

    # Simple test matrices (identity-ish for testing)
    M_BB = 0.1 * I(n_bdof)
    M_BI = 0.01 * ones(n_bdof, n_modes)
    M_II = I(n_modes)  # Mass-normalized modes
    K_BB = 100.0 * I(n_bdof)
    ω_modes = [10.0, 25.0]  # rad/s
    K_II = Diagonal(ω_modes.^2)
    ζ = [0.02, 0.02]  # 2% damping

    @named cb = CraigBampton(
        boundary_positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        M_BB = M_BB,
        M_BI = M_BI,
        M_II = M_II,
        K_BB = K_BB,
        K_II = K_II,
        ζ = ζ,
        ω = ω_modes,
    )

    @test cb isa System
    @test hasproperty(cb, :frame_a)
    @test hasproperty(cb, :frame_b)
end

@testset "Simulation with fixed boundary" begin
    # Test a simple case: CB component fixed at frame_a, free at frame_b
    n_boundaries = 2
    n_modes = 2
    n_bdof = 6 * n_boundaries

    # Create simple CB matrices for a cantilever-like beam
    # Boundary 1 (frame_a) at origin, Boundary 2 (frame_b) at x=1
    M_BB = Diagonal(vcat(
        [1.0, 1.0, 1.0, 0.1, 0.1, 0.1],  # Boundary 1 (fixed, so won't matter much)
        [0.5, 0.5, 0.5, 0.05, 0.05, 0.05]  # Boundary 2 mass
    ))
    M_BI = zeros(n_bdof, n_modes)
    # Coupling: boundary 2 y-displacement couples to mode 1 (bending mode)
    M_BI[8, 1] = 0.1  # y-displacement at boundary 2, mode 1
    M_BI[9, 2] = 0.05  # z-displacement at boundary 2, mode 2

    M_II = I(n_modes)

    # Stiffness
    K_BB = Diagonal(vcat(
        [1e6, 1e6, 1e6, 1e4, 1e4, 1e4],  # Boundary 1 (high stiffness for fixed)
        [100.0, 100.0, 100.0, 10.0, 10.0, 10.0]  # Boundary 2 flexible
    ))

    ω_modes = [10.0, 25.0]  # Fundamental frequencies
    K_II = Diagonal(ω_modes.^2)
    ζ = [0.05, 0.05]  # 5% damping

    @named cb = CraigBampton(
        boundary_positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        M_BB = M_BB,
        M_BI = M_BI,
        M_II = M_II,
        K_BB = K_BB,
        K_II = K_II,
        ζ = ζ,
        ω = ω_modes,
    )

    # Fix frame_a to world
    eqs = [
        connect(world.frame_b, cb.frame_a)
    ]

    @named model = System(eqs, t; systems=[world, cb])

    # Try to build the system
    ssys = multibody(model)

    @test ssys isa System

    # The system should have states for η, η̇, q_flex, q̇_flex
    # n_modes * 2 + (n_boundaries - 1) * 6 * 2 = 2*2 + 1*6*2 = 16 DOFs before simplification

    # Create and solve ODE problem with small initial deformation
    prob = ODEProblem(ssys, [
        cb.η[1] => 0.01,  # Small initial modal amplitude
        cb.η[2] => 0.0,
    ], (0.0, 1.0))

    sol = solve(prob, Rodas5P())
    @test SciMLBase.successful_retcode(sol)

    # The modal amplitude should decay due to damping
    η1_initial = sol[cb.η[1]][1]
    η1_final = sol[cb.η[1]][end]
    @test abs(η1_final) < abs(η1_initial)  # Should decay
end

@testset "Matrix size validation" begin
    n_boundaries = 2
    n_modes = 3
    n_bdof = 6 * n_boundaries

    # Correct sizes
    M_BB = zeros(n_bdof, n_bdof)
    M_BI = zeros(n_bdof, n_modes)
    M_II = I(n_modes)
    K_BB = zeros(n_bdof, n_bdof)
    K_II = I(n_modes)
    ζ = zeros(n_modes)
    ω = ones(n_modes)

    # Should succeed
    @named cb_ok = CraigBampton(
        boundary_positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        M_BB = M_BB,
        M_BI = M_BI,
        M_II = M_II,
        K_BB = K_BB,
        K_II = K_II,
        ζ = ζ,
        ω = ω,
    )
    @test cb_ok isa System

    # Wrong M_BB size should fail
    @test_throws AssertionError CraigBampton(
        name = :cb_bad,
        boundary_positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        M_BB = zeros(10, 10),  # Wrong size
        M_BI = M_BI,
        M_II = M_II,
        K_BB = K_BB,
        K_II = K_II,
        ζ = ζ,
        ω = ω,
    )

    # Wrong ζ length should fail
    @test_throws AssertionError CraigBampton(
        name = :cb_bad2,
        boundary_positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        M_BB = M_BB,
        M_BI = M_BI,
        M_II = M_II,
        K_BB = K_BB,
        K_II = K_II,
        ζ = [0.01, 0.02],  # Wrong length (2 instead of 3)
        ω = ω,
    )
end

@testset "Three boundary component" begin
    # Test with 3 boundaries
    n_boundaries = 3
    n_modes = 2
    n_bdof = 6 * n_boundaries

    M_BB = 0.1 * I(n_bdof)
    M_BI = zeros(n_bdof, n_modes)
    M_II = I(n_modes)
    K_BB = 100.0 * I(n_bdof)
    K_II = Diagonal([100.0, 400.0])
    ζ = [0.02, 0.02]
    ω = [10.0, 20.0]

    @named cb3 = CraigBampton(
        boundary_positions = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]],
        M_BB = M_BB,
        M_BI = M_BI,
        M_II = M_II,
        K_BB = K_BB,
        K_II = K_II,
        ζ = ζ,
        ω = ω,
    )

    @test cb3 isa System
    @test hasproperty(cb3, :frame_a)
    @test hasproperty(cb3, :frame_b)
    @test hasproperty(cb3, :frame_c)
end

# =========================================================================
# Analytical Solution Tests
# =========================================================================
# These tests verify the modal dynamics against known analytical solutions
# for the damped harmonic oscillator: η̈ + 2ζω*η̇ + ω²*η = 0

"""
Helper function to create a minimal CB component for analytical testing.
Uses M_BI = 0 to decouple modal dynamics from boundary dynamics.
"""
function create_decoupled_cb(; n_modes, ω_modes, ζ_modes, name)
    n_boundaries = 2
    n_bdof = 6 * n_boundaries

    # Decoupled system: M_BI = 0 means modes evolve independently
    M_BB = 0.1 * I(n_bdof)
    M_BI = zeros(n_bdof, n_modes)
    M_II = Matrix(1.0 * I(n_modes))  # Mass-normalized modes
    K_BB = 100.0 * I(n_bdof)
    K_II = Diagonal(ω_modes.^2)

    CraigBampton(;
        name,
        boundary_positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        M_BB = M_BB,
        M_BI = M_BI,
        M_II = M_II,
        K_BB = K_BB,
        K_II = K_II,
        ζ = ζ_modes,
        ω = ω_modes,
    )
end

"""
Analytical solution for damped harmonic oscillator with IC: η(0)=η₀, η̇(0)=0
"""
function analytical_damped_oscillator(t, η₀, ω, ζ)
    if ζ < 1  # Underdamped
        ωd = ω * sqrt(1 - ζ^2)
        return η₀ * exp(-ζ * ω * t) * (cos(ωd * t) + (ζ / sqrt(1 - ζ^2)) * sin(ωd * t))
    elseif ζ == 1  # Critically damped
        return η₀ * exp(-ω * t) * (1 + ω * t)
    else  # Overdamped
        r1 = -ω * (ζ + sqrt(ζ^2 - 1))
        r2 = -ω * (ζ - sqrt(ζ^2 - 1))
        A = η₀ * r2 / (r2 - r1)
        B = -η₀ * r1 / (r2 - r1)
        return A * exp(r1 * t) + B * exp(r2 * t)
    end
end

"""
Analytical solution for damped oscillator with IC: η(0)=0, η̇(0)=v₀
"""
function analytical_velocity_ic(t, v₀, ω, ζ)
    ωd = ω * sqrt(1 - ζ^2)
    return (v₀ / ωd) * exp(-ζ * ω * t) * sin(ωd * t)
end

@testset "Undamped single mode - frequency verification" begin
    # Test that undamped oscillator has correct frequency
    # η(t) = η₀ * cos(ω*t), period T = 2π/ω
    ω_test = 10.0  # rad/s
    ζ_test = 1e-8  # Essentially undamped (can't use exactly 0)
    η₀ = 0.1

    cb = create_decoupled_cb(n_modes=1, ω_modes=[ω_test], ζ_modes=[ζ_test], name=:cb_undamped)

    eqs = [connect(world.frame_b, cb.frame_a)]
    @named model = System(eqs, t; systems=[world, cb])
    ssys = multibody(model)

    # Simulate for multiple periods
    T_period = 2π / ω_test
    t_end = 5 * T_period
    prob = ODEProblem(ssys, [cb.η[1] => η₀], (0.0, t_end))
    sol = solve(prob, Rodas5P(); saveat=T_period/100)

    @test SciMLBase.successful_retcode(sol)

    # Compare to analytical solution at sample times
    for (i, ti) in enumerate(sol.t[1:min(100, end)])
        η_numerical = sol[cb.η[1]][i]
        η_analytical = η₀ * cos(ω_test * ti)
        @test isapprox(η_numerical, η_analytical, rtol=0.01, atol=0.01)
    end

    # Verify period by checking zero crossings
    η_vals = sol[cb.η[1]]
    # Find time of first positive-to-negative crossing (should be at T/4)
    for i in 2:length(η_vals)-1
        if η_vals[i] > 0 && η_vals[i+1] < 0
            # Linear interpolation for zero crossing
            t_cross = sol.t[i] + (0 - η_vals[i]) / (η_vals[i+1] - η_vals[i]) * (sol.t[i+1] - sol.t[i])
            # First zero crossing of cos(ωt) is at T/4
            expected_cross = T_period / 4
            @test isapprox(t_cross, expected_cross, rtol=0.02, atol=0.01)
            break
        end
    end
end

@testset "Damped oscillator - exponential envelope" begin
    # Verify that envelope decays as exp(-ζ*ω*t)
    ω_test = 10.0
    ζ_test = 0.05  # 5% damping
    η₀ = 0.1

    cb = create_decoupled_cb(n_modes=1, ω_modes=[ω_test], ζ_modes=[ζ_test], name=:cb_damped_env)

    eqs = [connect(world.frame_b, cb.frame_a)]
    @named model = System(eqs, t; systems=[world, cb])
    ssys = multibody(model)

    t_end = 3.0  # Long enough to see significant decay
    prob = ODEProblem(ssys, [cb.η[1] => η₀], (0.0, t_end))
    sol = solve(prob, Rodas5P(); saveat=0.01)

    @test SciMLBase.successful_retcode(sol)

    # Check envelope at specific times
    test_times = [0.5, 1.0, 1.5, 2.0, 2.5]
    for t_check in test_times
        idx = findfirst(t -> t >= t_check, sol.t)
        if idx !== nothing
            η_numerical = sol[cb.η[1]][idx]
            η_analytical = analytical_damped_oscillator(sol.t[idx], η₀, ω_test, ζ_test)
            @test isapprox(η_numerical, η_analytical, rtol=0.05)
        end
    end

    # Verify overall decay rate by comparing early and late values
    η_early = maximum(abs.(sol[cb.η[1]][1:50]))
    η_late = maximum(abs.(sol[cb.η[1]][end-50:end]))
    expected_decay = exp(-ζ_test * ω_test * (sol.t[end-25] - sol.t[25]))
    actual_decay = η_late / η_early
    @test actual_decay < 1.0  # Should definitely decay
    @test isapprox(actual_decay, expected_decay, rtol=0.3)  # Within 30% of expected
end

@testset "Logarithmic decrement" begin
    # δ = ln(η(t)/η(t+Td)) = 2πζ/√(1-ζ²)
    ω_test = 8.0  # Lower frequency for clearer peaks
    ζ_test = 0.08  # 8% damping (larger for clearer decay)
    η₀ = 0.1

    cb = create_decoupled_cb(n_modes=1, ω_modes=[ω_test], ζ_modes=[ζ_test], name=:cb_logdec)

    eqs = [connect(world.frame_b, cb.frame_a)]
    @named model = System(eqs, t; systems=[world, cb])
    ssys = multibody(model)

    ωd = ω_test * sqrt(1 - ζ_test^2)
    Td = 2π / ωd  # Damped period
    t_end = 6 * Td

    prob = ODEProblem(ssys, [cb.η[1] => η₀], (0.0, t_end))
    sol = solve(prob, Rodas5P(); saveat=Td/50)

    @test SciMLBase.successful_retcode(sol)

    # Find peaks (local maxima)
    η_vals = sol[cb.η[1]]
    peaks = Float64[]
    peak_times = Float64[]
    for i in 2:length(η_vals)-1
        if η_vals[i] > η_vals[i-1] && η_vals[i] > η_vals[i+1] && η_vals[i] > 0.01 * η₀
            push!(peaks, η_vals[i])
            push!(peak_times, sol.t[i])
        end
    end

    # Need at least 2 peaks to compute log decrement
    if length(peaks) >= 2
        # Compute logarithmic decrement from first two peaks
        δ_measured = log(peaks[1] / peaks[2])
        δ_expected = 2π * ζ_test / sqrt(1 - ζ_test^2)
        @test isapprox(δ_measured, δ_expected, rtol=0.15)
    end
end

@testset "Damped natural frequency" begin
    # Verify damped frequency ωd = ω*√(1-ζ²)
    ω_test = 12.0
    ζ_test = 0.1  # 10% damping
    η₀ = 0.1

    cb = create_decoupled_cb(n_modes=1, ω_modes=[ω_test], ζ_modes=[ζ_test], name=:cb_damped_freq)

    eqs = [connect(world.frame_b, cb.frame_a)]
    @named model = System(eqs, t; systems=[world, cb])
    ssys = multibody(model)

    ωd_expected = ω_test * sqrt(1 - ζ_test^2)
    Td_expected = 2π / ωd_expected
    t_end = 5 * Td_expected

    prob = ODEProblem(ssys, [cb.η[1] => η₀], (0.0, t_end))
    sol = solve(prob, Rodas5P(); saveat=Td_expected/100)

    @test SciMLBase.successful_retcode(sol)

    # Find zero crossings (positive to negative) to measure period
    η_vals = sol[cb.η[1]]
    zero_crossings = Float64[]
    for i in 1:length(η_vals)-1
        if η_vals[i] > 0 && η_vals[i+1] <= 0
            # Linear interpolation
            t_cross = sol.t[i] - η_vals[i] * (sol.t[i+1] - sol.t[i]) / (η_vals[i+1] - η_vals[i])
            push!(zero_crossings, t_cross)
        end
    end

    # Period is time between consecutive positive-to-negative crossings
    if length(zero_crossings) >= 2
        Td_measured = zero_crossings[2] - zero_crossings[1]  # One full period
        @test isapprox(Td_measured, Td_expected, rtol=0.02)
    end
end

@testset "Multi-mode independence" begin
    # Two modes with different frequencies should oscillate independently
    ω1 = 10.0
    ω2 = 25.0
    ζ1 = 0.02
    ζ2 = 0.03
    η₀_1 = 0.1
    η₀_2 = 0.05

    cb = create_decoupled_cb(n_modes=2, ω_modes=[ω1, ω2], ζ_modes=[ζ1, ζ2], name=:cb_multimode)

    eqs = [connect(world.frame_b, cb.frame_a)]
    @named model = System(eqs, t; systems=[world, cb])
    ssys = multibody(model)

    t_end = 2.0
    prob = ODEProblem(ssys, [
        cb.η[1] => η₀_1,
        cb.η[2] => η₀_2,
    ], (0.0, t_end))
    sol = solve(prob, Rodas5P(); saveat=0.005)

    @test SciMLBase.successful_retcode(sol)

    # Each mode should match its own analytical solution
    for (i, ti) in enumerate(sol.t)
        # Mode 1
        η1_numerical = sol[cb.η[1]][i]
        η1_analytical = analytical_damped_oscillator(ti, η₀_1, ω1, ζ1)
        @test isapprox(η1_numerical, η1_analytical, rtol=0.05, atol=1e-4)

        # Mode 2
        η2_numerical = sol[cb.η[2]][i]
        η2_analytical = analytical_damped_oscillator(ti, η₀_2, ω2, ζ2)
        @test isapprox(η2_numerical, η2_analytical, rtol=0.05, atol=1e-4)
    end
end

@testset "Initial velocity response" begin
    # IC: η(0) = 0, η̇(0) = v₀
    # Solution: η(t) = (v₀/ωd) * exp(-ζ*ω*t) * sin(ωd*t)
    ω_test = 10.0
    ζ_test = 0.05
    v₀ = 1.0  # Initial velocity

    cb = create_decoupled_cb(n_modes=1, ω_modes=[ω_test], ζ_modes=[ζ_test], name=:cb_vel_ic)

    eqs = [connect(world.frame_b, cb.frame_a)]
    @named model = System(eqs, t; systems=[world, cb])
    ssys = multibody(model)

    ωd = ω_test * sqrt(1 - ζ_test^2)
    Td = 2π / ωd
    t_end = 3 * Td

    prob = ODEProblem(ssys, [
        cb.η[1] => 0.0,
        cb.η̇[1] => v₀,
    ], (0.0, t_end))
    sol = solve(prob, Rodas5P(); saveat=0.01)

    @test SciMLBase.successful_retcode(sol)

    # Compare to analytical solution
    for (i, ti) in enumerate(sol.t)
        η_numerical = sol[cb.η[1]][i]
        η_analytical = analytical_velocity_ic(ti, v₀, ω_test, ζ_test)
        # Use absolute tolerance near zero, relative elsewhere
        if abs(η_analytical) > 0.01
            @test isapprox(η_numerical, η_analytical, rtol=0.05)
        else
            @test isapprox(η_numerical, η_analytical, atol=0.01)
        end
    end

    # Verify first peak location (should be at t = arctan(ωd/(ζω)) / ωd for this IC)
    # For small damping: peak ≈ Td/4
    η_vals = sol[cb.η[1]]
    peak_idx = argmax(η_vals)
    t_peak = sol.t[peak_idx]
    t_peak_expected = atan(ωd / (ζ_test * ω_test)) / ωd
    @test isapprox(t_peak, t_peak_expected, rtol=0.1)
end

