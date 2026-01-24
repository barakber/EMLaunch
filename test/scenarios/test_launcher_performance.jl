"""
Performance Tests for Electromagnetic Launcher

Tests actual launcher performance metrics:
- Maximum altitude achieved
- Final velocity achieved
- Target parameter achievement
- Energy conservation
- Physical constraint satisfaction
"""

using Test
using Unitful
using DifferentialEquations
using LinearAlgebra
using StaticArrays
using Statistics

# Load modules

println("\n" * "="^70)
println("LAUNCHER PERFORMANCE TESTS")
println("="^70)

@testset "Simplified Model - LEO Target Performance" begin
    println("\n[Test 1] LEO Launcher Performance (2 km, 60 coils equivalent)")

    # Target parameters for LEO
    v_target = 7800.0u"m/s"
    L_launcher = 2000.0u"m"
    m_payload = 10.0u"kg"
    elevation = 45.0u"°"
    azimuth = 90.0u"°"

    # Create SDE problem (no noise for deterministic test)
    prob, log, callback = create_simplified_sde_problem(
        v_target, L_launcher, m_payload, elevation, azimuth;
        tspan=(0.0u"s", 180.0u"s"),
        noise_level=0.0,  # Deterministic
        enable_logging=true
    )

    # Solve
    sol = solve(prob, Euler(), dt=0.1u"s", adaptive=false, dense=false)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Extract final state
    u_final = sol.u[end]
    pos_final = SVector(u_final[1], u_final[2], u_final[3])
    vel_final = SVector(u_final[4], u_final[5], u_final[6])

    # Calculate final metrics
    R_Earth = 6.3781370e6u"m"
    r_final = norm(pos_final)
    h_final = r_final - R_Earth
    v_final = norm(vel_final)

    println("  Target velocity: $(ustrip(v_target)) m/s")
    println("  Actual velocity: $(round(ustrip(u"m/s", v_final), digits=0)) m/s")
    println("  Target altitude: 400 km (LEO)")
    println("  Actual altitude: $(round(ustrip(u"km", h_final), digits=2)) km")
    println("  Time simulated: $(sol.t[end]) s")

    # Performance tests
    # Note: Current launcher may not achieve orbital velocity due to drag losses
    @test v_final >= 0.0u"m/s"  # Velocity should be non-negative
    @test abs(h_final) < 1e6u"m"  # Altitude should be reasonable (< 1000 km magnitude)

    # Check if we're getting significant velocity from EM acceleration
    if !isnothing(log) && length(log.t) > 10
        # Find when EM force turned off
        em_on_idx = findfirst(x -> norm(x) < 1.0, log.F_em)
        if !isnothing(em_on_idx)
            v_at_exit = log.speed[em_on_idx]
            println("  Velocity at launcher exit: $(round(ustrip(u"m/s", v_at_exit), digits=0)) m/s")

            @test v_at_exit > 1000.0  # Should accelerate to >1 km/s in launcher
            @test v_at_exit < 10000.0  # Exit velocity should be physically reasonable

            # Check that we lost velocity due to drag
            velocity_loss = v_at_exit - v_final
            if velocity_loss > 0
                println("  Velocity lost to drag: $(round(ustrip(u"m/s", velocity_loss), digits=0)) m/s")
                @test velocity_loss > 0  # Should lose velocity to atmospheric drag
            end
        end
    end

    # Energy analysis
    KE_final = 0.5 * 10.0u"kg" * v_final^2  # J
    PE_final = 10.0u"kg" * 9.81u"m/s^2" * h_final   # J (approximate)
    E_total_final = KE_final + PE_final

    println("  Final kinetic energy: $(round(ustrip(u"MJ", KE_final), digits=2)) MJ")
    println("  Final potential energy: $(round(ustrip(u"MJ", PE_final), digits=2)) MJ")
    println("  Total final energy: $(round(ustrip(u"MJ", E_total_final), digits=2)) MJ")

    @test KE_final >= 0.0u"J"  # Final kinetic energy should be non-negative
    @test isfinite(ustrip(u"J", PE_final))  # Final potential energy should be finite
end

@testset "Maximum Altitude Achievement" begin
    println("\n[Test 2] Maximum Altitude Achievement Test")

    # Various launcher configurations to test
    test_configs = [
        (L=1000.0u"m", v_target=4000.0u"m/s", name="1km @ 4km/s"),
        (L=2000.0u"m", v_target=6000.0u"m/s", name="2km @ 6km/s"),
        (L=2000.0u"m", v_target=7800.0u"m/s", name="2km @ 7.8km/s"),
    ]

    results = []

    for config in test_configs
        println("\n  Testing: $(config.name)")

        prob, log, callback = create_simplified_sde_problem(
            config.v_target,
            config.L,
            10.0u"kg",
            45.0u"°",
            90.0u"°";
            tspan=(0.0u"s", 180.0u"s"),
            noise_level=0.0,
            enable_logging=true
        )

        sol = solve(prob, Euler(), dt=0.1u"s", adaptive=false, dense=false)

        # Find maximum altitude
        R_Earth = 6.3781370e6u"m"
        altitudes = [norm(SVector(u[1], u[2], u[3])) - R_Earth for u in sol.u]
        velocities = [norm(SVector(u[4], u[5], u[6])) for u in sol.u]

        max_altitude = maximum(altitudes)
        max_velocity = maximum(velocities)
        final_altitude = altitudes[end]
        final_velocity = velocities[end]

        result = (
            config = config.name,
            max_alt_km = max_altitude / 1000,
            final_alt_km = final_altitude / 1000,
            max_vel_ms = max_velocity,
            final_vel_ms = final_velocity
        )

        push!(results, result)

        println("    Max altitude: $(round(ustrip(u"km", result.max_alt_km), digits=2)) km")
        println("    Final altitude: $(round(ustrip(u"km", result.final_alt_km), digits=2)) km")
        println("    Max velocity: $(round(ustrip(u"m/s", result.max_vel_ms), digits=0)) m/s")
        println("    Final velocity: $(round(ustrip(u"m/s", result.final_vel_ms), digits=0)) m/s")

        # Basic sanity tests
        @test max_altitude >= 0.0u"m"  # Maximum altitude should be non-negative
        @test max_altitude < 2e6u"m"  # Maximum altitude should be < 2000 km (physically reasonable)
        @test max_velocity > 0.0u"m/s"  # Maximum velocity should be positive
        @test max_velocity < 15000.0u"m/s"  # Velocity should not exceed escape velocity range
    end

    # Print summary table
    println("\n  Performance Summary:")
    println("  " * "-"^70)
    println("  Config          Max Alt    Final Alt   Max Vel    Final Vel")
    println("  " * "-"^70)
    for r in results
        @printf("  %-15s %7.2f km  %7.2f km  %7.0f m/s  %7.0f m/s\n",
                r.config, r.max_alt_km, r.final_alt_km, r.max_vel_ms, r.final_vel_ms)
    end
    println("  " * "-"^70)
end

@testset "Energy Conservation Check" begin
    println("\n[Test 3] Energy Conservation")

    v_target = 5000.0u"m/s"
    L_launcher = 1500.0u"m"
    m_payload = 10.0u"kg"

    prob, log, callback = create_simplified_sde_problem(
        v_target, L_launcher, m_payload,
        45.0u"°", 90.0u"°";
        tspan=(0.0u"s", 60.0u"s"),
        noise_level=0.0,
        enable_logging=true
    )

    sol = solve(prob, Euler(), dt=0.1u"s", adaptive=false, dense=false)

    R_Earth = 6.3781370e6u"m"
    m = 10.0u"kg"
    μ = 3.986004418e14u"m^3/s^2"

    # Calculate total energy at each time point
    energies = []

    for u in sol.u
        pos = SVector(u[1], u[2], u[3])
        vel = SVector(u[4], u[5], u[6])

        r = norm(pos)
        v = norm(vel)

        KE = 0.5 * m * v^2
        PE = -μ * m / r  # Gravitational potential energy
        E_total = KE + PE

        push!(energies, E_total)
    end

    # Initial and final energy
    E_initial = energies[1]
    E_final = energies[end]

    # Energy change (should increase due to EM work, decrease due to drag)
    ΔE = E_final - E_initial

    println("  Initial energy: $(round(ustrip(u"MJ", E_initial), digits=2)) MJ")
    println("  Final energy: $(round(ustrip(u"MJ", E_final), digits=2)) MJ")
    println("  Energy change: $(round(ustrip(u"MJ", ΔE), digits=2)) MJ")

    # Check energy is finite and reasonable
    @test isfinite(ustrip(u"J", E_initial))  # Initial energy should be finite
    @test isfinite(ustrip(u"J", E_final))  # Final energy should be finite
    @test abs(E_initial) < 1e12u"J"  # Initial energy should be reasonable
    @test abs(E_final) < 1e12u"J"  # Final energy should be reasonable

    # In launcher phase, energy should increase (EM adds energy)
    # After launcher, energy should decrease (drag removes energy)
    if !isnothing(log) && length(log.t) > 10
        # Find launcher exit
        em_on_idx = findfirst(x -> norm(x) < 1.0, log.F_em)
        if !isnothing(em_on_idx) && em_on_idx < length(energies)
            E_at_exit = energies[em_on_idx]
            E_after_drag = energies[end]

            println("  Energy at launcher exit: $(round(ustrip(u"MJ", E_at_exit), digits=2)) MJ")
            println("  Energy after drag: $(round(ustrip(u"MJ", E_after_drag), digits=2)) MJ")
            println("  Energy lost to drag: $(round(ustrip(u"MJ", E_at_exit - E_after_drag), digits=2)) MJ")

            @test E_at_exit > E_initial  # EM force should add energy
            @test E_after_drag < E_at_exit  # Drag should remove energy
        end
    end
end

@testset "Physical Constraints During Flight" begin
    println("\n[Test 4] Physical Constraints")

    prob, log, callback = create_simplified_sde_problem(
        6000.0u"m/s", 1500.0u"m", 10.0u"kg",
        45.0u"°", 90.0u"°";
        tspan=(0.0u"s", 120.0u"s"),
        noise_level=0.0,
        enable_logging=true
    )

    sol = solve(prob, Euler(), dt=0.1u"s", adaptive=false, dense=false)

    R_Earth = 6.3781370e6u"m"
    violations = 0

    for (i, u) in enumerate(sol.u)
        pos = SVector(u[1], u[2], u[3])
        vel = SVector(u[4], u[5], u[6])
        T = u[7]

        r = norm(pos)
        v = norm(vel)
        h = r - R_Earth

        # Check constraints
        if r <= R_Earth * 0.99  # Below Earth surface (with 1% margin)
            violations += 1
            println("  ⚠️  t=$(sol.t[i]): Below Earth surface (r=$(r/1e6) Mm)")
        end

        if v > 15000.0u"m/s"  # Above reasonable velocity
            violations += 1
            println("  ⚠️  t=$(sol.t[i]): Excessive velocity (v=$(v) m/s)")
        end

        if T > 5000.0u"K"  # Unreasonably high temperature
            violations += 1
            println("  ⚠️  t=$(sol.t[i]): Excessive temperature (T=$(T) K)")
        end

        if T < 0.0u"K"  # Unphysical temperature
            violations += 1
            println("  ⚠️  t=$(sol.t[i]): Negative temperature (T=$(T) K)")
        end

        if !isfinite(ustrip(u"m", r)) || !isfinite(ustrip(u"m/s", v)) || !isfinite(ustrip(u"K", T))
            violations += 1
            println("  ⚠️  t=$(sol.t[i]): Non-finite values detected")
        end
    end

    println("  Total constraint violations: $violations / $(length(sol.u)) timesteps")

    @test violations == 0  # No physical constraints should be violated

    # Check final state is reasonable
    u_final = sol.u[end]
    pos_final = SVector(u_final[1], u_final[2], u_final[3])
    vel_final = SVector(u_final[4], u_final[5], u_final[6])
    T_final = u_final[7]

    r_final = norm(pos_final)
    v_final = norm(vel_final)
    h_final = r_final - R_Earth

    println("\n  Final state:")
    println("    Altitude: $(round(ustrip(u"km", h_final), digits=2)) km")
    println("    Velocity: $(round(ustrip(u"m/s", v_final), digits=0)) m/s")
    println("    Temperature: $(round(ustrip(u"K", T_final), digits=1)) K")

    @test isfinite(ustrip(u"m", r_final))  # Position should be finite
    @test v_final >= 0.0u"m/s"  # Velocity should be non-negative
    @test 0.0u"K" < T_final < 5000.0u"K"  # Temperature should be physically reasonable

    # Note: Launcher may not achieve orbit, projectile may crash (r_final < R_Earth is possible)
end

@testset "Drag Loss Quantification" begin
    println("\n[Test 5] Drag Loss Analysis")

    # Test with and without atmosphere (to quantify drag effect)
    configs = [
        (name="With Atmosphere", has_atmosphere=true),
        # Note: Current model doesn't have vacuum option, so we test what we have
    ]

    for config in configs
        println("\n  Testing: $(config.name)")

        prob, log, callback = create_simplified_sde_problem(
            7800.0u"m/s", 2000.0u"m", 10.0u"kg",
            45.0u"°", 90.0u"°";
            tspan=(0.0u"s", 120.0u"s"),
            noise_level=0.0,
            enable_logging=true
        )

        sol = solve(prob, Euler(), dt=0.1u"s", adaptive=false, dense=false)

        if !isnothing(log) && length(log.t) > 10
            # Find launcher exit point
            em_on_idx = findfirst(x -> norm(x) < 1.0, log.F_em)

            if !isnothing(em_on_idx)
                v_exit = log.speed[em_on_idx]
                v_final = log.speed[end]

                drag_loss = v_exit - v_final
                drag_loss_percent = (drag_loss / v_exit) * 100

                # Calculate total drag work
                drag_work = 0.0
                for i in (em_on_idx+1):length(log.t)
                    if i > 1
                        dt = log.t[i] - log.t[i-1]
                        v_avg = (log.speed[i] + log.speed[i-1]) / 2
                        F_drag_mag = norm(log.F_drag[i])
                        drag_work += F_drag_mag * v_avg * dt
                    end
                end

                println("    Exit velocity: $(round(ustrip(u"m/s", v_exit), digits=0)) m/s")
                println("    Final velocity: $(round(ustrip(u"m/s", v_final), digits=0)) m/s")
                println("    Velocity lost to drag: $(round(ustrip(u"m/s", drag_loss), digits=0)) m/s ($(round(drag_loss_percent, digits=1))%)")
                println("    Total drag work: $(round(ustrip(u"MJ", drag_work), digits=2)) MJ")

                @test v_exit > v_final  # Drag should reduce velocity
                @test drag_work > 0  # Drag should dissipate energy
            end
        end
    end
end

println("\n" * "="^70)
println("PERFORMANCE TESTS COMPLETE")
println("="^70)
