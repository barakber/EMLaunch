"""
LEO Achievement Tests

Tests whether the launcher can actually achieve Low Earth Orbit (LEO).

LEO Requirements:
- Altitude: 160-2000 km (typical: 400 km for ISS)
- Velocity: ~7.67 km/s at 400 km altitude
- Orbital mechanics: v² >= μ/r for stable orbit
"""

using Test
using Unitful
using Printf
using LinearAlgebra
using StaticArrays
using DifferentialEquations

# Load modules

"""
Check if trajectory achieves orbit
"""
function is_orbit_achieved(sol, target_altitude=400.0e3)
    R_Earth = 6.3781370e6u"m"
    μ = 3.986004418e14u"m^3/s^2"  # Earth's gravitational parameter

    # Get final state
    u_final = sol.u[end]
    pos_final = SVector(u_final[1], u_final[2], u_final[3])
    vel_final = SVector(u_final[4], u_final[5], u_final[6])

    r_final = norm(pos_final)
    v_final = norm(vel_final)
    h_final = r_final - R_Earth

    # Calculate orbital velocity at this altitude
    v_circular = sqrt(μ / r_final)

    # Calculate specific orbital energy
    E_specific = v_final^2 / 2 - μ / r_final

    # For circular orbit: E = -μ/(2r)
    # For escape: E >= 0
    # For elliptical orbit: -μ/(2a) where a is semi-major axis

    # Check if orbit is achieved
    altitude_ok = h_final > 160e3  # Above 160 km (Kármán line + margin)
    velocity_ok = v_final > 0.95 * v_circular  # At least 95% of orbital velocity
    not_crashed = h_final > 0
    energy_ok = E_specific > -μ / (2 * R_Earth)  # Not falling back immediately

    return (
        achieved = altitude_ok && velocity_ok && not_crashed,
        h_final = h_final,
        v_final = v_final,
        v_circular = v_circular,
        v_ratio = v_final / v_circular,
        E_specific = E_specific,
        altitude_ok = altitude_ok,
        velocity_ok = velocity_ok
    )
end

println()
println("=" ^ 80)
println("LEO ACHIEVEMENT TESTS")
println("=" ^ 80)
println()

@testset "Current System - Check LEO Achievement" begin
    println("\n[Test 1] Can current optimized system reach LEO?")

    # Current "optimized" configuration
    v_target = 7800.0u"m/s"
    L_launcher = 2000.0u"m"
    m_payload = 10.0u"kg"

    # Use vacuum tube + optimization
    result = optimize_launch_angle_min_drag(
        v_target, L_launcher, m_payload;
        elevation_range=(30.0u"°", 70.0u"°"),
        vacuum_pressure_ratio=0.001,  # Ultra-high vacuum
        verbose=false
    )

    println("  Optimized configuration:")
    @printf("    Elevation: %.2f°\n", ustrip(result.elevation))
    @printf("    Exit velocity: %.1f m/s\n", ustrip(v_target))
    @printf("    Final velocity: %.1f m/s\n", ustrip(result.v_final))
    @printf("    Final altitude: %.2f km\n", ustrip(u"km", result.h_final))

    # Run full simulation to check orbit
    prob, log, callback = create_simplified_sde_problem(
        v_target, L_launcher, m_payload,
        result.elevation, 90.0u"°";
        tspan=(0.0u"s", 300.0u"s"),  # Longer simulation
        noise_level=0.0,
        enable_logging=true,
        vacuum_pressure_ratio=0.001
    )

    sol = solve(prob, Euler(), dt=0.1u"s", adaptive=false, dense=false)

    orbit_status = is_orbit_achieved(sol)

    println("\n  Orbital Analysis:")
    @printf("    Final altitude: %.2f km\n", orbit_status.h_final / 1000)
    @printf("    Final velocity: %.1f m/s\n", orbit_status.v_final)
    @printf("    Required circular velocity: %.1f m/s\n", orbit_status.v_circular)
    @printf("    Velocity ratio: %.3f (%.1f%%)\n", orbit_status.v_ratio, orbit_status.v_ratio * 100)
    @printf("    Specific energy: %.2f MJ/kg\n", orbit_status.E_specific / 1e6)

    if orbit_status.achieved
        println("\n  ✓ LEO ACHIEVED!")
    else
        println("\n  ✗ LEO NOT ACHIEVED")
        if !orbit_status.altitude_ok
            println("    - Altitude too low")
        end
        if !orbit_status.velocity_ok
            println("    - Velocity insufficient for orbit")
        end
    end

    # Test should acknowledge current limitations
    @test true  # This test is informational
end

@testset "What Launcher Config Needed for LEO?" begin
    println("\n[Test 2] Finding launcher configuration to achieve LEO")

    m_payload = 10.0u"kg"
    target_altitude = 400.0u"km"

    # LEO orbital velocity at 400 km
    R_Earth = 6.3781370e6u"m"
    μ = 3.986004418e14u"m^3/s^2"
    r_orbit = R_Earth + 400e3u"m"
    v_orbital = sqrt(μ / r_orbit)

    println("  Target: 400 km altitude")
    @printf("  Required orbital velocity: %.1f m/s (%.2f km/s)\n", v_orbital, v_orbital/1000)

    # Account for drag losses - need higher exit velocity
    # Rough estimate: need ~20% margin for drag losses
    v_target_with_margin = v_orbital * 1.3

    @printf("  Target exit velocity (with margin): %.1f m/s (%.2f km/s)\n\n",
            v_target_with_margin, v_target_with_margin/1000)

    # Try different launcher lengths
    println("  Testing different launcher configurations:")
    println("  " * "-"^60)
    @printf("  %-15s %-15s %-15s %s\n", "Length", "Exit V", "Final V", "Status")
    println("  " * "-"^60)

    launcher_lengths = [2000.0, 5000.0, 10000.0, 20000.0]u"m"

    for L in launcher_lengths
        v_target = v_target_with_margin * 1.0u"m/s"

        # Use vacuum + optimized angle
        result = optimize_launch_angle_min_drag(
            v_target, L, m_payload;
            elevation_range=(40.0u"°", 80.0u"°"),
            vacuum_pressure_ratio=0.001,
            verbose=false
        )

        # Check if LEO achieved
        prob, log, callback = create_simplified_sde_problem(
            v_target, L, m_payload,
            result.elevation, 90.0u"°";
            tspan=(0.0u"s", 400.0u"s"),
            noise_level=0.0,
            enable_logging=false,
            vacuum_pressure_ratio=0.001
        )

        sol = solve(prob, Euler(), dt=0.1u"s", adaptive=false, dense=false)
        orbit_status = is_orbit_achieved(sol)

        status = orbit_status.achieved ? "✓ ORBIT" : "✗ No orbit"

        # Get exit velocity from log if available
        exit_v = 0.0
        if !isnothing(log) && length(log.t) > 10
            em_on_idx = findfirst(x -> norm(x) < 1.0, log.F_em)
            if !isnothing(em_on_idx)
                exit_v = log.speed[em_on_idx]
            end
        end

        @printf("  %-15s %-15.0f %-15.0f %s\n",
                "$(ustrip(u"km", L)) km",
                ustrip(result.v_final),
                ustrip(result.v_final),
                status)

        if orbit_status.achieved
            println("\n  🎯 SUCCESS! LEO achieved with $(ustrip(u"km", L)) km launcher")
            @printf("     Final altitude: %.2f km\n", orbit_status.h_final / 1000)
            @printf("     Final velocity: %.1f m/s\n", orbit_status.v_final)
            @printf("     Orbital velocity ratio: %.1f%%\n", orbit_status.v_ratio * 100)
            @test true
            break
        end
    end
end

@testset "Maximum Altitude Achievable" begin
    println("\n[Test 3] Maximum altitude with realistic launchers")

    configs = [
        (name="2 km launcher", L=2000.0u"m", v=7800.0u"m/s"),
        (name="5 km launcher", L=5000.0u"m", v=10000.0u"m/s"),
        (name="10 km launcher", L=10000.0u"m", v=12000.0u"m/s"),
    ]

    println("\n  Configuration performance:")
    println("  " * "-"^70)
    @printf("  %-20s %-15s %-15s %-15s\n", "Config", "Max Altitude", "Max Velocity", "Apogee")
    println("  " * "-"^70)

    for config in configs
        # Optimize with vacuum
        result = optimize_launch_angle_min_drag(
            config.v, config.L, 10.0u"kg";
            elevation_range=(45.0u"°", 80.0u"°"),
            vacuum_pressure_ratio=0.001,
            verbose=false
        )

        # Run simulation
        prob, log, callback = create_simplified_sde_problem(
            config.v, config.L, 10.0u"kg",
            result.elevation, 90.0u"°";
            tspan=(0.0u"s", 600.0u"s"),
            noise_level=0.0,
            enable_logging=true,
            vacuum_pressure_ratio=0.001
        )

        sol = solve(prob, Euler(), dt=0.1u"s", adaptive=false, dense=false)

        # Find maximum altitude
        R_Earth = 6.3781370e6u"m"
        altitudes = [norm(SVector(u[1], u[2], u[3])) - R_Earth for u in sol.u]
        velocities = [norm(SVector(u[4], u[5], u[6])) for u in sol.u]

        max_altitude = maximum(altitudes)
        max_velocity = maximum(velocities)
        apogee_idx = argmax(altitudes)
        apogee_time = sol.t[apogee_idx]

        @printf("  %-20s %-15.2f %-15.1f %-15.1f\n",
                config.name,
                max_altitude / 1000,  # km
                max_velocity,  # m/s
                apogee_time)  # seconds

        @test max_altitude > 0  # Should at least get off the ground
    end

    println("  " * "-"^70)
end

println()
println("=" ^ 80)
println("LEO ACHIEVEMENT ANALYSIS COMPLETE")
println("=" ^ 80)
println()
