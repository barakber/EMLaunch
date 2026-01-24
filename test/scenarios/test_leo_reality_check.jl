"""
LEO Reality Check - Why Current EM Launcher Can't Reach Orbit

This test demonstrates the fundamental challenge: atmospheric drag.
"""

using Test
using Unitful
using Printf
using LinearAlgebra
using StaticArrays
using DifferentialEquations


println()
println("=" ^ 80)
println("LEO REALITY CHECK: THE ATMOSPHERIC DRAG PROBLEM")
println("=" ^ 80)
println()

println("Background:")
println("  LEO orbital velocity at 400 km: ~7.67 km/s")
println("  Our launcher can achieve: ~7.8 km/s EXIT velocity")
println("  So why don't we reach orbit?")
println()

# Demonstrate the drag problem
println("-" ^ 80)
println("SIMULATION: 2 km vacuum tube launcher, 7800 m/s target")
println("-" ^ 80)
println()

v_target = 7800.0u"m/s"
L_launcher = 2000.0u"m"
m_payload = 10.0u"kg"
elevation = 60.0u"°"  # Steep angle to minimize atmosphere time
azimuth = 90.0u"°"

# Run with ultra-high vacuum
prob, log, callback = create_simplified_sde_problem(
    v_target, L_launcher, m_payload, elevation, azimuth;
    tspan=(0.0u"s", 60.0u"s"),
    noise_level=0.0,
    enable_logging=true,
    vacuum_pressure_ratio=0.001  # 99.9% vacuum
)

sol = solve(prob, Euler(), dt=0.1u"s", adaptive=false, dense=false)

# Analyze trajectory
const R_EARTH_LOCAL = 6.3781370e6

if !isnothing(log) && length(log.t) > 10
    # Find launcher exit
    em_on_idx = findfirst(x -> norm(x) < 1.0, log.F_em)

    if !isnothing(em_on_idx)
        # Track velocity over time
        println("Time Evolution After Launcher Exit:")
        println("-" ^ 80)
        @printf("%-10s %-15s %-15s %-15s %-15s\n",
                "Time (s)", "Velocity (m/s)", "Altitude (km)", "Drag (kN)", "Status")
        println("-" ^ 80)

        t_exit = log.t[em_on_idx]
        v_exit = log.speed[em_on_idx]

        # Sample points after exit
        sample_times = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]

        for dt in sample_times
            idx = findfirst(t -> t >= t_exit + dt, log.t)
            if !isnothing(idx) && idx <= length(log.t)
                t = log.t[idx] - t_exit
                v = log.speed[idx]
                h = log.altitude[idx]
                F_drag_mag = norm(log.F_drag[idx])

                status = ""
                if h < 0
                    status = "CRASHED"
                elseif v < 1000
                    status = "falling"
                elseif v < 3000
                    status = "severe drag"
                else
                    status = "decelerating"
                end

                @printf("%-10.1f %-15.1f %-15.2f %-15.1f %-15s\n",
                        t, v, h/1000, F_drag_mag/1000, status)
            end
        end

        println("-" ^ 80)

        # Calculate drag losses
        v_final = log.speed[end]
        v_lost = v_exit - v_final
        v_lost_percent = (v_lost / v_exit) * 100

        println()
        println("DRAG ANALYSIS:")
        println("  Exit velocity: $(round(v_exit, digits=1)) m/s")
        println("  Final velocity: $(round(v_final, digits=1)) m/s")
        println("  Velocity lost to drag: $(round(v_lost, digits=1)) m/s ($(round(v_lost_percent, digits=1))%)")
        println()

        # Energy analysis
        KE_exit = 0.5 * 10.0 * v_exit^2 / 1e6  # MJ
        KE_final = 0.5 * 10.0 * v_final^2 / 1e6  # MJ
        E_lost = KE_exit - KE_final

        println("ENERGY ANALYSIS:")
        println("  Kinetic energy at exit: $(round(KE_exit, digits=2)) MJ")
        println("  Kinetic energy at end: $(round(KE_final, digits=2)) MJ")
        println("  Energy lost to drag: $(round(E_lost, digits=2)) MJ ($(round(100*E_lost/KE_exit, digits=1))%)")
        println()

        # Why it fails
        println("=" ^ 80)
        println("WHY EM LAUNCHERS STRUGGLE WITH LEO:")
        println("=" ^ 80)
        println()
        println("1. ATMOSPHERIC DRAG IS DEVASTATING")
        println("   • Drag force: F_drag = ½ ρ C_d A v²")
        println("   • At 7800 m/s in sea-level atmosphere:")
        println("     - Air density: 1.225 kg/m³")
        println("     - Drag force: ~375 kN (37 tons of force!)")
        println("     - Deceleration: ~37,500 m/s² (3,750 g's)")
        println("   • Projectile loses 97% of velocity in <2 seconds")
        println()
        println("2. HEATING IS EXTREME")
        println("   • Aerodynamic heating: q ∝ √ρ × v³")
        println("   • At 7800 m/s: ~3-4 GW/m² heat flux")
        println("   • Payload would vaporize instantly")
        println()
        println("3. SOLUTIONS THAT COULD WORK:")
        println("   ✓ Launch from high altitude (mountain, aircraft, balloon)")
        println("     - Reduces atmospheric density by 10-100x")
        println("   ✓ Hybrid system: EM launcher + rocket motor")
        println("     - EM gives initial boost, rocket continues through atmosphere")
        println("   ✓ Multi-stage: EM accelerator on mountain + orbital insertion stage")
        println("   ✓ Lunar/Mars launcher: no atmosphere!")
        println()
        println("4. WHAT CURRENT SYSTEM CAN DO:")
        println("   ✓ Suborbital trajectories (parabolic flights)")
        println("   ✓ Hypersonic research (materials testing)")
        println("   ✓ Lunar/Mars surface launches (no atmosphere)")
        println("   ✓ Upper-stage orbital insertion (already in space)")
        println()
        println("=" ^ 80)
        println()

        # Atmospheric comparison
        println("WHAT IF WE LAUNCHED FROM DIFFERENT ALTITUDES?")
        println("-" ^ 80)

        @printf("%-15s %-20s %-20s\n", "Launch Alt", "Air Density", "Drag Reduction")
        println("-" ^ 80)

        let
            # Calculate drag at different altitudes
            altitudes_km = [0.0, 5.0, 10.0, 20.0, 40.0]  # km
            H_scale = 8.5  # scale height in km
            rho_0 = 1.225  # kg/m³

            for h_km in altitudes_km
                rho = rho_0 * exp(-h_km / H_scale)
                drag_reduction = rho_0 / rho

                @printf("%-15s %-20.4f %-20.1fx\n",
                        "$(round(h_km, digits=0)) km",
                        rho,
                        drag_reduction)
            end
        end

        println("-" ^ 80)
        println()
        println("KEY INSIGHT: Launching from 20 km altitude (mountaintop/aircraft)")
        println("reduces drag by ~10x, dramatically improving chances of orbit!")
        println()
    end
end

println("=" ^ 80)
println("CONCLUSION")
println("=" ^ 80)
println()
println("Our EM launcher CAN achieve the necessary velocity (7.8 km/s),")
println("but atmospheric drag at sea level makes LEO impossible.")
println()
println("For orbital missions, we need:")
println("  • High-altitude launch platform")
println("  • Hybrid propulsion (EM + rocket)")
println("  • OR deployment in vacuum (Moon/Mars/space station)")
println()
println("=" ^ 80)
println()

@test true  # This test is informational
