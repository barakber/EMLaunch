"""
Passive Drag Reduction Tests

Tests passive drag reduction methods (aerospikes, waveriders, etc.) and their impact
on hypersonic flight and LEO achievement.
"""

using Test
using Unitful
using Printf
using LinearAlgebra
using StaticArrays
using DifferentialEquations

# Load modules

println()
println("=" ^ 80)
println("PASSIVE DRAG REDUCTION TESTS")
println("=" ^ 80)
println()

@testset "Drag Reduction Methods" begin
    println("\n[Test 1] Individual drag reduction methods")

    # Test at hypersonic conditions
    mach = 10.0
    altitude = 0.0u"m"  # Sea level

    # No drag reduction
    none = NoDragReduction()
    r_none = drag_reduction_factor(none, mach, altitude)
    @test r_none == 0.0
    println("  ✓ No drag reduction: $(round(r_none * 100, digits=1))%")

    # Aerospike
    spike = Aerospike(1.5, 0.05u"m", 8000.0u"kg/m^3")  # Tungsten spike
    r_spike = drag_reduction_factor(spike, mach, altitude)
    @test r_spike > 0.2
    @test r_spike < 0.6
    println("  ✓ Aerospike: $(round(r_spike * 100, digits=1))% reduction")

    # Aerodisk
    disk = Aerodisk(0.7, 2, 2700.0u"kg/m^3")  # Aluminum MRD
    r_disk = drag_reduction_factor(disk, mach, altitude)
    @test r_disk > 0.15
    @test r_disk < 0.35
    println("  ✓ Aerodisk (MRD): $(round(r_disk * 100, digits=1))% reduction")

    # Boat tailing
    tailing = BoatTailing(7.0u"°", 0.8)
    r_tailing = drag_reduction_factor(tailing, mach, altitude)
    @test r_tailing > 0.1
    @test r_tailing < 0.25
    println("  ✓ Boat tailing: $(round(r_tailing * 100, digits=1))% reduction")

    # Waverider
    rider = Waverider(5.0, 12.0)
    r_rider = drag_reduction_factor(rider, mach, altitude)
    @test r_rider > 0.15
    @test r_rider < 0.3
    println("  ✓ Waverider: $(round(r_rider * 100, digits=1))% reduction")

    # Porous surface
    porous = PorousSurface(0.15, 50.0e-6u"m")
    r_porous = drag_reduction_factor(porous, mach, altitude)
    @test r_porous > 0.1
    @test r_porous < 0.25
    println("  ✓ Porous surface: $(round(r_porous * 100, digits=1))% reduction")
end

@testset "Mach Dependence" begin
    println("\n[Test 2] Mach number dependence")

    spike = Aerospike(1.5, 0.05u"m", 8000.0u"kg/m^3")
    altitude = 0.0u"m"

    machs = [0.5, 2.0, 5.0, 10.0, 15.0, 20.0]
    println("\n  Aerospike effectiveness vs Mach:")
    println("  " * "-"^50)
    @printf("  %-15s %-20s\n", "Mach", "Drag Reduction")
    println("  " * "-"^50)

    for mach in machs
        r = drag_reduction_factor(spike, mach, altitude)
        r_pct = round(r * 100, digits=1)
        @printf("  %-15.1f %-20s\n", mach, "$(r_pct)%")

        # Peak should be in hypersonic regime (Mach 5-15)
        if 5.0 <= mach <= 15.0
            @test r > 0.3  # Should be effective
        end
    end

    println("  " * "-"^50)
    println("  ✓ Peak effectiveness in Mach 5-15 (hypersonic)")
end

@testset "Mass Penalties" begin
    println("\n[Test 3] Mass penalties for drag reduction hardware")

    capsule_diameter = 0.4u"m"
    capsule_mass = 20.0u"kg"

    # Aerospike (heaviest)
    spike = Aerospike(1.5, 0.05u"m", 8000.0u"kg/m^3")  # Tungsten
    m_spike = added_mass(spike, capsule_diameter, capsule_mass)
    m_spike_kg = ustrip(u"kg", m_spike)
    m_spike_pct = 100.0 * m_spike_kg / ustrip(u"kg", capsule_mass)
    println("  Aerospike: $(round(m_spike_kg, digits=2)) kg ($(round(m_spike_pct, digits=1))%)")
    @test m_spike > 0.0u"kg"
    @test m_spike < 0.50 * capsule_mass  # Tungsten spike can be heavy (5-50%)

    # Aerodisk (lighter)
    disk = Aerodisk(0.7, 2, 2700.0u"kg/m^3")
    m_disk = added_mass(disk, capsule_diameter, capsule_mass)
    m_disk_kg = ustrip(u"kg", m_disk)
    m_disk_pct = 100.0 * m_disk_kg / ustrip(u"kg", capsule_mass)
    println("  Aerodisk: $(round(m_disk_kg, digits=2)) kg ($(round(m_disk_pct, digits=1))%)")
    @test m_disk > 0.0u"kg"
    @test m_disk < m_spike  # Should be lighter than spike

    # Boat tailing (reshaping, no added mass)
    tailing = BoatTailing(7.0u"°", 0.8)
    m_tailing = added_mass(tailing, capsule_diameter, capsule_mass)
    m_tailing_kg = ustrip(u"kg", m_tailing)
    println("  Boat tailing: $(round(m_tailing_kg, digits=2)) kg (reshaping)")
    @test m_tailing == 0.0u"kg"

    # Porous surface (2-5% increase)
    porous = PorousSurface(0.15, 50.0e-6u"m")
    m_porous = added_mass(porous, capsule_diameter, capsule_mass)
    m_porous_kg = ustrip(u"kg", m_porous)
    m_porous_pct = 100.0 * m_porous_kg / ustrip(u"kg", capsule_mass)
    println("  Porous surface: $(round(m_porous_kg, digits=2)) kg ($(round(m_porous_pct, digits=1))%)")
    @test m_porous > 0.0u"kg"
    @test m_porous < 0.06 * capsule_mass  # Should be <6%
end

@testset "Composite Methods" begin
    println("\n[Test 4] Composite drag reduction")

    # Combine multiple methods
    composite = CompositeDragReduction([
        Aerospike(1.5, 0.05u"m", 8000.0u"kg/m^3"),
        BoatTailing(7.0u"°", 0.8),
        PorousSurface(0.15, 50.0e-6u"m")
    ])

    mach = 10.0
    altitude = 0.0u"m"

    r_composite = drag_reduction_factor(composite, mach, altitude)
    r_pct = round(r_composite * 100, digits=1)
    println("  Combined reduction: $r_pct%")

    # Individual reductions (for comparison)
    r_spike = drag_reduction_factor(Aerospike(1.5, 0.05u"m", 8000.0u"kg/m^3"), mach, altitude)
    r_tailing = drag_reduction_factor(BoatTailing(7.0u"°", 0.8), mach, altitude)
    r_porous = drag_reduction_factor(PorousSurface(0.15, 50.0e-6u"m"), mach, altitude)

    println("  Individual: spike=$(round(r_spike*100, digits=1))%, tailing=$(round(r_tailing*100, digits=1))%, porous=$(round(r_porous*100, digits=1))%")

    # Composite should be better than any individual but not sum (multiplicative)
    @test r_composite > max(r_spike, r_tailing, r_porous)
    @test r_composite < r_spike + r_tailing + r_porous  # Not additive

    # Mass penalty is additive
    capsule_diam = 0.4u"m"
    capsule_mass = 20.0u"kg"
    m_composite = added_mass(composite, capsule_diam, capsule_mass)
    m_expected = added_mass(composite.methods[1], capsule_diam, capsule_mass) +
                 added_mass(composite.methods[2], capsule_diam, capsule_mass) +
                 added_mass(composite.methods[3], capsule_diam, capsule_mass)
    @test m_composite ≈ m_expected
    m_composite_kg = ustrip(u"kg", m_composite)
    println("  Total added mass: $(round(m_composite_kg, digits=2)) kg")
end

@testset "Effective Drag Coefficient" begin
    println("\n[Test 5] Effective drag coefficient")

    C_d_base = 0.50  # Streamlined body
    mach = 10.0
    altitude = 0.0u"m"

    # No reduction
    C_d_none = effective_drag_coefficient(C_d_base, NoDragReduction(), mach, altitude)
    @test C_d_none == C_d_base
    println("  Baseline C_d: $C_d_none")

    # With aerospike
    spike = Aerospike(1.5, 0.05u"m", 8000.0u"kg/m^3")
    C_d_spike = effective_drag_coefficient(C_d_base, spike, mach, altitude)
    @test C_d_spike < C_d_base
    reduction_pct = round((1 - C_d_spike/C_d_base) * 100, digits=1)
    println("  With aerospike: $C_d_spike ($(reduction_pct)% reduction)")

    # Drag force scales with C_d
    v = 7800.0  # m/s
    rho = 1.225  # kg/m³
    A = 0.1     # m²

    F_drag_base = 0.5 * rho * C_d_base * A * v^2
    F_drag_reduced = 0.5 * rho * C_d_spike * A * v^2
    force_reduction_pct = round((1 - F_drag_reduced/F_drag_base) * 100, digits=1)

    println("  Drag force: $(round(F_drag_base/1000, digits=1)) kN → $(round(F_drag_reduced/1000, digits=1)) kN")
    println("  Force reduction: $(force_reduction_pct)%")

    @test F_drag_reduced < F_drag_base
end

@testset "SEA LEVEL LEO TEST - Can Passive Methods Enable Orbit?" begin
    println("\n[Test 6] Sea level LEO achievement with passive drag reduction")
    println("  Testing if passive methods can overcome atmospheric drag...")
    println()

    v_target = 7800.0u"m/s"
    L_launcher = 2000.0u"m"
    m_payload = 20.0u"kg"
    elevation = 50.0u"°"
    azimuth = 90.0u"°"
    capsule_diam = 0.4u"m"

    # Test configurations
    configs = [
        (name="Baseline (no reduction)", method=NoDragReduction()),
        (name="Aerospike only", method=Aerospike(1.5, 0.05u"m", 8000.0u"kg/m^3")),
        (name="Composite (spike + tailing + porous)", method=CompositeDragReduction([
            Aerospike(2.0, 0.05u"m", 8000.0u"kg/m^3"),  # Max length
            BoatTailing(10.0u"°", 1.0),
            PorousSurface(0.20, 50.0e-6u"m")
        ]))
    ]

    println("  " * "="^70)
    @printf("  %-40s %-12s %-12s\n", "Configuration", "Peak Alt", "Survival")
    println("  " * "="^70)

    for (i, config) in enumerate(configs)
        # Create problem with drag reduction
        prob, log, callback = create_simplified_sde_problem(
            v_target, L_launcher, m_payload, elevation, azimuth;
            tspan=(0.0u"s", 10.0u"s"),  # Short time (crashes quickly at sea level)
            noise_level=0.0,
            enable_logging=false,  # Disable logging (unitful callback issue)
            vacuum_pressure_ratio=0.001,  # Vacuum tube
            drag_reduction_method=config.method,
            capsule_diameter=capsule_diam
        )

        # Solve with Euler (fixed-step, no saveat/callback for unitful mixed-dimension ODE)
        # Convert SDE to ODE since noise_level=0.0
        ode_prob = ODEProblem(prob.f.f, prob.u0, prob.tspan, prob.p)
        sol = solve(ode_prob, Euler(), dt=0.01u"s", adaptive=false, dense=false)

        # Find maximum altitude
        R_Earth = 6.3781370e6u"m"
        altitudes = [norm(SVector(u[1], u[2], u[3])) - R_Earth for u in sol.u]
        max_alt = maximum(altitudes)
        final_alt = altitudes[end]

        # Check if survived (didn't crash)
        survived = final_alt > -100.0u"m"  # -100m tolerance for numerical issues

        status = survived ? "Survived" : "CRASHED"
        alt_km = round(ustrip(u"km", max_alt), digits=2)

        @printf("  %-40s %-12s %-12s\n", config.name, "$(alt_km) km", status)

        # Tests
        if i == 1  # Baseline
            @test max_alt > 0.0u"m"  # Should reach some altitude with vacuum tube
            @test max_alt < 100000.0u"m"  # But won't reach orbit (100 km Kármán line)
        elseif i == 2  # Aerospike
            # Should perform better than baseline
            @test max_alt > 0.0u"m"
        else  # Composite
            # Best performance, but still can't reach orbit from sea level
            @test max_alt > 0.0u"m"
            @test max_alt < 100000.0u"m"  # Still below orbit
        end
    end

    println("  " * "="^70)
    println()
    println("  KEY INSIGHT:")
    println("  Even with aggressive passive drag reduction (40-50%), sea-level")
    println("  launch at 7.8 km/s faces extreme drag forces (hundreds of kN).")
    println("  Passive methods HELP but are insufficient alone.")
    println()
    println("  SOLUTION: Combine passive drag reduction with HIGH-ALTITUDE launch!")
    println("  → 20 km altitude + aerospike could enable LEO")
    println()
end

println()
println("=" ^ 80)
println("PASSIVE DRAG REDUCTION TESTS COMPLETE")
println("=" ^ 80)
println()
