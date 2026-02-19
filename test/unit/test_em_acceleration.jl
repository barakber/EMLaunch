"""
Tests for em_acceleration.jl

Tests:
- Launcher configuration
- Electromagnetic force calculations
- Energy calculations
- Coil voltage control
- Physical reasonableness
"""


@testset "Launcher Configuration" begin
    # Create simple launcher
    launcher = create_uniform_launcher(
        length = 1000.0u"m",
        num_coils = 100,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 10.0u"F",
        voltage_per_coil = 5000.0u"V",
        gradient_per_coil = 0.002u"H/m"
    )

    @test launcher.length == 1000.0u"m"
    @test launcher.num_coils == 100
    @test length(launcher.coils) == 100

    # Check first coil
    coil1 = launcher.coils[1]
    @test coil1.position == 0.0u"m"
    @test coil1.inductance == 0.001u"H"
    @test coil1.resistance == 0.01u"Ω"
    @test coil1.capacitance == 10.0u"F"
    @test coil1.voltage == 5000.0u"V"

    # Check spacing
    coil2 = launcher.coils[2]
    spacing = coil2.position - coil1.position
    @test spacing ≈ 10.0u"m"  # 1000m / 100 coils

    # Check last coil position
    coil_last = launcher.coils[end]
    @test coil_last.position ≈ 990.0u"m"  # (99/100) * 1000m
end

@testset "Inductance Gradient" begin
    L = 0.001u"H"
    L_coil = 0.2u"m"

    dL_dx = inductance_gradient(L, L_coil)

    @test dimension(dL_dx) == dimension(1.0u"H/m")
    @test dL_dx > 0.0u"H/m"
    @test dL_dx ≈ 0.01u"H/m"  # 2 * 0.001 / 0.2
end

@testset "Electromagnetic Force" begin
    # Create a coil
    coil = CoilConfig(
        100.0u"m",          # position
        0.001u"H",          # inductance
        0.01u"Ω",           # resistance
        10.0u"F",           # capacitance
        5000.0u"V",         # voltage
        0.0u"s",            # trigger_time
        0.002u"H/m"         # gradient
    )

    # Current in coil
    I = 10000.0u"A"

    # Force calculation
    x = 99.0u"m"  # Just before coil
    F = em_force(x, I, coil)

    @test dimension(F) == dimension(1.0u"N")
    @test F > 0.0u"N"  # Should be attractive

    # Force proportional to I²
    I2 = 20000.0u"A"
    F2 = em_force(x, I2, coil)
    @test F2 ≈ 4.0 * F rtol=0.01  # I² scaling

    # Force should be zero far from coil
    F_far = em_force(50.0u"m", I, coil)
    @test F_far == 0.0u"N"

    F_past = em_force(105.0u"m", I, coil)
    @test F_past ≈ 0.0u"N" atol=1.0u"N"
end

@testset "Total EM Force from Multiple Coils" begin
    launcher = create_uniform_launcher(
        length = 100.0u"m",
        num_coils = 10,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 10.0u"F",
        voltage_per_coil = 5000.0u"V",
        gradient_per_coil = 0.002u"H/m"
    )

    # All coils at same current
    currents = fill(10000.0u"A", 10)

    # Position between coils
    x = 15.0u"m"
    F_total = total_em_force(x, currents, launcher)

    @test dimension(F_total) == dimension(1.0u"N")
    @test F_total > 0.0u"N"

    # Zero current → zero force
    currents_zero = fill(0.0u"A", 10)
    F_zero = total_em_force(x, currents_zero, launcher)
    @test F_zero == 0.0u"N"
end

@testset "Coil Voltage Control" begin
    coil = CoilConfig(
        100.0u"m",
        0.001u"H",
        0.01u"Ω",
        10.0u"F",
        5000.0u"V",
        0.5u"s",  # trigger at 0.5s
        0.002u"H/m"
    )

    # Before trigger: external source maintains capacitor charge
    V_before = coil_voltage(0.3u"s", coil)
    @test V_before == 5000.0u"V"

    # After trigger: source disconnected, capacitor discharges through RLC
    V_after = coil_voltage(0.6u"s", coil)
    @test V_after == 0.0u"V"

    # At trigger: discharge begins
    V_at = coil_voltage(0.5u"s", coil)
    @test V_at == 0.0u"V"
end

@testset "Energy Calculations" begin
    launcher = create_uniform_launcher(
        length = 2000.0u"m",
        num_coils = 8,  # Adjusted for 64% efficiency
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 10.0u"F",
        voltage_per_coil = 5000.0u"V",
        gradient_per_coil = 0.002u"H/m"
    )

    # Final velocity and payload mass
    v_final = 8000.0u"m/s"
    m_payload = 20.0u"kg"

    # Efficiency
    η = energy_efficiency(v_final, m_payload, launcher)

    @test 0.0 < η < 1.0  # Should be between 0 and 100%
    @test η ≈ 0.64 rtol=0.1  # Expected ~64% efficiency with 8 coils

    # Required energy
    E_req = required_capacitor_energy(v_final, m_payload, 0.60)

    @test dimension(E_req) == dimension(1.0u"J")
    @test E_req > 0.0u"J"

    # Kinetic energy
    KE = 0.5 * m_payload * v_final^2
    @test KE ≈ 640.0u"MJ" rtol=0.01

    # Required energy should be larger (due to efficiency)
    @test E_req > KE
    @test E_req ≈ 1067.0u"MJ" rtol=0.1
end

@testset "Capacitor Energy Storage" begin
    # Single coil
    C = 10.0u"F"
    V = 5000.0u"V"
    E = 0.5 * C * V^2

    @test dimension(E) == dimension(1.0u"J")
    @test E ≈ 125.0u"MJ" rtol=0.01

    # Total launcher energy
    launcher = create_uniform_launcher(
        length = 2000.0u"m",
        num_coils = 200,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 10.0u"F",
        voltage_per_coil = 5000.0u"V",
        gradient_per_coil = 0.002u"H/m"
    )

    E_total = sum(0.5 * coil.capacitance * coil.voltage^2
                  for coil in launcher.coils)

    @test E_total ≈ 200 * 125.0u"MJ" rtol=0.01
    @test E_total ≈ 25000.0u"MJ" rtol=0.01  # 25 GJ
end

@testset "Physical Reasonableness" begin
    # Force should be positive (accelerating)
    coil = CoilConfig(
        100.0u"m", 0.001u"H", 0.01u"Ω", 10.0u"F",
        5000.0u"V", 0.0u"s", 0.002u"H/m"
    )

    for I in [1000.0, 5000.0, 10000.0, 50000.0] .* 1.0u"A"
        F = em_force(99.0u"m", I, coil)
        @test F >= 0.0u"N"
    end

    # Higher current should give higher force
    F1 = em_force(99.0u"m", 10000.0u"A", coil)
    F2 = em_force(99.0u"m", 20000.0u"A", coil)
    @test F2 > F1

    # Force should scale as I²
    @test F2 / F1 ≈ 4.0 rtol=0.1
end

@testset "Coil Spacing Consistency" begin
    # Different coil counts should give consistent spacing
    for N in [50, 100, 200, 500]
        launcher = create_uniform_launcher(
            length = 1000.0u"m",
            num_coils = N,
            inductance_per_coil = 0.001u"H",
            resistance_per_coil = 0.01u"Ω",
            capacitance_per_coil = 10.0u"F",
            voltage_per_coil = 5000.0u"V",
            gradient_per_coil = 0.002u"H/m"
        )

        expected_spacing = 1000.0u"m" / N

        # Check spacing between consecutive coils
        for i in 1:(N-1)
            spacing = launcher.coils[i+1].position - launcher.coils[i].position
            @test spacing ≈ expected_spacing rtol=1e-6
        end
    end
end

@testset "Unit Consistency" begin
    # All electromagnetic quantities should have correct units
    coil = CoilConfig(
        100.0u"m", 0.001u"H", 0.01u"Ω", 10.0u"F",
        5000.0u"V", 0.0u"s", 0.002u"H/m"
    )

    @test dimension(coil.position) == dimension(1.0u"m")
    @test dimension(coil.inductance) == dimension(1.0u"H")
    @test dimension(coil.resistance) == dimension(1.0u"Ω")
    @test dimension(coil.capacitance) == dimension(1.0u"F")
    @test dimension(coil.voltage) == dimension(1.0u"V")
    @test dimension(coil.trigger_time) == dimension(1.0u"s")
    @test dimension(coil.gradient) == dimension(1.0u"H/m")

    # Force calculation units
    F = em_force(99.0u"m", 10000.0u"A", coil)
    @test dimension(F) == dimension(1.0u"N")

    # Energy units
    E = 0.5 * coil.capacitance * coil.voltage^2
    @test dimension(E) == dimension(1.0u"J")
end

# =============================================================================
# COMPREHENSIVE COIL VOLTAGE TESTS
# These tests guard against the critical coil_voltage inversion bug.
# Physical invariant: V=V0 maintains capacitor charge (switch open),
#                     V=0  allows RLC discharge (switch closed).
# =============================================================================

@testset "Coil Voltage — Semantic Invariant" begin
    # The CRITICAL invariant: triggered coils must return 0V to allow
    # capacitor discharge, untriggered must return V0 to maintain charge.
    for V0 in [100.0, 2000.0, 5000.0, 18000.0] .* 1.0u"V"
        coil = CoilConfig(50.0u"m", 0.001u"H", 0.01u"Ω", 1.0u"F", V0, 0.5u"s", 0.002u"H/m")

        # Time-based: before trigger = V0 (maintain), after trigger = 0 (discharge)
        @test coil_voltage(0.0u"s", coil) == V0    # Before trigger
        @test coil_voltage(0.5u"s", coil) == 0.0u"V"  # At trigger
        @test coil_voltage(1.0u"s", coil) == 0.0u"V"  # After trigger

        # Position-based (unitful): within range = 0 (discharge), outside = V0
        @test coil_voltage(50.0u"m", coil) == 0.0u"V"  # At coil (triggered)
        @test coil_voltage(300.0u"m", coil) == V0       # Far away (not triggered)

        # Position-based (dimensionless): same invariant
        @test coil_voltage(50.0, coil) == 0.0u"V"   # At coil (triggered)
        @test coil_voltage(300.0, coil) == V0        # Far away (not triggered)
    end
end

@testset "Coil Voltage — Position-Based Trigger Boundaries" begin
    # Coil at 200m, trigger_distance = 100m → active range [100m, 300m]
    coil = CoilConfig(200.0u"m", 0.001u"H", 0.01u"Ω", 1.0u"F", 5000.0u"V", 0.0u"s", 0.002u"H/m")

    # Unitful: outside trigger range → V0 (maintain charge)
    @test coil_voltage(0.0u"m", coil) == 5000.0u"V"    # Way before
    @test coil_voltage(99.0u"m", coil) == 5000.0u"V"    # Just before range
    @test coil_voltage(301.0u"m", coil) == 5000.0u"V"   # Just past range
    @test coil_voltage(500.0u"m", coil) == 5000.0u"V"   # Way past

    # Unitful: inside trigger range → 0V (allow discharge)
    @test coil_voltage(100.0u"m", coil) == 0.0u"V"   # At lower boundary
    @test coil_voltage(150.0u"m", coil) == 0.0u"V"   # Approaching coil
    @test coil_voltage(200.0u"m", coil) == 0.0u"V"   # At coil center
    @test coil_voltage(250.0u"m", coil) == 0.0u"V"   # Past coil
    @test coil_voltage(300.0u"m", coil) == 0.0u"V"   # At upper boundary

    # Dimensionless: same boundaries
    @test coil_voltage(99.0, coil) == 5000.0u"V"
    @test coil_voltage(100.0, coil) == 0.0u"V"
    @test coil_voltage(200.0, coil) == 0.0u"V"
    @test coil_voltage(300.0, coil) == 0.0u"V"
    @test coil_voltage(301.0, coil) == 5000.0u"V"
end

@testset "Coil Voltage — Time-Based Trigger Edge Cases" begin
    # Multiple trigger times
    for trigger_time in [0.0u"s", 0.1u"s", 1.0u"s", 10.0u"s"]
        coil = CoilConfig(0.0u"m", 0.001u"H", 0.01u"Ω", 1.0u"F", 3000.0u"V", trigger_time, 0.002u"H/m")

        if trigger_time > 0.0u"s"
            # Before trigger: maintaining charge
            @test coil_voltage(trigger_time - 0.001u"s", coil) == 3000.0u"V"
        end
        # At trigger: discharge begins
        @test coil_voltage(trigger_time, coil) == 0.0u"V"
        # After trigger: discharge continues
        @test coil_voltage(trigger_time + 1.0u"s", coil) == 0.0u"V"
    end
end

@testset "Coil Voltage — Sequential Triggering Along Launcher" begin
    # Simulate a projectile moving through a 5-coil launcher
    # Coils at 0m, 10m, 20m, 30m, 40m; trigger_distance = 100m
    coils = [CoilConfig(i * 10.0u"m", 0.001u"H", 0.01u"Ω", 1.0u"F", 2000.0u"V", 0.0u"s", 0.003u"H/m") for i in 0:4]

    # At x=0: all coils within 100m trigger range (positions 0-40m)
    for coil in coils
        @test coil_voltage(0.0u"m", coil) == 0.0u"V"  # All triggered
    end

    # At x=200m: no coils within 100m range (all coils at 0-40m)
    for coil in coils
        @test coil_voltage(200.0u"m", coil) == 2000.0u"V"  # None triggered
    end
end

@testset "Coil Voltage — Coil at Position Zero" begin
    # Special case: first coil at position 0
    coil = CoilConfig(0.0u"m", 0.001u"H", 0.01u"Ω", 0.2u"F", 2000.0u"V", 0.0u"s", 0.003u"H/m")

    # Projectile at launcher start (x=0): coil triggered
    @test coil_voltage(0.0u"m", coil) == 0.0u"V"
    @test coil_voltage(0.0, coil) == 0.0u"V"

    # Projectile far away: coil not triggered
    @test coil_voltage(200.0u"m", coil) == 2000.0u"V"
    @test coil_voltage(200.0, coil) == 2000.0u"V"
end

# =============================================================================
# RLC CIRCUIT DYNAMICS TESTS
# Verifies that the circuit model produces correct current when triggered.
# =============================================================================

@testset "RLC Dynamics — Triggered Coil Produces Current" begin
    # When a coil is triggered, V_applied=0 and Q/C=V0, so
    # dI/dt = (0 - V0 - 0) / L = -V0/L (non-zero!)
    for (V0, L, C) in [
        (2000.0, 0.001, 0.2),   # Moon launcher
        (5000.0, 0.001, 10.0),  # Standard launcher
        (18000.0, 0.001, 2.5),  # High-power launcher
    ]
        coil = CoilConfig(50.0u"m", L*u"H", 0.01u"Ω", C*u"F", V0*u"V", 0.0u"s", 0.003u"H/m")
        Q0 = coil.capacitance * coil.voltage   # Fully charged
        I0 = 0.0u"A"

        # Triggered: V_applied = 0
        V_triggered = coil_voltage(50.0u"m", coil)
        @test V_triggered == 0.0u"V"

        dI_dt = (V_triggered - Q0 / coil.capacitance - coil.resistance * I0) / coil.inductance
        @test abs(ustrip(u"A/s", dI_dt)) > 1e4  # Significant current change
        @test ustrip(u"A/s", dI_dt) ≈ -V0 / L   # Exact: -V0/L
    end
end

@testset "RLC Dynamics — Untriggered Coil Holds Charge" begin
    # When not triggered, V_applied=V0 and Q/C=V0, so dI/dt = 0
    coil = CoilConfig(50.0u"m", 0.001u"H", 0.01u"Ω", 0.2u"F", 2000.0u"V", 0.0u"s", 0.003u"H/m")
    Q0 = coil.capacitance * coil.voltage
    I0 = 0.0u"A"

    # Not triggered: V_applied = V0
    V_untriggered = coil_voltage(500.0u"m", coil)  # Far from coil
    @test V_untriggered == 2000.0u"V"

    dI_dt = (V_untriggered - Q0 / coil.capacitance - coil.resistance * I0) / coil.inductance
    @test ustrip(u"A/s", dI_dt) ≈ 0.0 atol=1e-10  # No current change
end

@testset "RLC Dynamics — Current Produces EM Force" begin
    # After discharge begins, current builds → F = 0.5 * dL/dx * I² > 0
    coil = CoilConfig(50.0u"m", 0.001u"H", 0.01u"Ω", 0.2u"F", 2000.0u"V", 0.0u"s", 0.003u"H/m")

    # Simulate one Euler step of RLC discharge
    Q = coil.capacitance * coil.voltage  # 400 C
    I = 0.0u"A"
    dt = 1e-6u"s"  # 1 μs step

    V = coil_voltage(50.0u"m", coil)  # Triggered: 0V
    dI_dt = (V - Q / coil.capacitance - coil.resistance * I) / coil.inductance
    I_new = I + dI_dt * dt
    Q_new = Q + I * dt

    @test abs(I_new) > 0.0u"A"  # Current has begun flowing

    # After a few more steps, force should be significant
    for _ in 1:100
        V = coil_voltage(50.0u"m", coil)
        dI_dt = (V - Q_new / coil.capacitance - coil.resistance * I_new) / coil.inductance
        I_new = I_new + dI_dt * dt
        Q_new = Q_new + I_new * dt
    end

    # Current should now be large enough for meaningful force
    F = em_force(47.0u"m", I_new, coil)  # 3m before coil center
    @test abs(ustrip(u"N", F)) > 1.0  # At least 1 N of force
    @test F > 0.0u"N"  # Force should be in the forward direction
end

@testset "RLC Dynamics — Peak Current Estimate" begin
    # For underdamped RLC: I_peak ≈ V0 * sqrt(C/L)
    coil = CoilConfig(0.0u"m", 0.001u"H", 0.01u"Ω", 0.2u"F", 2000.0u"V", 0.0u"s", 0.003u"H/m")

    V0 = ustrip(u"V", coil.voltage)
    L = ustrip(u"H", coil.inductance)
    C = ustrip(u"F", coil.capacitance)
    R = ustrip(u"Ω", coil.resistance)

    I_peak_theory = V0 * sqrt(C / L)  # Underdamped approximation
    @test I_peak_theory > 10000  # Should be thousands of amps

    # Simulate RLC discharge and find peak |I|
    Q = C * V0
    I = 0.0
    dt = 1e-6
    I_max = 0.0

    for _ in 1:200000  # 0.2s of simulation at 1μs steps
        dI = (-Q / C - R * I) / L  # V=0 (triggered)
        I += dI * dt
        Q += I * dt
        I_max = max(I_max, abs(I))
    end

    # Peak current should be within 20% of theoretical (damping reduces it)
    @test I_max > 0.8 * I_peak_theory
    @test I_max < 1.2 * I_peak_theory
end

@testset "RLC Dynamics — Energy Conservation" begin
    # Total energy E = 0.5*L*I² + 0.5*Q²/C should decrease monotonically (resistive losses)
    coil = CoilConfig(0.0u"m", 0.001u"H", 0.01u"Ω", 0.2u"F", 2000.0u"V", 0.0u"s", 0.003u"H/m")

    V0 = ustrip(u"V", coil.voltage)
    L = ustrip(u"H", coil.inductance)
    C = ustrip(u"F", coil.capacitance)
    R = ustrip(u"Ω", coil.resistance)

    Q = C * V0
    I = 0.0
    dt = 1e-6

    E_initial = 0.5 * Q^2 / C  # All energy in capacitor initially
    E_prev = E_initial

    energy_increased = false
    for step in 1:100000
        dI = (-Q / C - R * I) / L
        I += dI * dt
        Q += I * dt

        E_current = 0.5 * L * I^2 + 0.5 * Q^2 / C

        if E_current > E_prev + 1e-6 * E_initial  # Allow tiny numerical noise
            energy_increased = true
            break
        end
        E_prev = E_current
    end

    @test !energy_increased  # Energy must never increase (2nd law)
    @test E_prev < E_initial  # Energy should be lost to resistance
    @test E_prev > 0.1 * E_initial  # Significant energy remains after 0.1s
end

# =============================================================================
# FULL ODE INTEGRATION TESTS
# Verifies that trajectory_ode_dimensionless! actually accelerates the projectile.
# =============================================================================

@testset "ODE Produces Acceleration — Earth Launcher" begin
    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 5,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 5000.0u"V",
        gradient_per_coil = 0.002u"H/m"
    )
    payload = PayloadConfig(10.0u"kg", π*(0.1u"m")^2, 0.05u"m", 0.9, 900.0u"J/(kg*K)", 288.15u"K", 2000.0u"K")
    mission = MissionProfile(0.0u"m", 32.5u"°", 34.9u"°", 90.0u"°", 85.0u"°", 1000.0u"m/s", 50.0u"km")

    u0_unitful = initial_state(launcher, payload, mission)
    n_coils = launcher.num_coils
    u0, units_ref = strip_units_state(u0_unitful, n_coils)
    p = (launcher, payload, mission, units_ref)

    # Evaluate derivative at t=0
    du = zeros(length(u0))
    trajectory_ode_dimensionless!(du, u0, p, 0.0)

    # Velocity derivatives (du[4:6]) should be non-zero — gravity at minimum
    accel_mag = sqrt(du[4]^2 + du[5]^2 + du[6]^2)
    @test accel_mag > 1.0  # At least 1 m/s² (gravity alone is ~9.8)

    # Current derivatives should be large for triggered coils
    # (capacitors discharging through RLC)
    max_dI = maximum(abs.(du[8:(7+n_coils)]))
    @test max_dI > 1e4  # At least 10,000 A/s

    # Run a few Euler steps and verify velocity increases
    dt = 1e-5
    u = copy(u0)
    for _ in 1:1000  # 0.01s
        trajectory_ode_dimensionless!(du, u, p, 0.0)
        u .+= du * dt
    end

    v_after = sqrt(u[4]^2 + u[5]^2 + u[6]^2)
    @test v_after > 10.0  # Projectile should have gained significant velocity
end

@testset "ODE Produces Acceleration — Moon Launcher" begin
    launcher = create_uniform_launcher(
        length = 500.0u"m", num_coils = 50,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 0.2u"F", voltage_per_coil = 2000.0u"V",
        gradient_per_coil = 0.003u"H/m", vacuum_pressure_ratio = 0.0001
    )
    payload = PayloadConfig(5.0u"kg", π*(0.05u"m")^2, 0.02u"m", 0.95, 900.0u"J/(kg*K)", 288.15u"K", 2500.0u"K")
    mission = MissionProfile(0.0u"m", 26.3u"°", -47.5u"°", 90.0u"°", 88.0u"°", 1800.0u"m/s", 100.0u"km", MOON)

    u0_unitful = initial_state(launcher, payload, mission)
    n_coils = launcher.num_coils
    u0, units_ref = strip_units_state(u0_unitful, n_coils)
    p = (launcher, payload, mission, units_ref)

    # Evaluate derivative at t=0
    du = zeros(length(u0))
    trajectory_ode_dimensionless!(du, u0, p, 0.0)

    # Gravity on Moon: ~1.625 m/s²
    accel_mag = sqrt(du[4]^2 + du[5]^2 + du[6]^2)
    @test accel_mag > 1.0  # At least gravity

    # Current derivatives must be non-zero for triggered coils
    max_dI = maximum(abs.(du[8:(7+n_coils)]))
    @test max_dI > 1e5  # V0/L = 2000/0.001 = 2e6

    # Run Euler steps and verify velocity increases
    dt = 1e-5
    u = copy(u0)
    for _ in 1:1000  # 0.01s
        trajectory_ode_dimensionless!(du, u, p, 0.0)
        u .+= du * dt
    end

    v_after = sqrt(u[4]^2 + u[5]^2 + u[6]^2)
    @test v_after > 10.0  # Must gain speed (not just sit there!)
end

@testset "ODE Exit Velocity Is Physically Reasonable" begin
    # Small test launcher: known energy, check exit velocity is in right ballpark
    launcher = create_uniform_launcher(
        length = 50.0u"m", num_coils = 5,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 0.5u"F", voltage_per_coil = 1000.0u"V",
        gradient_per_coil = 0.003u"H/m"
    )
    payload = PayloadConfig(5.0u"kg", π*(0.05u"m")^2, 0.02u"m", 0.9, 900.0u"J/(kg*K)", 288.15u"K", 2000.0u"K")
    mission = MissionProfile(0.0u"m", 0.0u"°", 0.0u"°", 90.0u"°", 85.0u"°", 500.0u"m/s", 10.0u"km")

    # Total stored energy
    E_stored = sum(0.5 * c.capacitance * c.voltage^2 for c in launcher.coils)
    E_stored_J = ustrip(u"J", E_stored)
    m = ustrip(u"kg", payload.mass)

    # Theoretical max exit velocity (100% efficiency, no gravity)
    v_max = sqrt(2 * E_stored_J / m)
    @test v_max > 0.0  # Sanity

    # Run ODE simulation
    u0_unitful = initial_state(launcher, payload, mission)
    n_coils = launcher.num_coils
    u0, units_ref = strip_units_state(u0_unitful, n_coils)
    p = (launcher, payload, mission, units_ref)

    du = zeros(length(u0))
    dt = 1e-6
    u = copy(u0)

    # Simulate until projectile exits launcher or 1s elapsed
    launch_pos = u[1:3]
    for step in 1:1000000
        trajectory_ode_dimensionless!(du, u, p, step * dt)
        u .+= du * dt

        x_launcher = sqrt(sum((u[1:3] .- launch_pos).^2))
        if x_launcher > ustrip(u"m", launcher.length)
            break
        end
    end

    v_exit = sqrt(u[4]^2 + u[5]^2 + u[6]^2)

    # Exit velocity must be positive and less than theoretical max
    @test v_exit > 1.0   # Must have some velocity
    @test v_exit < v_max  # Can't exceed energy limit
end

@testset "Altitude Stays Positive After Launch — Moon" begin
    # The original bug: Moon trajectories going to negative altitude.
    # This test ensures the projectile gains altitude after EM acceleration.
    launcher = create_uniform_launcher(
        length = 500.0u"m", num_coils = 50,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 0.2u"F", voltage_per_coil = 2000.0u"V",
        gradient_per_coil = 0.003u"H/m", vacuum_pressure_ratio = 0.0001
    )
    payload = PayloadConfig(5.0u"kg", π*(0.05u"m")^2, 0.02u"m", 0.95, 900.0u"J/(kg*K)", 288.15u"K", 2500.0u"K")
    mission = MissionProfile(0.0u"m", 26.3u"°", -47.5u"°", 90.0u"°", 88.0u"°", 1800.0u"m/s", 100.0u"km", MOON)

    u0_unitful = initial_state(launcher, payload, mission)
    n_coils = launcher.num_coils
    u0, units_ref = strip_units_state(u0_unitful, n_coils)
    p = (launcher, payload, mission, units_ref)

    du = zeros(length(u0))
    dt = 1e-5
    u = copy(u0)

    # Simulate 0.1s (enough for initial acceleration phase)
    for _ in 1:10000
        trajectory_ode_dimensionless!(du, u, p, 0.0)
        u .+= du * dt
    end

    # Compute altitude: distance from Moon center minus Moon radius
    r = sqrt(u[1]^2 + u[2]^2 + u[3]^2)
    R_moon = ustrip(u"m", MOON.radius)
    altitude = r - R_moon

    @test altitude >= -1.0  # Must not go significantly below surface
    # After 0.1s of EM acceleration, should have positive upward velocity
    v_mag = sqrt(u[4]^2 + u[5]^2 + u[6]^2)
    @test v_mag > 50.0  # Must have gained substantial velocity
end

println("  ✓ EM acceleration tests passed")
