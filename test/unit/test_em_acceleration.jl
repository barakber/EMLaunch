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

    # Before trigger
    V_before = coil_voltage(0.3u"s", coil)
    @test V_before == 0.0u"V"

    # After trigger
    V_after = coil_voltage(0.6u"s", coil)
    @test V_after == 5000.0u"V"

    # At trigger
    V_at = coil_voltage(0.5u"s", coil)
    @test V_at == 5000.0u"V"
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

println("  ✓ EM acceleration tests passed")
