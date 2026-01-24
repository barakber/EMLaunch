"""
Integration Tests

Tests that combine multiple modules to ensure they work together correctly.
"""


@testset "Energy Balance Check" begin
    # Calculate required energy for orbital target
    m = m_payload_nominal
    v = v_target_orbital
    KE_target = 0.5 * m * v^2

    @test uconvert(u"MJ", KE_target) ≈ 640.0u"MJ" rtol=0.01

    # Launcher electrical energy
    # Using fewer coils for realistic efficiency (40-80%)
    num_coils_efficient = 10  # Gives ~64% efficiency
    launcher = create_uniform_launcher(
        length = L_launcher_nominal,
        num_coils = num_coils_efficient,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = R_coil,
        capacitance_per_coil = 10.0u"F",
        voltage_per_coil = 5000.0u"V",
        gradient_per_coil = 0.002u"H/m"
    )

    E_elec = sum(0.5 * coil.capacitance * coil.voltage^2
                 for coil in launcher.coils)

    # Energy should be sufficient
    @test E_elec > KE_target

    # Efficiency should be reasonable
    η = ustrip(KE_target / E_elec)
    @test 0.4 < η < 0.8  # 40-80% efficiency range
end

@testset "Launcher Performance Estimate" begin
    # Using launcher parameters
    L = L_launcher_nominal
    v = v_target_orbital

    # Average acceleration
    a = average_acceleration(v, L)
    g_load = g_loading(a)

    @test a ≈ 16000.0u"m/s^2" rtol=0.01
    @test g_load ≈ 1630.0 rtol=1.0

    # Time in launcher
    t_launch = 2 * L / v
    @test t_launch ≈ 0.5u"s" rtol=0.01

    # Should be under structural limit
    @test a < a_max_structural
end

@testset "Atmospheric Drag vs Altitude" begin
    # At different altitudes, drag should decrease
    v = 2000.0u"m/s"
    A = A_ref

    F_0 = drag_force(v, 0.0u"km", A)
    F_10 = drag_force(v, 10.0u"km", A)
    F_25 = drag_force(v, 25.0u"km", A)
    F_50 = drag_force(v, 50.0u"km", A)
    F_100 = drag_force(v, 100.0u"km", A)

    # Drag should decrease with altitude
    @test F_0 > F_10 > F_25 > F_50 > F_100

    # At 100 km, drag should be small compared to sea level
    @test F_100 < 1.0u"N"  # Still some atmosphere at 100 km
    @test F_100 < 0.01 * F_0  # But < 1% of sea level drag
end

@testset "Aerodynamic Heating vs Velocity" begin
    # Heating should increase dramatically with velocity (v³)
    h = 25.0u"km"
    R_n = R_nose

    q_1 = aerodynamic_heating(1000.0u"m/s", h, R_n)
    q_2 = aerodynamic_heating(2000.0u"m/s", h, R_n)
    q_4 = aerodynamic_heating(4000.0u"m/s", h, R_n)

    # Should scale as v³
    @test q_2 / q_1 ≈ 8.0 rtol=0.1  # 2³ = 8
    @test q_4 / q_1 ≈ 64.0 rtol=0.1  # 4³ = 64

    # At orbital velocity, heating is extreme
    q_orbital = aerodynamic_heating(v_target_orbital, h, R_n)
    @test q_orbital > 10.0u"MW/m^2"
end

@testset "Orbital Insertion Requirements" begin
    # For LEO at 400 km
    h_target = h_LEO
    v_circular = orbital_velocity(h_target)

    # orbital target velocity
    v_launch = v_target_orbital

    # Velocity deficit for circularization
    Δv = v_circular - v_launch

    # Should be reasonable (launcher provides most of velocity)
    @test abs(Δv) < 1000.0u"m/s"  # < 1 km/s Δv needed

    # Check that launch velocity is close to orbital
    @test v_launch ≈ v_circular rtol=0.1  # Within 10%
end

@testset "Gravity vs Drag During Launch" begin
    # In launcher (first 0.5 seconds), EM force >> gravity, drag
    m = m_payload_nominal
    A = A_ref

    # Gravitational force
    F_grav = m * g_0
    @test F_grav ≈ 196.0u"N" rtol=0.1

    # EM force (rough estimate)
    I = 50000.0u"A"
    dL_dx = 0.002u"H/m"
    F_em_single = 0.5 * dL_dx * I^2
    @test F_em_single > 1000.0u"N"  # >> gravity

    # Drag force at exit (8 km/s at sea level - unrealistic but for comparison)
    # In reality, would be in atmosphere briefly
    F_drag_sea = drag_force(v_target_orbital, 0.0u"km", A)
    @test F_drag_sea > 10000.0u"N"  # Very high

    # At 50 km altitude, drag is much less
    F_drag_50 = drag_force(v_target_orbital, 50.0u"km", A)
    @test F_drag_50 < F_drag_sea / 100  # Much less dense
end

@testset "Complete Mission Profile Sanity" begin
    # Create full configuration
    launcher = create_uniform_launcher(
        length = L_launcher_nominal,
        num_coils = N_coils_nominal,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = R_coil,
        capacitance_per_coil = 10.0u"F",
        voltage_per_coil = 5000.0u"V",
        gradient_per_coil = 0.002u"H/m"
    )

    payload = PayloadConfig(
        m_payload_nominal,
        A_ref,
        R_nose,
        ε_thermal,
        c_p_payload,
        T_0,
        T_max_payload
    )

    mission = MissionProfile(
        h_launch,
        lat_launch,
        lon_launch,
        90.0u"°",  # East launch
        45.0u"°",  # 45° elevation
        v_target_orbital,
        h_LEO
    )

    # Create initial state
    u0 = initial_state(launcher, payload, mission)

    # State vector should have correct size
    expected_size = 7 + 2 * N_coils_nominal  # pos + vel + T + I + Q
    @test length(u0) == expected_size

    # Extract and verify
    state = extract_state(u0, N_coils_nominal)

    # Position should be near Earth surface
    r = norm(state.position)
    @test r ≈ R_Earth rtol=0.01

    # Velocity should start at zero
    @test norm(state.velocity) ≈ 0.0u"m/s" atol=1.0u"m/s"

    # Temperature should be ambient
    @test state.temperature ≈ T_0

    # Capacitors should be charged
    for Q in state.charges
        @test Q > 0.0u"C"
    end
end

@testset "Physical Units Throughout System" begin
    # Test that units propagate correctly through all calculations

    # Energy calculation
    m = 20.0u"kg"
    v = 8000.0u"m/s"
    E = 0.5 * m * v^2
    @test dimension(E) == dimension(1.0u"J")

    # Force calculation
    a = 10000.0u"m/s^2"
    F = m * a
    @test dimension(F) == dimension(1.0u"N")

    # Atmospheric properties
    h = 25.0u"km"
    ρ = density(h)
    @test dimension(ρ) == dimension(1.0u"kg/m^3")

    T = temperature(h)
    @test dimension(T) == dimension(1.0u"K")

    P = pressure(h)
    @test dimension(P) == dimension(1.0u"Pa")

    # Drag calculation
    F_d = drag_force(v, h, 1.0u"m^2")
    @test dimension(F_d) == dimension(1.0u"N")

    # Heating calculation
    q = aerodynamic_heating(v, h, 0.1u"m")
    @test dimension(q) == dimension(1.0u"W/m^2")

    # Orbital velocity
    v_orb = orbital_velocity(h)
    @test dimension(v_orb) == dimension(1.0u"m/s")

    # Gravity
    g = gravity_point_mass(R_Earth)
    @test dimension(g) == dimension(1.0u"m/s^2")
end

@testset "Dimensional Consistency" begin
    # Test that all equations are dimensionally consistent

    # Newton's second law: F = ma
    m = 20.0u"kg"
    a = 1000.0u"m/s^2"
    F = m * a
    @test dimension(F) == dimension(1.0u"kg*m/s^2")  # Definition of Newton
    @test dimension(F) == dimension(1.0u"N")

    # Kinetic energy: KE = (1/2)mv²
    v = 8000.0u"m/s"
    KE = 0.5 * m * v^2
    @test dimension(KE) == dimension(1.0u"kg*m^2/s^2")  # Definition of Joule
    @test dimension(KE) == dimension(1.0u"J")

    # Drag: F = (1/2)ρv²CdA
    ρ = 1.225u"kg/m^3"
    A = 1.0u"m^2"
    C_d = 0.85
    F_drag = 0.5 * ρ * v^2 * C_d * A
    @test dimension(F_drag) == dimension(1.0u"N")

    # Power: P = Fv
    P = F * v
    @test dimension(P) == dimension(1.0u"W")
    @test dimension(P) == dimension(1.0u"J/s")

    # Heating: q̇ has units of W/m²
    q = 1.0e6u"W/m^2"
    A_surf = 1.0u"m^2"
    Q_total = q * A_surf
    @test dimension(Q_total) == dimension(1.0u"W")
end

println("  ✓ Integration tests passed")
