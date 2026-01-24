"""
Edge Cases and Boundary Condition Tests

Tests:
- Extreme parameter values
- Zero and infinity handling
- Boundary conditions
- Invalid inputs
- Numerical edge cases
"""


@testset "Atmosphere Edge Cases" begin
    # Sea level
    @test density(0.0u"m") > 0.0u"kg/m^3"
    @test temperature(0.0u"m") > 0.0u"K"
    @test pressure(0.0u"m") > 0.0u"Pa"

    # Very high altitude (space)
    h_space = 500.0u"km"
    @test density(h_space) > 0.0u"kg/m^3"  # Should be very small but positive
    @test density(h_space) < 1e-10u"kg/m^3"

    # Karman line
    h_karman = 100.0u"km"
    ρ_karman = density(h_karman)
    @test ρ_karman > 0.0u"kg/m^3"
    @test ρ_karman < density(50.0u"km")

    # Extreme altitude
    h_extreme = 1000.0u"km"
    @test density(h_extreme) > 0.0u"kg/m^3"
    @test isfinite(ustrip(u"kg/m^3", density(h_extreme)))
end

@testset "Drag Force Edge Cases" begin
    # Zero velocity
    F_zero = drag_force(0.0u"m/s", 10.0u"km", 1.0u"m^2")
    @test F_zero == 0.0u"N"

    # Very high velocity
    F_high = drag_force(10000.0u"m/s", 10.0u"km", 1.0u"m^2")
    @test F_high > 0.0u"N"
    @test isfinite(ustrip(u"N", F_high))

    # Zero area
    F_zero_area = drag_force(1000.0u"m/s", 10.0u"km", 0.0u"m^2")
    @test F_zero_area == 0.0u"N"

    # Very small area
    F_small = drag_force(1000.0u"m/s", 10.0u"km", 1e-6u"m^2")
    @test F_small > 0.0u"N"
    @test F_small < 1.0u"N"

    # High altitude (thin atmosphere)
    F_thin = drag_force(5000.0u"m/s", 200.0u"km", 1.0u"m^2")
    @test F_thin >= 0.0u"N"
    @test F_thin < drag_force(5000.0u"m/s", 10.0u"km", 1.0u"m^2")
end

@testset "Heating Edge Cases" begin
    # Zero velocity
    q_zero = aerodynamic_heating(0.0u"m/s", 25.0u"km", 0.1u"m")
    @test q_zero == 0.0u"W/m^2"

    # Very high velocity (reentry speeds)
    q_reentry = aerodynamic_heating(11000.0u"m/s", 50.0u"km", 0.1u"m")
    @test q_reentry > 1.0u"MW/m^2"
    @test isfinite(ustrip(u"W/m^2", q_reentry))

    # Very small nose radius (sharp nose)
    q_sharp = aerodynamic_heating(5000.0u"m/s", 25.0u"km", 0.01u"m")
    @test q_sharp > aerodynamic_heating(5000.0u"m/s", 25.0u"km", 0.1u"m")

    # Space (no atmosphere)
    q_space = aerodynamic_heating(8000.0u"m/s", 500.0u"km", 0.1u"m")
    @test q_space >= 0.0u"W/m^2"
    @test q_space < 1000.0u"W/m^2"  # Should be very small
end

@testset "Gravity Edge Cases" begin
    # At Earth surface
    a_surface = gravity_point_mass(R_Earth)
    @test abs(a_surface) ≈ g_0 rtol=0.01

    # Very close to Earth (but not inside)
    a_close = gravity_point_mass(R_Earth + 1.0u"m")
    @test abs(a_close) > 0.0u"m/s^2"
    @test abs(a_close) ≈ abs(a_surface) rtol=0.001

    # Very far from Earth (Moon distance)
    r_moon = 384400.0u"km"
    a_moon = gravity_point_mass(r_moon)
    @test abs(a_moon) > 0.0u"m/s^2"
    @test abs(a_moon) < 0.01u"m/s^2"

    # Geostationary orbit
    r_geo = R_Earth + 35786.0u"km"
    a_geo = gravity_point_mass(r_geo)
    @test abs(a_geo) > 0.0u"m/s^2"
    @test abs(a_geo) < 1.0u"m/s^2"
end

@testset "Orbital Velocity Edge Cases" begin
    # Very low orbit (just above surface)
    v_low = orbital_velocity(1.0u"km")
    @test v_low > 7.0u"km/s"
    @test v_low < 8.0u"km/s"

    # Very high orbit (lunar orbit)
    v_high = orbital_velocity(384400.0u"km")
    @test v_high > 0.0u"km/s"
    @test v_high < 2.0u"km/s"

    # At surface (theoretical)
    v_surface = orbital_velocity(0.0u"m")
    @test v_surface > 7.0u"km/s"
    @test v_surface < 8.0u"km/s"
end

@testset "Mach Number Edge Cases" begin
    # Subsonic
    M_sub = mach_number(100.0u"m/s", 0.0u"km")
    @test M_sub < 1.0
    @test M_sub > 0.0

    # Sonic
    M_sonic = mach_number(340.0u"m/s", 0.0u"km")
    @test M_sonic ≈ 1.0 rtol=0.05

    # Supersonic
    M_super = mach_number(1000.0u"m/s", 0.0u"km")
    @test M_super > 1.0
    @test M_super < 5.0

    # Hypersonic
    M_hyper = mach_number(8000.0u"m/s", 50.0u"km")
    @test M_hyper > 10.0

    # Zero velocity
    M_zero = mach_number(0.0u"m/s", 0.0u"km")
    @test M_zero == 0.0
end

@testset "Temperature Bounds" begin
    # Temperature should never be negative
    for h in [0.0, 11.0, 25.0, 50.0, 100.0, 200.0] .* 1.0u"km"
        T = temperature(h)
        @test T > 0.0u"K"
        @test T < 400.0u"K"  # Reasonable upper bound for atmosphere
    end

    # Temperature should be finite everywhere
    for h in range(0.0u"km", 500.0u"km", length=50)
        T = temperature(h)
        @test isfinite(ustrip(u"K", T))
    end
end

@testset "Density Monotonicity" begin
    # Density should decrease monotonically
    altitudes = range(0.0u"km", 100.0u"km", length=100)
    densities = density.(altitudes)

    for i in 2:length(densities)
        @test densities[i] < densities[i-1]
    end

    # All densities should be positive
    for ρ in densities
        @test ρ > 0.0u"kg/m^3"
        @test isfinite(ustrip(u"kg/m^3", ρ))
    end
end

@testset "Pressure Monotonicity" begin
    # Pressure should decrease monotonically
    altitudes = range(0.0u"km", 100.0u"km", length=100)
    pressures = pressure.(altitudes)

    for i in 2:length(pressures)
        @test pressures[i] < pressures[i-1]
    end

    # All pressures should be positive
    for P in pressures
        @test P > 0.0u"Pa"
        @test isfinite(ustrip(u"Pa", P))
    end
end

@testset "Escape Velocity Bounds" begin
    # Escape velocity should be sqrt(2) times orbital velocity
    for h in [0.0, 400.0, 1000.0, 10000.0] .* 1.0u"km"
        v_esc = escape_velocity(h)
        v_orb = orbital_velocity(h)

        @test v_esc ≈ sqrt(2) * v_orb rtol=0.001
        @test v_esc > v_orb
    end
end

@testset "Specific Energy Signs" begin
    # Circular orbit - negative energy
    h = 400.0u"km"
    v_circ = orbital_velocity(h)
    r = R_Earth + h
    ε_circ = specific_orbital_energy(v_circ, r)
    @test ε_circ < 0.0u"m^2/s^2"

    # Escape trajectory - zero energy
    v_esc = escape_velocity(h)
    ε_esc = specific_orbital_energy(v_esc, r)
    @test abs(ε_esc) < 100.0u"m^2/s^2"  # Near zero

    # Hyperbolic - positive energy
    v_hyper = 1.5 * v_esc
    ε_hyper = specific_orbital_energy(v_hyper, r)
    @test ε_hyper > 0.0u"m^2/s^2"
end

@testset "Reynolds Number Ranges" begin
    # Low Reynolds (laminar)
    Re_low = reynolds_number(1.0u"m/s", 0.0u"km", 0.01u"m")
    @test Re_low < 1e5

    # High Reynolds (turbulent)
    Re_high = reynolds_number(1000.0u"m/s", 0.0u"km", 1.0u"m")
    @test Re_high > 1e6

    # Space (very low density)
    Re_space = reynolds_number(8000.0u"m/s", 500.0u"km", 1.0u"m")
    @test Re_space >= 0.0
    @test Re_space < 1e3
end

@testset "Altitude Conversion Consistency" begin
    # Round trip: position -> altitude -> position
    test_altitudes = [0.0, 100.0, 400.0, 1000.0, 10000.0] .* 1.0u"km"

    for h in test_altitudes
        # Create position
        pos = position_from_altitude_angle(h, 0.0u"°", 0.0u"°")

        # Convert back to altitude
        h_back = altitude_from_position(pos)

        # Should match
        @test h_back ≈ h rtol=1e-6
    end
end

@testset "G-Loading Ranges" begin
    # Normal gravity
    g_1 = g_loading(g_0)
    @test g_1 ≈ 1.0 rtol=1e-10

    # High-g (fighter jet)
    g_9 = g_loading(9.0 * g_0)
    @test g_9 ≈ 9.0 rtol=1e-10

    # EM launcher average
    g_moon = g_loading(a_avg_nominal)
    @test g_moon > 1000.0
    @test g_moon < 2000.0

    # Extreme (railgun)
    g_extreme = g_loading(100000.0u"m/s^2")
    @test g_extreme > 10000.0
end

@testset "Coil Spacing Edge Cases" begin
    # Many coils
    spacing_many = coil_spacing(1000, 1000.0u"m")
    @test spacing_many == 1.0u"m"

    # Few coils
    spacing_few = coil_spacing(10, 1000.0u"m")
    @test spacing_few == 100.0u"m"

    # One coil
    spacing_one = coil_spacing(1, 1000.0u"m")
    @test spacing_one == 1000.0u"m"
end

@testset "Energy Calculation Bounds" begin
    # Small mass, low velocity
    E_small = 0.5 * 1.0u"kg" * (100.0u"m/s")^2
    @test E_small > 0.0u"J"
    @test E_small < 10.0u"kJ"

    # nominal
    E_nominal = 0.5 * m_payload_nominal * v_target_orbital^2
    @test E_nominal > 600.0u"MJ"
    @test E_nominal < 700.0u"MJ"

    # Large mass, high velocity
    E_large = 0.5 * 1000.0u"kg" * (10000.0u"m/s")^2
    @test uconvert(u"GJ", E_large) >= 50.0u"GJ"
end

@testset "Physical Constant Precision" begin
    # Check that constants have reasonable precision
    @test μ_Earth ≈ 3.986004418e14u"m^3/s^2" rtol=1e-10
    @test R_Earth ≈ 6.3781370e6u"m" rtol=1e-8
    @test g_0 ≈ 9.80665u"m/s^2" rtol=1e-10

    # Standard atmosphere
    @test ρ_0 ≈ 1.225u"kg/m^3" rtol=0.01
    @test P_0 ≈ 101325.0u"Pa" rtol=0.001
    @test T_0 ≈ 288.15u"K" rtol=0.01
end

println("  ✓ Edge cases and boundary condition tests passed")
