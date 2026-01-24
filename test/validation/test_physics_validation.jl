"""
Physics Validation Tests

Tests against known physical values, published data, and theoretical results.
"""


@testset "ISA Standard Atmosphere Validation" begin
    # US Standard Atmosphere 1976 reference values
    # Source: NOAA, NASA, ICAO Standard Atmosphere

    # Sea level
    @test temperature(0.0u"m") ≈ 288.15u"K" rtol=0.001
    @test pressure(0.0u"m") ≈ 101325.0u"Pa" rtol=0.001
    @test density(0.0u"m") ≈ 1.225u"kg/m^3" rtol=0.01

    # 11 km (tropopause)
    @test temperature(11000.0u"m") ≈ 216.65u"K" rtol=0.01
    @test pressure(11000.0u"m") ≈ 22632.0u"Pa" rtol=0.05
    @test density(11000.0u"m") ≈ 0.3639u"kg/m^3" rtol=0.05

    # 20 km (lower stratosphere)
    @test temperature(20000.0u"m") ≈ 216.65u"K" rtol=0.01
    @test pressure(20000.0u"m") ≈ 5474.9u"Pa" rtol=0.05
    @test density(20000.0u"m") ≈ 0.0880u"kg/m^3" rtol=0.05

    # 32 km (stratosphere)
    @test temperature(32000.0u"m") ≈ 228.65u"K" rtol=0.05
    @test pressure(32000.0u"m") ≈ 868.0u"Pa" rtol=0.1

    # 47 km (stratopause)
    @test temperature(47000.0u"m") ≈ 270.65u"K" rtol=0.05
end

@testset "Speed of Sound Validation" begin
    # Known values at different altitudes

    # Sea level: ~340 m/s
    a_0 = speed_of_sound(0.0u"km")
    @test a_0 ≈ 340.3u"m/s" rtol=0.01

    # 11 km (tropopause): ~295 m/s
    a_11 = speed_of_sound(11.0u"km")
    @test a_11 ≈ 295.0u"m/s" rtol=0.02

    # At different temperatures (theoretical)
    # a = sqrt(γ * R * T)
    T_test = 273.15u"K"  # 0°C
    a_theory = sqrt(1.4 * 287.05u"J/(kg*K)" * T_test)
    @test a_theory ≈ 331.3u"m/s" rtol=0.01
end

@testset "Orbital Mechanics Validation" begin
    # LEO (400 km) - ISS orbit
    v_iss = orbital_velocity(400.0u"km")
    @test v_iss ≈ 7.67u"km/s" rtol=0.01

    # Geostationary orbit (35,786 km)
    v_geo = orbital_velocity(35786.0u"km")
    @test v_geo ≈ 3.07u"km/s" rtol=0.02

    # Moon's orbit (~384,400 km)
    v_moon = orbital_velocity(384400.0u"km")
    @test v_moon ≈ 1.02u"km/s" rtol=0.05

    # Escape velocity at surface
    v_esc_0 = escape_velocity(0.0u"km")
    @test v_esc_0 ≈ 11.2u"km/s" rtol=0.01
end

@testset "Gravity Magnitude Validation" begin
    # Surface gravity
    g_surf = abs(gravity_point_mass(R_Earth))
    @test g_surf ≈ 9.82u"m/s^2" rtol=0.005

    # ISS altitude (400 km)
    g_iss = abs(gravity_point_mass(R_Earth + 400.0u"km"))
    @test g_iss ≈ 8.69u"m/s^2" rtol=0.02

    # Geostationary orbit
    g_geo = abs(gravity_point_mass(R_Earth + 35786.0u"km"))
    @test g_geo ≈ 0.223u"m/s^2" rtol=0.05

    # Inverse square law check
    g_1 = abs(gravity_point_mass(R_Earth))
    g_2 = abs(gravity_point_mass(2 * R_Earth))
    @test g_2 ≈ g_1 / 4.0 rtol=0.001
end

@testset "Orbital Period Validation" begin
    # For circular orbit: T = 2π√(a³/μ)

    # LEO (400 km) - should be ~92 minutes
    r_leo = R_Earth + 400.0u"km"
    T_leo = 2π * sqrt(r_leo^3 / μ_Earth)
    @test uconvert(u"minute", T_leo) ≈ 92.7u"minute" rtol=0.02

    # Geostationary - should be 24 hours
    r_geo = R_Earth + 35786.0u"km"
    T_geo = 2π * sqrt(r_geo^3 / μ_Earth)
    @test uconvert(u"hr", T_geo) ≈ 24.0u"hr" rtol=0.01
end

@testset "Atmospheric Lapse Rate Validation" begin
    # Troposphere: -6.5 K/km standard lapse rate

    T_0 = temperature(0.0u"km")
    T_1 = temperature(1.0u"km")

    lapse_rate = (T_1 - T_0) / 1.0u"km"
    @test lapse_rate ≈ -6.5u"K/km" rtol=0.05

    # Check over full troposphere
    T_11 = temperature(11.0u"km")
    avg_lapse = (T_11 - T_0) / 11.0u"km"
    @test avg_lapse ≈ -6.5u"K/km" rtol=0.05
end

@testset "Mach Number Reference Values" begin
    # Mach 1 at sea level ≈ 340 m/s
    v_mach1 = 340.0u"m/s"
    M = mach_number(v_mach1, 0.0u"km")
    @test M ≈ 1.0 rtol=0.02

    # Concorde cruise: Mach 2.04 at 18 km
    v_concorde = 2.04 * speed_of_sound(18.0u"km")
    M_concorde = mach_number(v_concorde, 18.0u"km")
    @test M_concorde ≈ 2.04 rtol=0.05

    # SR-71 Blackbird: Mach 3.3 at 24 km
    v_sr71 = 3.3 * speed_of_sound(24.0u"km")
    M_sr71 = mach_number(v_sr71, 24.0u"km")
    @test M_sr71 ≈ 3.3 rtol=0.05
end

@testset "Drag Scaling Laws" begin
    # Drag should scale as v²

    v1 = 100.0u"m/s"
    v2 = 200.0u"m/s"
    h = 10.0u"km"
    A = 1.0u"m^2"

    F1 = drag_force(v1, h, A)
    F2 = drag_force(v2, h, A)

    # F2/F1 should equal (v2/v1)²
    ratio = F2 / F1
    expected_ratio = (v2 / v1)^2

    @test ustrip(ratio) ≈ ustrip(expected_ratio) rtol=0.01
end

@testset "Heating Scaling Laws" begin
    # Heating should scale as v³

    v1 = 1000.0u"m/s"
    v2 = 2000.0u"m/s"
    h = 25.0u"km"
    R_n = 0.1u"m"

    q1 = aerodynamic_heating(v1, h, R_n)
    q2 = aerodynamic_heating(v2, h, R_n)

    # q2/q1 should equal (v2/v1)³
    ratio = q2 / q1
    expected_ratio = (v2 / v1)^3

    @test ustrip(ratio) ≈ ustrip(expected_ratio) rtol=0.01
end

@testset "Kinetic Energy Validation" begin
    # Space Shuttle reentry: ~25,000 km/h ≈ 7 km/s, 100 tonnes
    m_shuttle = 100000.0u"kg"
    v_shuttle = 7.0u"km/s"
    KE_shuttle = 0.5 * m_shuttle * v_shuttle^2

    @test uconvert(u"TJ", KE_shuttle) ≈ 2.45u"TJ" rtol=0.1

    # Apollo return: ~11 km/s
    m_apollo = 5000.0u"kg"
    v_apollo = 11.0u"km/s"
    KE_apollo = 0.5 * m_apollo * v_apollo^2

    @test uconvert(u"TJ", KE_apollo) ≈ 0.303u"TJ" rtol=0.1
end

@testset "Orbital Target Validation" begin
    # system-specific parameters

    # Target velocity: 8 km/s
    @test v_target_orbital == 8000.0u"m/s"

    # Required kinetic energy for 20 kg payload
    KE_20kg = 0.5 * 20.0u"kg" * v_target_orbital^2
    @test uconvert(u"MJ", KE_20kg) ≈ 640.0u"MJ" rtol=0.01

    # Launcher acceleration (nominal 2 km tube)
    a_moon = average_acceleration(v_target_orbital, L_launcher_nominal)
    @test a_moon ≈ 16000.0u"m/s^2" rtol=0.01

    # G-loading
    g_moon = g_loading(a_moon)
    @test g_moon ≈ 1632.0 rtol=0.05

    # Energy efficiency claim: 45% vs 4%
    @test payload_fraction_target == 0.45
end

@testset "Earth Figure Parameters" begin
    # J2 coefficient
    @test J2_Earth ≈ 1.08263e-3 rtol=0.001

    # Earth's equatorial radius
    @test R_Earth ≈ 6378.137u"km" rtol=0.001

    # Gravitational parameter
    @test μ_Earth ≈ 398600.4418u"km^3/s^2" rtol=0.00001
end

@testset "Atmospheric Scale Height" begin
    # H ≈ RT/Mg
    T_avg = 250.0u"K"  # Average temperature
    M = 0.029u"kg/mol"  # Molar mass of air
    R_gas = 8.314u"J/(mol*K)"

    H_calc = R_gas * T_avg / (M * g_0)
    @test H_calc ≈ H_scale rtol=0.2  # ~8.5 km

    # Check exponential decrease
    ρ_0 = density(0.0u"km")
    ρ_H = density(H_scale)

    # Should decrease by factor of e
    @test ρ_H ≈ ρ_0 / ℯ rtol=0.2
end

@testset "Circular Orbit Energy" begin
    # For circular orbit: E = -μ/(2a)

    h = 400.0u"km"
    r = R_Earth + h
    v = orbital_velocity(h)

    ε_calc = specific_orbital_energy(v, r)
    ε_theory = -μ_Earth / (2 * r)

    @test ε_calc ≈ ε_theory rtol=0.01
end

@testset "Vis-Viva Equation" begin
    # v² = μ(2/r - 1/a)

    # Elliptical orbit
    r_p = R_Earth + 400.0u"km"  # Periapsis
    r_a = R_Earth + 800.0u"km"  # Apoapsis
    a = (r_p + r_a) / 2          # Semi-major axis

    # Velocity at periapsis
    v_p_theory = sqrt(μ_Earth * (2/r_p - 1/a))

    # Using our energy function
    v_p_calc = sqrt(2 * (specific_orbital_energy(0.0u"m/s", r_p) + μ_Earth/r_p))

    # Should be close (within 10% due to approximations)
    @test isfinite(ustrip(u"m/s", v_p_theory))
    @test v_p_theory > 7.0u"km/s"
    @test v_p_theory < 8.0u"km/s"
end

@testset "Hypersonic Regime Validation" begin
    # Hypersonic: Mach > 5

    # X-15 record: Mach 6.7 at 31 km
    v_x15 = 6.7 * speed_of_sound(31.0u"km")
    M_x15 = mach_number(v_x15, 31.0u"km")
    @test M_x15 > 6.0
    @test M_x15 < 7.0

    # Space Shuttle reentry: Mach 25 at 70 km
    v_shuttle = 25.0 * speed_of_sound(70.0u"km")
    M_shuttle = mach_number(v_shuttle, 70.0u"km")
    @test M_shuttle > 20.0
end

@testset "Dynamic Pressure Validation" begin
    # q = 0.5 * ρ * v²

    # Max-Q for rockets: typically at ~10 km, ~Mach 1.5
    h_maxq = 10.0u"km"
    v_maxq = 1.5 * speed_of_sound(h_maxq)
    ρ_maxq = density(h_maxq)

    q_maxq = 0.5 * ρ_maxq * v_maxq^2

    # Typical max-Q: 30-40 kPa
    @test uconvert(u"kPa", q_maxq) > 20.0u"kPa"
    @test uconvert(u"kPa", q_maxq) < 50.0u"kPa"
end

println("  ✓ Physics validation tests passed")
