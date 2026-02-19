"""
Tests for gravity.jl with Moon body parameter

Tests:
- Lunar gravity values
- Lunar orbital/escape velocities
- Body parameter forwarding
"""

@testset "Moon Point Mass Gravity" begin
    # At Moon surface
    r = MOON.radius
    a_g = gravity_point_mass(r; body=MOON)

    @test dimension(a_g) == dimension(1.0u"m/s^2")
    @test abs(a_g) ≈ 1.625u"m/s^2" rtol=0.01  # Lunar surface gravity

    # Inverse square law still holds
    a_g_2r = gravity_point_mass(2*r; body=MOON)
    @test abs(a_g_2r) ≈ abs(a_g) / 4.0 rtol=0.01
end

@testset "Moon J2 Gravity" begin
    r = MOON.radius
    pos = SVector(r, 0.0u"m", 0.0u"m")
    a_j2 = gravity_with_J2(pos; body=MOON)

    @test dimension(a_j2[1]) == dimension(1.0u"m/s^2")

    # J2 is a small perturbation
    a_point = gravity_point_mass(pos; body=MOON)
    @test norm(a_j2 - a_point) / norm(a_point) < 0.01
end

@testset "Moon Altitude Calculations" begin
    r = MOON.radius
    pos = SVector(r, 0.0u"m", 0.0u"m")
    h = altitude_from_position(pos; body=MOON)
    @test h ≈ 0.0u"m" atol=1.0u"m"

    # At 100 km altitude
    pos_100 = SVector(r + 100.0u"km", 0.0u"m", 0.0u"m")
    h_100 = altitude_from_position(pos_100; body=MOON)
    @test h_100 ≈ 100.0u"km" rtol=1e-6
end

@testset "Moon Position from Altitude" begin
    h = 100.0u"km"
    lat = 0.0u"°"
    lon = 0.0u"°"
    pos = position_from_altitude_angle(h, lat, lon; body=MOON)

    @test length(pos) == 3
    @test dimension(pos[1]) == dimension(1.0u"m")

    # Check altitude roundtrip
    h_check = altitude_from_position(pos; body=MOON)
    @test h_check ≈ h rtol=1e-4
end

@testset "Low Lunar Orbit Velocity" begin
    # LLO at 100 km altitude
    h = 100.0u"km"
    v_llo = orbital_velocity(h; body=MOON)

    @test dimension(v_llo) == dimension(1.0u"m/s")
    # v_LLO ≈ 1634 m/s (well-known value)
    @test v_llo ≈ 1634.0u"m/s" rtol=0.02
end

@testset "Moon Escape Velocity" begin
    # At surface
    v_esc = escape_velocity(0.0u"km"; body=MOON)

    @test dimension(v_esc) == dimension(1.0u"m/s")
    # v_escape ≈ 2376 m/s (well-known value)
    @test v_esc ≈ 2376.0u"m/s" rtol=0.02

    # v_escape = sqrt(2) * v_orbital at same altitude
    v_orb = orbital_velocity(0.0u"km"; body=MOON)
    @test v_esc ≈ sqrt(2) * v_orb rtol=0.01

    # Lower at higher altitude
    v_esc_100 = escape_velocity(100.0u"km"; body=MOON)
    @test v_esc_100 < v_esc
end

@testset "Moon Specific Orbital Energy" begin
    h = 100.0u"km"
    v = orbital_velocity(h; body=MOON)
    r = MOON.radius + h

    ε = specific_orbital_energy(v, r; body=MOON)

    @test dimension(ε) == dimension(1.0u"m^2/s^2")
    @test ε < 0.0u"m^2/s^2"  # Bound orbit

    # Escape velocity → zero energy
    v_esc = escape_velocity(h; body=MOON)
    ε_esc = specific_orbital_energy(v_esc, r; body=MOON)
    @test ε_esc ≈ 0.0u"m^2/s^2" atol=1.0u"m^2/s^2"
end

@testset "Moon Orbital Elements" begin
    h = 100.0u"km"
    r = MOON.radius + h
    v_orb = orbital_velocity(h; body=MOON)

    position = SVector(r, 0.0u"m", 0.0u"m")
    velocity = SVector(0.0u"m/s", v_orb, 0.0u"m/s")

    elements = orbital_elements(position, velocity; body=MOON)

    @test elements.a ≈ r rtol=1e-3
    @test elements.e < 0.01  # Circular orbit
end

@testset "Earth vs Moon Comparison" begin
    # Same altitude, different body → different results
    h = 100.0u"km"

    v_earth = orbital_velocity(h; body=EARTH)
    v_moon = orbital_velocity(h; body=MOON)

    @test v_earth > v_moon  # Earth requires faster orbital velocity

    # Escape velocity comparison
    v_esc_earth = escape_velocity(0.0u"km"; body=EARTH)
    v_esc_moon = escape_velocity(0.0u"km"; body=MOON)

    @test v_esc_earth > v_esc_moon
    @test v_esc_earth ≈ 11.2u"km/s" rtol=0.01
    @test v_esc_moon ≈ 2.376u"km/s" rtol=0.02
end

println("  ✓ Lunar gravity tests passed")
