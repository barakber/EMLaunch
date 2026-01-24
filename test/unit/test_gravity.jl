"""
Tests for gravity.jl

Tests:
- Point mass gravity calculations
- J2 perturbation
- Orbital velocities
- Escape velocities
- Orbital elements
- Position/altitude conversions
"""


@testset "Point Mass Gravity" begin
    # At Earth surface
    r = 6.3781370e6u"m"
    a_g = gravity_point_mass(r)

    @test dimension(a_g) == dimension(1.0u"m/s^2")
    @test abs(a_g) ≈ 9.82u"m/s^2" rtol=0.01  # Should be ≈ g_0

    # Inverse square law
    a_g_2r = gravity_point_mass(2*r)
    @test abs(a_g_2r) ≈ abs(a_g) / 4.0 rtol=0.01

    # Vector position
    pos = SVector(r, 0.0u"m", 0.0u"m")
    a_vec = gravity_point_mass(pos)

    @test dimension(a_vec[1]) == dimension(1.0u"m/s^2")
    @test norm(a_vec) ≈ abs(a_g) rtol=0.01
    @test a_vec[1] < 0.0u"m/s^2"  # Points toward center
end

@testset "J2 Perturbation" begin
    # At equator (z=0), J2 should reduce radial acceleration
    r = 6.3781370e6u"m"
    pos_equator = SVector(r, 0.0u"m", 0.0u"m")
    a_j2 = gravity_with_J2(pos_equator)

    @test dimension(a_j2[1]) == dimension(1.0u"m/s^2")

    # At pole (x=y=0), J2 effect is different
    pos_pole = SVector(0.0u"m", 0.0u"m", r)
    a_j2_pole = gravity_with_J2(pos_pole)

    # J2 should cause difference between equator and pole
    @test norm(a_j2) != norm(a_j2_pole)

    # J2 is small perturbation (should be close to point mass)
    a_point = gravity_point_mass(pos_equator)
    @test norm(a_j2 - a_point) / norm(a_point) < 0.01  # < 1% difference
end

@testset "Altitude from Position" begin
    # At surface
    r = 6.3781370e6u"m"
    pos = SVector(r, 0.0u"m", 0.0u"m")
    h = altitude_from_position(pos)

    @test h ≈ 0.0u"m" atol=1.0u"m"

    # At 400 km altitude
    pos_leo = SVector(r + 400.0u"km", 0.0u"m", 0.0u"m")
    h_leo = altitude_from_position(pos_leo)
    @test h_leo ≈ 400.0u"km" rtol=1e-6

    # Scalar radius
    h_scalar = altitude_from_position(r + 100.0u"km")
    @test h_scalar ≈ 100.0u"km" rtol=1e-6
end

@testset "Position from Altitude and Angles" begin
    # Equator, prime meridian
    h = 400.0u"km"
    lat = 0.0u"°"
    lon = 0.0u"°"

    pos = position_from_altitude_angle(h, lat, lon)

    @test length(pos) == 3
    @test dimension(pos[1]) == dimension(1.0u"m")

    # Should be on x-axis
    @test abs(pos[1]) > abs(pos[2])
    @test abs(pos[1]) > abs(pos[3])

    # Check altitude
    h_check = altitude_from_position(pos)
    @test h_check ≈ h rtol=1e-4

    # North pole
    pos_pole = position_from_altitude_angle(h, 90.0u"°", 0.0u"°")
    @test abs(pos_pole[3]) > abs(pos_pole[1])  # Should be on z-axis
end

@testset "Circular Orbital Velocity" begin
    # LEO at 400 km
    h = 400.0u"km"
    v_orb = orbital_velocity(h)

    @test dimension(v_orb) == dimension(1.0u"m/s")
    @test v_orb ≈ 7.67u"km/s" rtol=0.01

    # GEO at ~35,786 km
    h_geo = 35786.0u"km"
    v_geo = orbital_velocity(h_geo)
    @test v_geo ≈ 3.07u"km/s" rtol=0.01

    # Higher orbit → slower velocity
    @test v_orb > v_geo
end

@testset "Escape Velocity" begin
    # At surface
    h = 0.0u"km"
    v_esc = escape_velocity(h)

    @test dimension(v_esc) == dimension(1.0u"m/s")
    @test v_esc ≈ 11.2u"km/s" rtol=0.01

    # Escape velocity = √2 × orbital velocity
    v_orb = orbital_velocity(h)
    @test v_esc ≈ sqrt(2) * v_orb rtol=0.01

    # At LEO altitude
    h_leo = 400.0u"km"
    v_esc_leo = escape_velocity(h_leo)
    @test v_esc_leo < v_esc  # Lower at higher altitude
    @test v_esc_leo ≈ 10.8u"km/s" rtol=0.02
end

@testset "Specific Orbital Energy" begin
    # Circular orbit at 400 km
    h = 400.0u"km"
    v = orbital_velocity(h)
    r = 6.3781370e6u"m" + h

    ε = specific_orbital_energy(v, r)

    @test dimension(ε) == dimension(1.0u"m^2/s^2")  # Or J/kg
    @test ε < 0.0u"m^2/s^2"  # Bound orbit

    # Escape velocity → zero energy
    v_esc = escape_velocity(h)
    ε_esc = specific_orbital_energy(v_esc, r)
    @test ε_esc ≈ 0.0u"m^2/s^2" atol=1.0u"m^2/s^2"

    # Hyperbolic → positive energy
    v_hyp = 1.5 * v_esc
    ε_hyp = specific_orbital_energy(v_hyp, r)
    @test ε_hyp > 0.0u"m^2/s^2"
end

@testset "Orbital Elements - Circular Orbit" begin
    # Perfect circular orbit at 400 km
    h = 400.0u"km"
    r = 6.3781370e6u"m" + h
    v_orb = orbital_velocity(h)

    # Position and velocity for circular equatorial orbit
    position = SVector(r, 0.0u"m", 0.0u"m")
    velocity = SVector(0.0u"m/s", v_orb, 0.0u"m/s")

    elements = orbital_elements(position, velocity)

    # Semi-major axis should equal radius for circular orbit
    @test elements.a ≈ r rtol=1e-4

    # Eccentricity should be ~0 for circular orbit
    @test elements.e < 0.01

    # Inclination should be ~0 for equatorial orbit
    @test elements.i ≈ 0.0u"rad" atol=0.1u"rad"
end

@testset "Orbital Elements - Elliptical Orbit" begin
    # Elliptical orbit
    h_periapsis = 400.0u"km"
    h_apoapsis = 800.0u"km"

    r_p = 6.3781370e6u"m" + h_periapsis
    r_a = 6.3781370e6u"m" + h_apoapsis

    a = (r_p + r_a) / 2  # Semi-major axis

    # Velocity at periapsis
    μ = 3.986004418e14u"m^3/s^2"
    v_p = sqrt(μ * (2/r_p - 1/a))

    position = SVector(r_p, 0.0u"m", 0.0u"m")
    velocity = SVector(0.0u"m/s", v_p, 0.0u"m/s")

    elements = orbital_elements(position, velocity)

    # Semi-major axis
    @test elements.a ≈ a rtol=0.01

    # Eccentricity
    e_expected = (r_a - r_p) / (r_a + r_p)
    @test elements.e ≈ e_expected rtol=0.01

    # Should be elliptical
    @test 0.0 < elements.e < 1.0
end

@testset "Physical Reasonableness" begin
    # Gravity should decrease with altitude
    a_0 = abs(gravity_point_mass(6.378e6u"m"))
    a_400 = abs(gravity_point_mass(6.778e6u"m"))  # +400 km
    @test a_0 > a_400

    # Orbital velocity should decrease with altitude
    v_200 = orbital_velocity(200.0u"km")
    v_400 = orbital_velocity(400.0u"km")
    v_800 = orbital_velocity(800.0u"km")
    @test v_200 > v_400 > v_800

    # Escape velocity should decrease with altitude
    v_esc_0 = escape_velocity(0.0u"km")
    v_esc_400 = escape_velocity(400.0u"km")
    @test v_esc_0 > v_esc_400
end

@testset "Conservation Laws" begin
    # Specific energy should be conserved in circular orbit
    h = 400.0u"km"
    r = 6.3781370e6u"m" + h
    v = orbital_velocity(h)

    position1 = SVector(r, 0.0u"m", 0.0u"m")
    velocity1 = SVector(0.0u"m/s", v, 0.0u"m/s")

    # At different position in orbit (90 degrees later)
    position2 = SVector(0.0u"m", r, 0.0u"m")
    velocity2 = SVector(-v, 0.0u"m/s", 0.0u"m/s")

    ε1 = specific_orbital_energy(velocity1, position1)
    ε2 = specific_orbital_energy(velocity2, position2)

    @test ε1 ≈ ε2 rtol=1e-6  # Energy conserved
end

println("  ✓ Gravity model tests passed")
