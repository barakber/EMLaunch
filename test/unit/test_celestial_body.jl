"""
Tests for celestial_body.jl

Tests:
- EARTH and MOON constants
- has_atmosphere dispatch
- CelestialBody struct fields
- AtmosphereModel type hierarchy
"""

@testset "EARTH Constants" begin
    @test EARTH.name == "Earth"
    @test EARTH.mu ≈ 3.986004418e14u"m^3/s^2" rtol=1e-6
    @test EARTH.radius ≈ 6.3781370e6u"m" rtol=1e-6
    @test EARTH.J2 ≈ 1.08263e-3 rtol=1e-3
    @test EARTH.rotation_rate ≈ 7.2921159e-5u"rad/s" rtol=1e-6
    @test EARTH.surface_gravity ≈ 9.80665u"m/s^2" rtol=1e-6
    @test EARTH.atmosphere isa StandardAtmosphere1976
end

@testset "MOON Constants" begin
    @test MOON.name == "Moon"
    @test MOON.mu ≈ 4.9048695e12u"m^3/s^2" rtol=1e-3
    @test MOON.radius ≈ 1.7374e6u"m" rtol=1e-3
    @test MOON.J2 ≈ 2.033e-4 rtol=0.01
    @test MOON.surface_gravity ≈ 1.625u"m/s^2" rtol=0.01
    @test MOON.atmosphere isa VacuumAtmosphere
end

@testset "has_atmosphere" begin
    @test has_atmosphere(EARTH) == true
    @test has_atmosphere(MOON) == false
end

@testset "AtmosphereModel Types" begin
    @test VacuumAtmosphere <: AtmosphereModel
    @test StandardAtmosphere1976 <: AtmosphereModel

    # EARTH_ATMOSPHERE should have valid layers
    @test length(EARTH_ATMOSPHERE.layers) == 7
    @test EARTH_ATMOSPHERE.sea_level_density ≈ 1.225u"kg/m^3" rtol=1e-3
    @test EARTH_ATMOSPHERE.sea_level_temperature ≈ 288.15u"K" rtol=1e-6
    @test EARTH_ATMOSPHERE.heat_ratio ≈ 1.4 rtol=1e-6
end

@testset "CelestialBody Construction" begin
    # Create a custom body
    mars = CelestialBody(
        "Mars",
        4.2828e13u"m^3/s^2",
        3.3895e6u"m",
        1.96045e-3,
        7.088e-5u"rad/s",
        3.72076u"m/s^2",
        VacuumAtmosphere()  # Simplified (Mars has thin atmosphere)
    )
    @test mars.name == "Mars"
    @test mars.surface_gravity ≈ 3.72076u"m/s^2" rtol=1e-3
    @test !has_atmosphere(mars)
end

println("  ✓ Celestial body tests passed")
