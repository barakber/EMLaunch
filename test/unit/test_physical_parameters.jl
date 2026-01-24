"""
Tests for physical_parameters.jl

Tests:
- Physical constants have correct values
- Units are properly assigned
- Helper functions work correctly
- No dimensional errors
"""


@testset "Earth Parameters" begin
    # Test gravitational parameter
    @test μ_Earth ≈ 3.986004418e14u"m^3/s^2"
    @test dimension(μ_Earth) == dimension(1.0u"m^3/s^2")

    # Test Earth radius
    @test R_Earth ≈ 6.3781370e6u"m" rtol=1e-6
    @test dimension(R_Earth) == dimension(1.0u"m")

    # Test standard gravity
    @test g_0 ≈ 9.80665u"m/s^2"
    @test dimension(g_0) == dimension(1.0u"m/s^2")

    # Test J2 coefficient (dimensionless)
    @test J2_Earth ≈ 1.08263e-3
    @test typeof(J2_Earth) <: Real
end

@testset "Atmospheric Parameters" begin
    # Sea level conditions
    @test ρ_0 ≈ 1.225u"kg/m^3"
    @test dimension(ρ_0) == dimension(1.0u"kg/m^3")

    @test P_0 ≈ 101325.0u"Pa"
    @test dimension(P_0) == dimension(1.0u"Pa")

    @test T_0 ≈ 288.15u"K"
    @test dimension(T_0) == dimension(1.0u"K")

    # Scale height
    @test H_scale ≈ 8500.0u"m"
    @test dimension(H_scale) == dimension(1.0u"m")
end

@testset "EM Launcher Parameters" begin
    # Target velocities
    @test v_target_orbital ≈ 8000.0u"m/s"
    @test dimension(v_target_orbital) == dimension(1.0u"m/s")

    @test v_target_hypersonic ≈ 2000.0u"m/s"
    @test dimension(v_target_hypersonic) == dimension(1.0u"m/s")

    # Launcher dimensions
    @test L_launcher_nominal ≈ 2000.0u"m"
    @test L_launcher_min < L_launcher_nominal < L_launcher_max

    # Payload mass
    @test m_payload_nominal ≈ 20.0u"kg"
    @test m_payload_min < m_payload_nominal < m_payload_max

    # Acceleration
    @test dimension(a_avg_nominal) == dimension(1.0u"m/s^2")
    @test a_avg_nominal > 0.0u"m/s^2"
end

@testset "Helper Functions" begin
    # Test coil_spacing
    spacing = coil_spacing(200, 2000.0u"m")
    @test spacing ≈ 10.0u"m"
    @test dimension(spacing) == dimension(1.0u"m")

    # Test average_acceleration
    a = average_acceleration(8000.0u"m/s", 2000.0u"m")
    @test a ≈ 16000.0u"m/s^2"
    @test dimension(a) == dimension(1.0u"m/s^2")

    # Test g_loading
    g_load = g_loading(9.80665u"m/s^2")
    @test g_load ≈ 1.0 rtol=1e-10

    g_load_high = g_loading(98.0665u"m/s^2")
    @test g_load_high ≈ 10.0 rtol=1e-6

    # Note: mach_number is now in atmosphere.jl and takes altitude instead of temperature
    # It is tested in test_atmosphere.jl
end

@testset "Unit Consistency" begin
    # Test that energy calculations have correct units
    m = 20.0u"kg"
    v = 8000.0u"m/s"
    KE = 0.5 * m * v^2

    @test dimension(KE) == dimension(1.0u"J")  # Dimensionally equivalent
    @test uconvert(u"MJ", KE) ≈ 640.0u"MJ"

    # Test that force has correct units
    a = 10.0u"m/s^2"
    F = m * a
    @test dimension(F) == dimension(1.0u"N")  # Dimensionally equivalent

    # Test dimensional consistency
    L = 2000.0u"m"
    t = L / v
    @test dimension(t) == dimension(1.0u"s")
    @test t ≈ 0.25u"s"
end

@testset "Physical Constant Ranges" begin
    # Ensure physical values are in reasonable ranges
    @test 6.0e6u"m" < R_Earth < 7.0e6u"m"
    @test 3.9e14u"m^3/s^2" < μ_Earth < 4.0e14u"m^3/s^2"
    @test 9.7u"m/s^2" < g_0 < 9.9u"m/s^2"
    @test 1.0u"kg/m^3" < ρ_0 < 1.3u"kg/m^3"
    @test 100000.0u"Pa" < P_0 < 102000.0u"Pa"
    @test 287.0u"K" < T_0 < 290.0u"K"
end

println("  ✓ Physical parameters tests passed")
