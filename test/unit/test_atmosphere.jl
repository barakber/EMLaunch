"""
Tests for atmosphere.jl

Tests:
- US Standard Atmosphere 1976 values at known altitudes
- Density, temperature, pressure calculations
- Drag force calculations
- Aerodynamic heating
- Mach number calculations
"""

using StaticArrays

@testset "Sea Level Conditions" begin
    h = 0.0u"km"

    # Temperature at sea level
    T = temperature(h)
    @test T ≈ 288.15u"K" rtol=1e-6

    # Pressure at sea level
    P = pressure(h)
    @test P ≈ 101325.0u"Pa" rtol=1e-4

    # Density at sea level
    ρ = density(h)
    @test ρ ≈ 1.225u"kg/m^3" rtol=1e-3

    # Speed of sound at sea level
    a = speed_of_sound(h)
    @test a ≈ 340.3u"m/s" rtol=0.01
end

@testset "Standard Atmosphere at Known Altitudes" begin
    # 11 km (tropopause)
    h = 11.0u"km"
    T = temperature(h)
    @test T ≈ 216.65u"K" rtol=1e-3

    # 20 km
    h = 20.0u"km"
    T = temperature(h)
    @test T ≈ 216.65u"K" rtol=1e-3  # Isothermal layer

    # 25 km
    h = 25.0u"km"
    T = temperature(h)
    @test T ≈ 221.65u"K" rtol=0.01

    ρ = density(h)
    @test ρ ≈ 0.04u"kg/m^3" rtol=0.1  # Order of magnitude

    # 50 km
    h = 50.0u"km"
    T = temperature(h)
    @test T > 200.0u"K"
    @test T < 300.0u"K"

    # Density should decrease with altitude
    ρ_0 = density(0.0u"km")
    ρ_10 = density(10.0u"km")
    ρ_25 = density(25.0u"km")
    ρ_50 = density(50.0u"km")

    @test ρ_0 > ρ_10 > ρ_25 > ρ_50
end

@testset "Exponential Atmosphere Model" begin
    # Test simplified exponential model
    h = 8500.0u"m"  # One scale height
    ρ = density_exponential(h)

    @test ρ ≈ 1.225u"kg/m^3" * exp(-1) rtol=0.01
    @test ρ < density_exponential(0.0u"m")
    @test ρ > density_exponential(17000.0u"m")
end

@testset "Speed of Sound" begin
    # Sea level
    a_0 = speed_of_sound(0.0u"km")
    @test a_0 ≈ 340.3u"m/s" rtol=0.01

    # High altitude (colder, slower)
    a_11 = speed_of_sound(11.0u"km")
    @test a_11 < a_0
    @test a_11 ≈ 295.0u"m/s" rtol=0.02
end

@testset "Mach Number" begin
    # Mach 1 at sea level
    v = 340.0u"m/s"
    M = mach_number(v, 0.0u"km")
    @test M ≈ 1.0 rtol=0.01

    # Mach 6 at altitude
    v = 1800.0u"m/s"
    M = mach_number(v, 25.0u"km")
    @test M ≈ 6.0 rtol=0.1

    # Hypersonic regime
    v = 8000.0u"m/s"
    M = mach_number(v, 50.0u"km")
    @test M > 20.0  # Highly hypersonic
end

@testset "Drag Coefficient" begin
    # Subsonic
    C_d_sub = drag_coefficient(0.5, regime=:subsonic)
    @test C_d_sub ≈ 0.4 rtol=0.01

    # Supersonic
    C_d_sup = drag_coefficient(3.0, regime=:supersonic)
    @test C_d_sup ≈ 0.85 rtol=0.01

    # Hypersonic
    C_d_hyp = drag_coefficient(10.0, regime=:hypersonic)
    @test C_d_hyp ≈ 0.90 rtol=0.01

    # Auto-detect regime
    C_d_auto_sub = drag_coefficient(0.5, regime=:auto)
    @test C_d_auto_sub ≈ 0.4 rtol=0.05

    C_d_auto_sup = drag_coefficient(3.0, regime=:auto)
    @test C_d_auto_sup ≈ 0.85 rtol=0.05
end

@testset "Drag Force" begin
    # Sea level, low speed
    v = 100.0u"m/s"
    h = 0.0u"km"
    A = 1.0u"m^2"

    F_d = drag_force(v, h, A)
    @test dimension(F_d) == dimension(1.0u"N")
    @test F_d > 0.0u"N"

    # Higher velocity → higher drag (quadratic)
    F_d_200 = drag_force(200.0u"m/s", h, A)
    @test F_d_200 ≈ 4.0 * F_d rtol=0.1  # v² scaling

    # Higher altitude → lower drag (lower density)
    F_d_high = drag_force(v, 25.0u"km", A)
    @test F_d_high < F_d

    # Vector velocity
    v_vec = SVector(100.0u"m/s", 0.0u"m/s", 0.0u"m/s")
    F_d_vec = drag_force(v_vec, h, A)
    @test norm(F_d_vec) ≈ F_d rtol=0.1
end

@testset "Aerodynamic Heating" begin
    # Low speed → low heating
    v = 100.0u"m/s"
    h = 25.0u"km"
    R_n = 0.1u"m"

    q_dot = aerodynamic_heating(v, h, R_n)
    @test dimension(q_dot) == dimension(1.0u"W/m^2")
    @test q_dot > 0.0u"W/m^2"

    # Hypersonic heating (v³ scaling)
    v_hyper = 2000.0u"m/s"
    q_dot_hyper = aerodynamic_heating(v_hyper, h, R_n)
    @test q_dot_hyper > q_dot
    @test uconvert(u"MW/m^2", q_dot_hyper) > 0.5u"MW/m^2"  # Significant heating

    # Orbital velocity heating
    v_orbital = 8000.0u"m/s"
    q_dot_orbital = aerodynamic_heating(v_orbital, 50.0u"km", R_n)
    @test uconvert(u"MW/m^2", q_dot_orbital) > 5.0u"MW/m^2"  # Very high heating

    # Smaller nose radius → higher heating
    q_dot_small = aerodynamic_heating(v_hyper, h, 0.05u"m")
    @test q_dot_small > q_dot_hyper
end

@testset "Reynolds Number" begin
    v = 100.0u"m/s"
    h = 0.0u"km"
    L = 1.0u"m"

    Re = reynolds_number(v, h, L)
    @test Re > 0.0
    @test Re ≈ 6.8e6 rtol=0.1  # Order of magnitude

    # Higher velocity → higher Re
    Re_high = reynolds_number(1000.0u"m/s", h, L)
    @test Re_high > Re

    # Higher altitude (lower density) → lower Re
    Re_alt = reynolds_number(v, 25.0u"km", L)
    @test Re_alt < Re
end

@testset "Atmospheric Property Monotonicity" begin
    # Density should decrease monotonically with altitude
    altitudes = [0.0, 10.0, 20.0, 30.0, 50.0, 70.0, 100.0] .* 1.0u"km"
    densities = density.(altitudes)

    for i in 1:(length(densities)-1)
        @test densities[i] > densities[i+1]
    end

    # Pressure should decrease monotonically
    pressures = pressure.(altitudes)
    for i in 1:(length(pressures)-1)
        @test pressures[i] > pressures[i+1]
    end
end

@testset "Physical Reasonableness" begin
    # Temperature should be within reasonable bounds
    for h in [0.0, 25.0, 50.0, 75.0, 100.0] .* 1.0u"km"
        T = temperature(h)
        @test 150.0u"K" < T < 350.0u"K"  # Atmospheric range
    end

    # Density should be positive
    for h in [0.0, 25.0, 50.0, 100.0] .* 1.0u"km"
        ρ = density(h)
        @test ρ > 0.0u"kg/m^3"
    end

    # Pressure should be positive
    for h in [0.0, 25.0, 50.0, 100.0] .* 1.0u"km"
        P = pressure(h)
        @test P > 0.0u"Pa"
    end
end

println("  ✓ Atmospheric model tests passed")
