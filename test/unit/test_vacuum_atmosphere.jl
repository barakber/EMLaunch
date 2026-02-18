"""
Tests for atmosphere.jl with vacuum body (Moon)

Tests:
- Zero drag force
- Zero aerodynamic heating
- Zero density and pressure
- Vacuum-safe Mach number
"""

@testset "Vacuum Density" begin
    h = 0.0u"km"
    ρ = density(h; body=MOON)
    @test ρ == 0.0u"kg/m^3"

    ρ_exp = density_exponential(h; body=MOON)
    @test ρ_exp == 0.0u"kg/m^3"
end

@testset "Vacuum Pressure" begin
    h = 0.0u"km"
    P = pressure(h; body=MOON)
    @test P == 0.0u"Pa"
end

@testset "Vacuum Temperature" begin
    h = 0.0u"km"
    T = temperature(h; body=MOON)
    @test T == 200.0u"K"  # Constant for vacuum bodies
end

@testset "Vacuum Speed of Sound" begin
    h = 0.0u"km"
    a = speed_of_sound(h; body=MOON)
    @test isnan(ustrip(a))  # NaN for vacuum
end

@testset "Vacuum Mach Number" begin
    v = 1000.0u"m/s"
    h = 0.0u"km"
    M = mach_number(v, h; body=MOON)
    @test M == 0.0  # Zero Mach for vacuum
end

@testset "Vacuum Drag Force" begin
    # Scalar velocity
    v = 1000.0u"m/s"
    h = 0.0u"km"
    A = 0.1u"m^2"
    F = drag_force(v, h, A; body=MOON)
    @test F == 0.0u"N"

    # Vector velocity
    v_vec = SVector(1000.0u"m/s", 0.0u"m/s", 0.0u"m/s")
    F_vec = drag_force(v_vec, h, A; body=MOON)
    @test norm(F_vec) == 0.0u"N"
end

@testset "Vacuum Aerodynamic Heating" begin
    v = 2000.0u"m/s"
    h = 0.0u"km"
    R_n = 0.1u"m"
    q = aerodynamic_heating(v, h, R_n; body=MOON)
    @test q == 0.0u"W/m^2"
end

@testset "Vacuum Reynolds Number" begin
    v = 1000.0u"m/s"
    h = 0.0u"km"
    L = 1.0u"m"
    Re = reynolds_number(v, h, L; body=MOON)
    @test Re == 0.0
end

@testset "Earth Atmosphere Unchanged" begin
    # Verify Earth functions still work the same with explicit body=EARTH
    h = 0.0u"km"

    T = temperature(h; body=EARTH)
    @test T ≈ 288.15u"K" rtol=1e-6

    P = pressure(h; body=EARTH)
    @test P ≈ 101325.0u"Pa" rtol=1e-4

    ρ = density(h; body=EARTH)
    @test ρ ≈ 1.225u"kg/m^3" rtol=1e-3

    a = speed_of_sound(h; body=EARTH)
    @test a ≈ 340.3u"m/s" rtol=0.01
end

println("  ✓ Vacuum atmosphere tests passed")
