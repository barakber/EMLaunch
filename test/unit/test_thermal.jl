"""
Thermal Protection & Ablation Tests

Tests ablative heat shield modeling, thermal survival, and heat shield sizing.
"""

using Test
using Unitful
using Printf

# Load thermal module

println()
println("=" ^ 80)
println("THERMAL PROTECTION & ABLATION TESTS")
println("=" ^ 80)
println()

@testset "Ablative Materials" begin
    println("\n[Test 1] Ablative material properties")

    # Check all materials are defined
    @test haskey(ABLATIVE_MATERIALS, :carbon_phenolic)
    @test haskey(ABLATIVE_MATERIALS, :pica)
    @test haskey(ABLATIVE_MATERIALS, :avcoat)
    @test haskey(ABLATIVE_MATERIALS, :cork)

    # Check PICA properties (NASA's high-performance material)
    pica = ABLATIVE_MATERIALS[:pica]
    @test pica.name == "PICA (Phenolic Impregnated Carbon Ablator)"
    @test pica.H_abl ≈ 12.5e6u"J/kg"
    @test pica.density ≈ 250.0u"kg/m^3"  # Lightweight
    @test pica.T_max ≈ 3000.0u"K"
    @test pica.emissivity ≈ 0.90

    println("  ✓ All ablative materials defined")
    println("  ✓ PICA properties validated")
end

@testset "Heat Shield Creation" begin
    println("\n[Test 2] Heat shield creation and initialization")

    # Create PICA heat shield
    thickness = 0.05u"m"  # 5 cm
    area = 1.0u"m^2"
    shield = HeatShield(:pica, thickness, area)

    @test shield.material.name == "PICA (Phenolic Impregnated Carbon Ablator)"
    @test shield.initial_thickness ≈ 0.05u"m"
    @test shield.thickness ≈ 0.05u"m"
    @test shield.area ≈ 1.0u"m^2"

    # Check mass calculation
    expected_mass = 250.0u"kg/m^3" * 0.05u"m" * 1.0u"m^2"
    @test shield.mass ≈ expected_mass
    @test shield.total_heat_absorbed ≈ 0.0u"J"

    println("  ✓ Heat shield created successfully")
    mass_kg = ustrip(u"kg", shield.mass)
    println("  ✓ Mass calculation correct: $mass_kg kg")
end

@testset "Ablation Rate Calculations" begin
    println("\n[Test 3] Ablation rate calculations")

    material = ABLATIVE_MATERIALS[:pica]

    # Test at moderate heat flux (similar to Space Shuttle reentry)
    q_dot = 1.0e6u"W/m^2"  # 1 MW/m²

    m_dot = ablation_rate(q_dot, material)
    h_dot = ablation_thickness_rate(q_dot, material)

    # Ablation rate = q / H_abl
    expected_m_dot = q_dot / material.H_abl
    @test m_dot ≈ expected_m_dot

    # Thickness rate = m_dot / ρ
    expected_h_dot = m_dot / material.density
    @test h_dot ≈ expected_h_dot

    m_dot_val = round(ustrip(u"kg/(m^2*s)", m_dot), digits=6)
    h_dot_val = round(ustrip(u"m/s", h_dot)*1000, digits=3)
    println("  ✓ Ablation rate at 1 MW/m²: $m_dot_val kg/(m²·s)")
    println("  ✓ Thickness loss rate: $h_dot_val mm/s")

    # Zero heat flux should give zero ablation
    @test ablation_rate(0.0u"W/m^2", material) ≈ 0.0u"kg/(m^2*s)"
    @test ablation_rate(-100.0u"W/m^2", material) ≈ 0.0u"kg/(m^2*s)"
end

@testset "Heat Shield Sizing" begin
    println("\n[Test 4] Heat shield sizing for missions")

    # Apollo reentry: ~3 GW total, ~100 s duration → ~50 MJ/m² total heat load
    apollo_heat_load = 50.0e6u"J/m^2"  # 50 MJ/m²
    area = 10.0u"m^2"  # Rough Apollo heat shield area

    # Carbon phenolic (Apollo used this)
    mass_required = required_heat_shield_mass(apollo_heat_load, area, :carbon_phenolic)
    thickness_required = required_heat_shield_thickness(apollo_heat_load, :carbon_phenolic)

    println("  Apollo-class reentry (50 MJ/m²):")
    mass_val = round(ustrip(u"kg", mass_required), digits=1)
    thick_val = round(ustrip(u"cm", thickness_required), digits=1)
    println("    Required mass: $mass_val kg")
    println("    Required thickness: $thick_val cm")

    # Check reasonableness (should be on order of mm to cm)
    # Apollo actually had 2-5 cm, but this includes additional safety margins
    @test thickness_required > 0.001u"m"  # At least 1 mm
    @test thickness_required < 0.20u"m"  # Less than 20 cm

    # Mars entry (Curiosity/Perseverance PICA heat shield): ~100 MJ/m²
    mars_heat_load = 100.0e6u"J/m^2"
    thickness_mars = required_heat_shield_thickness(mars_heat_load, :pica)

    println("  Mars entry (100 MJ/m²) with PICA:")
    thick_mars_val = round(ustrip(u"cm", thickness_mars), digits=1)
    println("    Required thickness: $thick_mars_val cm")

    @test thickness_mars > 0.02u"m"  # > 2 cm
    @test thickness_mars < 0.50u"m"  # < 50 cm
end

@testset "Heat Shield Dynamics" begin
    println("\n[Test 5] Heat shield ablation over time")

    # Create PICA shield
    shield = HeatShield(:pica, 0.10u"m", 1.0u"m^2")  # 10 cm thick

    initial_mass = shield.mass
    initial_thickness = shield.thickness

    # Simulate high heat flux for 10 seconds
    q_dot = 2.0e6u"W/m^2"  # 2 MW/m²
    dt = 0.1u"s"
    n_steps = 100  # Total 10 seconds

    total_mass_lost = 0.0u"kg"
    total_heat = 0.0u"J"

    for i in 1:n_steps
        result = update_heat_shield!(shield, q_dot, dt)
        total_mass_lost += result.mass_lost
        total_heat += result.heat_absorbed
    end

    println("  10 s exposure at 2 MW/m²:")
    init_thick = round(ustrip(u"cm", initial_thickness), digits=2)
    final_thick = round(ustrip(u"cm", shield.thickness), digits=2)
    lost_thick = round(ustrip(u"mm", initial_thickness - shield.thickness), digits=2)
    lost_mass = round(ustrip(u"kg", total_mass_lost), digits=3)
    heat_absorbed = round(ustrip(u"MJ", total_heat), digits=1)
    margin_pct = round(heat_shield_margin(shield) * 100, digits=1)
    println("    Initial thickness: $init_thick cm")
    println("    Final thickness: $final_thick cm")
    println("    Thickness lost: $lost_thick mm")
    println("    Mass lost: $lost_mass kg")
    println("    Heat absorbed: $heat_absorbed MJ")
    println("    Shield margin: $margin_pct%")

    # Check shield is still intact
    @test heat_shield_intact(shield)
    @test shield.thickness < initial_thickness
    @test shield.mass < initial_mass
    @test shield.total_heat_absorbed ≈ total_heat

    # Margin should be reasonable
    margin = heat_shield_margin(shield)
    @test margin > 0.5  # Should have >50% left for this short duration
    @test margin < 1.0  # Should have ablated some
end

@testset "Thermal Survival Criteria" begin
    println("\n[Test 6] Thermal survival checks")

    # Case 1: Intact heat shield, normal temperature
    shield = HeatShield(:pica, 0.05u"m", 1.0u"m^2")
    q_dot = 1.0e6u"W/m^2"
    T = 1500.0u"K"

    status = thermal_survival_check(q_dot, T, shield, 1.0u"s")

    @test status.survived == true
    @test status.shield_intact == true
    @test status.temp_ok == true
    @test status.margin ≈ 1.0
    @test isnothing(status.failure_mode)

    println("  ✓ Case 1: Intact shield, normal temp → SURVIVED")

    # Case 2: Heat shield burned through
    shield_thin = HeatShield(:pica, 0.001u"m", 1.0u"m^2")  # Only 1 mm
    shield_thin.thickness = 0.0u"m"  # Simulate burn-through

    status = thermal_survival_check(q_dot, T, shield_thin, 1.0u"s")

    @test status.survived == false
    @test status.shield_intact == false
    @test status.failure_mode == "Heat shield burned through"

    println("  ✓ Case 2: Burned through → FAILED")

    # Case 3: Temperature exceeded
    shield_hot = HeatShield(:pica, 0.05u"m", 1.0u"m^2")
    T_hot = 3500.0u"K"  # Above PICA's 3000 K limit

    status = thermal_survival_check(q_dot, T_hot, shield_hot, 1.0u"s")

    @test status.survived == false
    @test status.temp_ok == false

    println("  ✓ Case 3: Temperature exceeded → FAILED")

    # Case 4: No heat shield (radiative cooling only)
    status = thermal_survival_check(q_dot, 1200.0u"K", nothing, 1.0u"s")

    @test status.shield_intact == false
    @test status.temp_ok == true  # Below 1500 K limit

    println("  ✓ Case 4: No shield, low temp → SURVIVED (radiative)")

    status_hot = thermal_survival_check(q_dot, 2000.0u"K", nothing, 1.0u"s")
    @test status_hot.survived == false

    println("  ✓ Case 5: No shield, high temp → FAILED")
end

@testset "EM Launcher Heat Loads" begin
    println("\n[Test 7] Heat shield sizing for EM launcher exit")

    # From LEO analysis: 7800 m/s at sea level → ~3-4 GW/m² for ~2 seconds
    # That's ~6-8 GJ/m² total heat load!

    println("\n  EM Launcher Hypersonic Exit Scenarios:")
    println("  " * "-"^60)

    # Scenario 1: Sea level exit (worst case)
    v = 7800.0u"m/s"
    h = 0.0u"m"
    R_nose = 0.2u"m"

    q_dot_sea = aerodynamic_heating(v, h, R_nose)

    println("  Sea Level Exit (v=7800 m/s):")
    q_sea_val = round(ustrip(u"GW/m^2", q_dot_sea), digits=2)
    println("    Heat flux: $q_sea_val GW/m²")

    # Assume 2 seconds in atmosphere before crash
    duration = 2.0u"s"
    heat_load_sea = q_dot_sea * duration

    println("    Duration: 2 s (before crash)")
    heat_sea_val = round(ustrip(u"GJ/m^2", heat_load_sea), digits=1)
    println("    Total heat load: $heat_sea_val GJ/m²")

    # Calculate required shield
    thickness_sea = required_heat_shield_thickness(heat_load_sea, :pica)
    mass_sea = required_heat_shield_mass(heat_load_sea, 1.0u"m^2", :pica)

    thick_sea_val = round(ustrip(u"m", thickness_sea), digits=2)
    mass_sea_val = round(ustrip(u"kg", mass_sea), digits=1)
    println("    Required PICA thickness: $thick_sea_val m")
    println("    Required PICA mass (1 m²): $mass_sea_val kg")

    # This should be HUGE (meters thick!)
    @test thickness_sea > 0.1u"m"  # >10 cm minimum

    # Scenario 2: 20 km altitude exit (10x less drag)
    h_20km = 20000.0u"m"
    q_dot_20km = aerodynamic_heating(v, h_20km, R_nose)

    println("\n  20 km Altitude Exit (v=7800 m/s):")
    q_20km_val = round(ustrip(u"MW/m^2", q_dot_20km), digits=1)
    println("    Heat flux: $q_20km_val MW/m²")

    # Assume 30 seconds to get through remaining atmosphere
    duration_20km = 30.0u"s"
    heat_load_20km = q_dot_20km * duration_20km

    println("    Duration: 30 s")
    heat_20km_val = round(ustrip(u"MJ/m^2", heat_load_20km), digits=1)
    println("    Total heat load: $heat_20km_val MJ/m²")

    thickness_20km = required_heat_shield_thickness(heat_load_20km, :pica)
    mass_20km = required_heat_shield_mass(heat_load_20km, 1.0u"m^2", :pica)

    thick_20km_val = round(ustrip(u"cm", thickness_20km), digits=1)
    mass_20km_val = round(ustrip(u"kg", mass_20km), digits=2)
    println("    Required PICA thickness: $thick_20km_val cm")
    println("    Required PICA mass (1 m²): $mass_20km_val kg")

    # Key insight: 20 km requires MORE total heat shield because of longer duration!
    # Even though heat flux is lower, total heat load is higher (1.7 GJ vs 0.4 GJ)
    @test thickness_20km < 2.0u"m"  # Less than 2 meters
    @test mass_20km > 50.0u"kg"  # Significant mass due to longer duration

    # Heat flux ratio should match density reduction (~10x), but total load is different
    heat_flux_ratio = q_dot_sea / q_dot_20km  # Keep units for safety
    heat_flux_ratio_val = round(ustrip(heat_flux_ratio), digits=1)
    println("\n  Heat flux ratio (sea level / 20 km): $(heat_flux_ratio_val)x")
    @test q_dot_sea > 2.0 * q_dot_20km  # Should be significantly higher at sea level

    # Total heat load ratio (considering duration)
    total_heat_ratio = heat_load_sea / heat_load_20km  # Unit-aware
    total_heat_ratio_val = round(ustrip(total_heat_ratio), digits=2)
    println("  Total heat load ratio: $(total_heat_ratio_val)x")
    println("  → 20 km needs MORE shield due to longer exposure (30s vs 2s)")
end

@testset "Physical Validation" begin
    println("\n[Test 8] Physical reasonableness checks")

    # Heat of ablation should be positive
    for (key, mat) in ABLATIVE_MATERIALS
        @test mat.H_abl > 0.0u"J/kg"
        @test mat.density > 0.0u"kg/m^3"
        @test mat.T_max > 0.0u"K"
        @test 0.0 <= mat.emissivity <= 1.0
        @test 0.0 <= mat.char_fraction <= 1.0
    end

    # Ablation rate should increase with heat flux
    mat = ABLATIVE_MATERIALS[:pica]
    m_dot_1 = ablation_rate(1.0e6u"W/m^2", mat)
    m_dot_2 = ablation_rate(2.0e6u"W/m^2", mat)
    @test m_dot_2 ≈ 2.0 * m_dot_1

    # Heat shield mass should scale with area
    heat_load = 100.0e6u"J/m^2"
    mass_1m2 = required_heat_shield_mass(heat_load, 1.0u"m^2", :pica)
    mass_10m2 = required_heat_shield_mass(heat_load, 10.0u"m^2", :pica)
    @test mass_10m2 ≈ 10.0 * mass_1m2

    println("  ✓ All physical properties valid")
    println("  ✓ Ablation rate scales linearly with heat flux")
    println("  ✓ Heat shield mass scales with area")
end

println()
println("=" ^ 80)
println("THERMAL PROTECTION TESTS COMPLETE")
println("=" ^ 80)
println()
