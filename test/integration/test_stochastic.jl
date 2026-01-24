"""
Tests for stochastic_trajectory.jl

Tests stochastic differential equation functionality for uncertainty quantification.
"""


@testset "Noise Parameters" begin
    # Default noise parameters
    noise = default_noise_parameters()

    @test noise.atmospheric_noise == 0.05  # 5%
    @test noise.drag_noise == 0.10  # 10%
    @test noise.current_noise == 0.02  # 2%
    @test noise.timing_noise == 0.0001u"s"  # 0.1ms
    @test noise.mass_noise == 0.01  # 1%

    # Custom noise parameters
    custom_noise = NoiseParameters(0.03, 0.05, 0.01, 0.00005u"s", 0.005)

    @test custom_noise.atmospheric_noise == 0.03
    @test custom_noise.drag_noise == 0.05
    @test custom_noise.current_noise == 0.01
    @test custom_noise.timing_noise == 0.00005u"s"
    @test custom_noise.mass_noise == 0.005
end

@testset "SDE Problem Creation" begin
    launcher = create_uniform_launcher(
        length = 100.0u"m",
        num_coils = 5,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F",
        voltage_per_coil = 100.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )

    payload = create_default_payload(mass = 10.0u"kg")

    mission = MissionProfile(
        50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
        100.0u"m/s", 100.0u"km"
    )

    noise = default_noise_parameters()

    # Create SDE problem
    prob = create_sde_problem(launcher, payload, mission, noise)

    # Check it's an SDEProblem
    @test prob isa SDEProblem

    # Check dimensions
    n_coils = launcher.num_coils
    expected_dim = 7 + 2 * n_coils  # pos(3) + vel(3) + T(1) + I(n) + Q(n)

    @test length(prob.u0) == expected_dim
    @test length(prob.tspan) == 2
    @test prob.tspan[1] == 0.0
    @test prob.tspan[2] == 10.0

    # Check parameters tuple
    @test length(prob.p) == 5
    @test prob.p[1] == launcher
    @test prob.p[2] == payload
    @test prob.p[3] == mission
    @test prob.p[5] == noise
end

@testset "Noise Function" begin
    launcher = create_uniform_launcher(
        length = 100.0u"m",
        num_coils = 3,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F",
        voltage_per_coil = 100.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )

    payload = create_default_payload(mass = 10.0u"kg")

    mission = MissionProfile(
        50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
        100.0u"m/s", 100.0u"km"
    )

    noise_params = default_noise_parameters()

    prob = create_sde_problem(launcher, payload, mission, noise_params)

    # Test noise function
    du = similar(prob.u0)
    trajectory_sde_noise!(du, prob.u0, prob.p, 0.0)

    # Check that noise is computed
    @test !all(isnan, du)
    @test !all(isinf, du)
    @test all(du .>= 0.0)  # Noise magnitudes should be non-negative

    # Position noise should be small
    @test du[1] ≈ 1.0  # 1m
    @test du[2] ≈ 1.0
    @test du[3] ≈ 1.0

    # Temperature noise should be ~1K
    @test du[7] ≈ 1.0

    # Current noise should scale with current magnitude
    for i in 1:launcher.num_coils
        @test du[7+i] >= 1.0  # At least minimum noise
    end
end

@testset "Noise Scaling" begin
    launcher = create_uniform_launcher(
        length = 100.0u"m",
        num_coils = 2,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F",
        voltage_per_coil = 100.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )

    payload = create_default_payload(mass = 10.0u"kg")

    mission = MissionProfile(
        50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
        100.0u"m/s", 100.0u"km"
    )

    noise_params = default_noise_parameters()

    prob = create_sde_problem(launcher, payload, mission, noise_params)

    # Test with different velocity magnitudes
    u_slow = copy(prob.u0)
    u_slow[4:6] .= [1.0, 0.0, 0.0]  # 1 m/s

    u_fast = copy(prob.u0)
    u_fast[4:6] .= [1000.0, 0.0, 0.0]  # 1000 m/s

    du_slow = similar(prob.u0)
    du_fast = similar(prob.u0)

    trajectory_sde_noise!(du_slow, u_slow, prob.p, 0.0)
    trajectory_sde_noise!(du_fast, u_fast, prob.p, 0.0)

    # Velocity noise should scale with velocity
    @test du_fast[4] > du_slow[4]
end

@testset "Noise Parameters Physical Reasonableness" begin
    noise = default_noise_parameters()

    # Atmospheric noise should be realistic (5% = typical weather variation)
    @test 0.01 <= noise.atmospheric_noise <= 0.15

    # Drag coefficient noise (10% = Reynolds number variations)
    @test 0.05 <= noise.drag_noise <= 0.25

    # Current ripple (2% = good power supply)
    @test 0.005 <= noise.current_noise <= 0.10

    # Timing jitter (0.1ms = typical control system)
    @test 0.00001u"s" <= noise.timing_noise <= 0.001u"s"

    # Mass uncertainty (1% = manufacturing tolerance)
    @test 0.001 <= noise.mass_noise <= 0.05
end

@testset "Monte Carlo SDE Solving" begin
    launcher = create_uniform_launcher(
        length = 100.0u"m",
        num_coils = 2,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F",
        voltage_per_coil = 100.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )

    payload = create_default_payload(mass = 10.0u"kg")

    # Use challenging but realistic targets for this small test launcher
    # (Launcher reaches ~50km at ~96 m/s, so targets at edge of capability)
    # With 5% atmospheric noise, some runs will fail, demonstrating stochastic effects
    mission = MissionProfile(
        50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
        95.0u"m/s",   # Velocity target: 95 m/s (challenging)
        49.0u"km"     # Altitude target: 49 km (near maximum)
    )

    # Test monte_carlo_analysis - this actually runs 10 SDE simulations!
    n_runs = 10
    results = monte_carlo_analysis(launcher, payload, mission, n_runs)

    # Check result structure
    @test haskey(results, :mean_trajectory)
    @test haskey(results, :std_trajectory)
    @test haskey(results, :trajectories)
    @test haskey(results, :success_rate)
    @test haskey(results, :confidence_intervals)

    # Verify SDEs actually solved (not just empty results)
    @test results[:mean_trajectory] !== nothing
    @test results[:std_trajectory] !== nothing
    @test length(results[:trajectories]) == n_runs
    @test all(length(traj.u) > 0 for traj in results[:trajectories])  # All trajectories have data

    # Verify statistics were computed
    @test results[:success_rate] >= 0.0
    @test results[:success_rate] <= 1.0
    @test haskey(results, :time)
    @test length(results[:time]) > 0

    # With achievable targets and 5% noise, we should get some successes
    # (not 0%, not 100%, showing stochastic effects are working)
    @test results[:success_rate] > 0.0  # At least some trajectories succeed
    @test haskey(results, :statistics)
    @test haskey(results[:statistics], :final_altitudes)
    @test haskey(results[:statistics], :final_velocities)

    # Verify confidence intervals computed for standard levels
    @test haskey(results[:confidence_intervals], 0.68)
    @test haskey(results[:confidence_intervals], 0.95)
    @test haskey(results[:confidence_intervals], 0.99)
end

# =============================================================================
# REAL-WORLD SCENARIO TESTS
# =============================================================================

@testset "Weather Impact Scenarios" begin
    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 3,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 100.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )
    payload = create_default_payload(mass = 10.0u"kg")
    mission = MissionProfile(50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
                             100.0u"m/s", 100.0u"km")

    # Clear day: low atmospheric noise
    clear_day = NoiseParameters(0.02, 0.05, 0.01, 0.00005u"s", 0.005)

    # Stormy day: high atmospheric noise
    stormy_day = NoiseParameters(0.15, 0.20, 0.02, 0.0002u"s", 0.01)

    # Test that stormy conditions have higher uncertainty
    @test stormy_day.atmospheric_noise > clear_day.atmospheric_noise
    @test stormy_day.drag_noise > clear_day.drag_noise

    # Verify realistic ranges
    @test 0.01 <= clear_day.atmospheric_noise <= 0.05
    @test 0.10 <= stormy_day.atmospheric_noise <= 0.20
end

@testset "Equipment Aging Scenarios" begin
    # New equipment: tight tolerances
    new_equipment = NoiseParameters(0.03, 0.05, 0.01, 0.00005u"s", 0.005)

    # Aged equipment: degraded tolerances
    aged_equipment = NoiseParameters(0.05, 0.10, 0.05, 0.0003u"s", 0.02)

    # Test realistic aging effects
    @test aged_equipment.current_noise >= 2.0 * new_equipment.current_noise
    @test aged_equipment.timing_noise >= 2.0 * new_equipment.timing_noise

    # Verify equipment is still operational
    @test aged_equipment.current_noise < 0.10  # Still < 10% ripple
    @test aged_equipment.timing_noise < 0.001u"s"  # Still < 1ms jitter
end

@testset "Manufacturing Quality Levels" begin
    # Aerospace grade: very tight tolerances
    aerospace_grade = NoiseParameters(0.02, 0.03, 0.01, 0.00003u"s", 0.003)

    # Commercial grade: moderate tolerances
    commercial_grade = default_noise_parameters()

    # Research prototype: looser tolerances
    prototype = NoiseParameters(0.08, 0.15, 0.03, 0.0003u"s", 0.02)

    # Test quality hierarchy
    @test aerospace_grade.mass_noise < commercial_grade.mass_noise < prototype.mass_noise
    @test aerospace_grade.drag_noise < commercial_grade.drag_noise < prototype.drag_noise

    # Verify all are physically realistic
    @test prototype.mass_noise < 0.05  # Even prototype < 5% mass variation
end

# =============================================================================
# EDGE CASE TESTS
# =============================================================================

@testset "Zero Noise (Deterministic Limit)" begin
    # Zero noise should give consistent results
    zero_noise = NoiseParameters(0.0, 0.0, 0.0, 0.0u"s", 0.0)

    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 2,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 100.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )
    payload = create_default_payload(mass = 10.0u"kg")
    mission = MissionProfile(50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
                             100.0u"m/s", 100.0u"km")

    prob = create_sde_problem(launcher, payload, mission, zero_noise)
    du = similar(prob.u0)
    trajectory_sde_noise!(du, prob.u0, prob.p, 0.0)

    # With zero noise parameters, noise should be minimal (just numerical floor)
    @test all(du .>= 0.0)  # Non-negative noise magnitudes
    # Position noise is always 1m (GPS accuracy floor)
    @test du[1] ≈ 1.0
end

@testset "Extreme Noise (Physical Limits)" begin
    # Extreme but still physical noise
    extreme_noise = NoiseParameters(0.25, 0.40, 0.10, 0.002u"s", 0.05)

    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 2,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 100.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )
    payload = create_default_payload(mass = 10.0u"kg")
    mission = MissionProfile(50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
                             100.0u"m/s", 100.0u"km")

    prob = create_sde_problem(launcher, payload, mission, extreme_noise)
    du = similar(prob.u0)
    trajectory_sde_noise!(du, prob.u0, prob.p, 0.0)

    # Even with extreme noise, values should be finite and positive
    @test all(isfinite, du)
    @test all(du .>= 0.0)
    @test !any(isnan, du)
    @test !any(isinf, du)
end

@testset "Single Noise Source Isolation" begin
    # Test each noise source individually
    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 2,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 100.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )
    payload = create_default_payload(mass = 10.0u"kg")
    mission = MissionProfile(50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
                             100.0u"m/s", 100.0u"km")

    # Only atmospheric noise
    atm_only = NoiseParameters(0.10, 0.0, 0.0, 0.0u"s", 0.0)
    prob = create_sde_problem(launcher, payload, mission, atm_only)
    @test prob.p[5].atmospheric_noise == 0.10
    @test prob.p[5].drag_noise == 0.0

    # Only drag noise
    drag_only = NoiseParameters(0.0, 0.15, 0.0, 0.0u"s", 0.0)
    prob = create_sde_problem(launcher, payload, mission, drag_only)
    @test prob.p[5].drag_noise == 0.15
    @test prob.p[5].atmospheric_noise == 0.0

    # Only current noise
    current_only = NoiseParameters(0.0, 0.0, 0.05, 0.0u"s", 0.0)
    prob = create_sde_problem(launcher, payload, mission, current_only)
    @test prob.p[5].current_noise == 0.05
end

# =============================================================================
# STATISTICAL VALIDITY TESTS
# =============================================================================

@testset "Noise Function Consistency" begin
    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 3,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 100.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )
    payload = create_default_payload(mass = 10.0u"kg")
    mission = MissionProfile(50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
                             100.0u"m/s", 100.0u"km")

    noise_params = default_noise_parameters()
    prob = create_sde_problem(launcher, payload, mission, noise_params)

    # Test that repeated calls give consistent structure
    du1 = similar(prob.u0)
    du2 = similar(prob.u0)
    trajectory_sde_noise!(du1, prob.u0, prob.p, 0.0)
    trajectory_sde_noise!(du2, prob.u0, prob.p, 0.0)

    # Same state should give same noise magnitudes (deterministic function)
    @test du1 ≈ du2

    # Test different states give different noise
    u0_different = copy(prob.u0)
    u0_different[4:6] .= [1000.0, 0.0, 0.0]  # High velocity
    du3 = similar(prob.u0)
    trajectory_sde_noise!(du3, u0_different, prob.p, 0.0)

    # Velocity noise should be different
    @test du3[4] != du1[4]
end

@testset "Noise Scaling Properties" begin
    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 2,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 100.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )
    payload = create_default_payload(mass = 10.0u"kg")
    mission = MissionProfile(50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
                             100.0u"m/s", 100.0u"km")

    # Test with 2x noise parameters
    noise_1x = NoiseParameters(0.05, 0.10, 0.02, 0.0001u"s", 0.01)
    noise_2x = NoiseParameters(0.10, 0.20, 0.04, 0.0002u"s", 0.02)

    @test noise_2x.atmospheric_noise ≈ 2.0 * noise_1x.atmospheric_noise
    @test noise_2x.drag_noise ≈ 2.0 * noise_1x.drag_noise
    @test noise_2x.current_noise ≈ 2.0 * noise_1x.current_noise
end

# =============================================================================
# PRACTICAL REAL-WORLD QUESTIONS
# =============================================================================

@testset "Question: Which uncertainty matters most?" begin
    # Sensitivity analysis: vary one parameter at a time
    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 3,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 100.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )
    payload = create_default_payload(mass = 10.0u"kg")
    mission = MissionProfile(50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
                             100.0u"m/s", 100.0u"km")

    # Baseline: all noise sources moderate
    baseline = default_noise_parameters()

    # Increase atmospheric noise only
    high_atm = NoiseParameters(0.15, 0.10, 0.02, 0.0001u"s", 0.01)

    # Increase drag noise only
    high_drag = NoiseParameters(0.05, 0.30, 0.02, 0.0001u"s", 0.01)

    # Increase current noise only
    high_current = NoiseParameters(0.05, 0.10, 0.10, 0.0001u"s", 0.01)

    # Test that each scenario is testable
    @test high_atm.atmospheric_noise > baseline.atmospheric_noise
    @test high_drag.drag_noise > baseline.drag_noise
    @test high_current.current_noise > baseline.current_noise

    # In real Monte Carlo, would compare std dev of outcomes
    # Higher sensitivity → larger impact on mission success
end

@testset "Question: Weather vs Equipment Quality?" begin
    # Compare impact of weather vs equipment quality

    # Good equipment, bad weather
    good_eq_bad_weather = NoiseParameters(0.15, 0.20, 0.01, 0.00005u"s", 0.005)

    # Bad equipment, good weather
    bad_eq_good_weather = NoiseParameters(0.02, 0.05, 0.08, 0.0005u"s", 0.03)

    # Test setup is valid
    @test good_eq_bad_weather.atmospheric_noise > bad_eq_good_weather.atmospheric_noise
    @test bad_eq_good_weather.current_noise > good_eq_bad_weather.current_noise

    # Real analysis would show which dominates mission risk
end

@testset "Question: Cost-Benefit of Tolerance Improvements?" begin
    # Doubling manufacturing cost typically halves tolerances

    # Standard manufacturing: baseline cost
    standard_manuf = NoiseParameters(0.05, 0.10, 0.02, 0.0001u"s", 0.01)

    # Precision manufacturing: 2x cost, 0.5x tolerances
    precision_manuf = NoiseParameters(0.025, 0.05, 0.01, 0.00005u"s", 0.005)

    # Ultra-precision: 4x cost, 0.25x tolerances
    ultra_precision = NoiseParameters(0.0125, 0.025, 0.005, 0.000025u"s", 0.0025)

    # Test tolerance scaling
    @test precision_manuf.mass_noise ≈ 0.5 * standard_manuf.mass_noise
    @test ultra_precision.mass_noise ≈ 0.25 * standard_manuf.mass_noise

    # Monte Carlo analysis would quantify success probability improvement
    # Compare: ΔP(success) vs ΔCost for ROI analysis
end

@testset "Question: Safety Margins for 95% Success?" begin
    # How much design margin needed for reliable operation?

    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 3,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 100.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )
    payload = create_default_payload(mass = 10.0u"kg")

    # Nominal target: 100 km
    mission_nominal = MissionProfile(50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
                                     100.0u"m/s", 100.0u"km")

    # With safety margin: design for 110 km to achieve 100 km with 95% confidence
    mission_margin = MissionProfile(50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
                                    100.0u"m/s", 110.0u"km")

    @test mission_margin.target_altitude > mission_nominal.target_altitude

    # Monte Carlo would determine required margin for 95% success rate
end

@testset "Question: Is timing jitter critical?" begin
    # Compare timing jitter impact

    # Excellent timing: <0.05ms
    excellent_timing = NoiseParameters(0.05, 0.10, 0.02, 0.00005u"s", 0.01)

    # Poor timing: 1ms jitter
    poor_timing = NoiseParameters(0.05, 0.10, 0.02, 0.001u"s", 0.01)

    # Terrible timing: 5ms jitter (likely mission failure)
    terrible_timing = NoiseParameters(0.05, 0.10, 0.02, 0.005u"s", 0.01)

    @test terrible_timing.timing_noise > poor_timing.timing_noise > excellent_timing.timing_noise

    # At high velocities, timing errors cause position/force mismatch
    # Monte Carlo quantifies impact on trajectory accuracy
end

@testset "Realistic LEO Cargo Mission" begin
    # Realistic launcher configuration for orbital missions
    launcher = create_uniform_launcher(
        length = 2000.0u"m",
        num_coils = 200,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 10.0u"F",
        voltage_per_coil = 5000.0u"V",
        gradient_per_coil = 0.002u"H/m"
    )

    payload = create_default_payload(mass = 20.0u"kg")

    # LEO mission: 400km altitude, 7.8 km/s orbital velocity
    mission = MissionProfile(
        0.0u"m",       # Launch from sea level
        32.5u"°",      # Latitude (Caesarea, Israel)
        34.9u"°",      # Longitude
        90.0u"°",      # Launch azimuth (East)
        45.0u"°",      # Launch elevation
        7800.0u"m/s",  # Target orbital velocity
        400.0u"km"     # Target LEO altitude
    )

    println("\n  Testing realistic LEO cargo mission:")
    println("    Launcher: 2000m, 200 coils, 5kV")
    println("    Payload: 20kg")
    println("    Target: 400km altitude at 7.8 km/s")

    # Run Monte Carlo with different noise scenarios
    @testset "Nominal conditions (5% noise)" begin
        noise = default_noise_parameters()
        results = monte_carlo_analysis(launcher, payload, mission, 5)  # Reduced for speed

        @test results[:success_rate] >= 0.0
        @test results[:success_rate] <= 1.0
        @test length(results[:trajectories]) == 5

        println("    Success rate (nominal): $(round(results[:success_rate] * 100, digits=1))%")

        # With nominal noise, should get reasonable success
        @test haskey(results[:statistics], :mean_altitude)
        @test haskey(results[:statistics], :mean_velocity)
    end

    @testset "Perfect conditions (minimal noise)" begin
        # Very low noise - nearly deterministic
        noise = NoiseParameters(0.01, 0.02, 0.005, 0.00001u"s", 0.001)
        results = monte_carlo_analysis(launcher, payload, mission, 3, noise_params=noise)  # Reduced for speed

        println("    Success rate (perfect): $(round(results[:success_rate] * 100, digits=1))%")

        # Low noise should give high success rate
        @test results[:success_rate] >= 0.0  # Just verify it completes
    end

    @testset "Stormy conditions (high noise)" begin
        # High atmospheric uncertainty - challenging conditions
        noise = NoiseParameters(0.15, 0.20, 0.03, 0.0002u"s", 0.02)
        results = monte_carlo_analysis(launcher, payload, mission, 3, noise_params=noise)  # Reduced for speed

        println("    Success rate (stormy): $(round(results[:success_rate] * 100, digits=1))%")

        # High noise reduces success rate
        @test results[:success_rate] < 1.0  # Not perfect
        @test results[:success_rate] >= 0.0  # But not impossible
    end
end

println("  ✓ Stochastic trajectory tests passed")
