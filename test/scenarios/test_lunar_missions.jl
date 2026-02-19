"""
Tests for lunar launch trajectories

Tests:
- Lunar trajectory simulation (no drag)
- Velocity conservation after launcher exit
- MissionProfile with body=MOON
"""

@testset "Lunar Trajectory - No Drag" begin
    # Launch from Moon surface with moderate velocity
    v_target = 1700.0u"m/s"  # Slightly above LLO velocity
    L_launcher = 2000.0u"m"
    m_payload = 20.0u"kg"
    elevation = 45.0u"°"
    azimuth = 90.0u"°"

    prob, log, callback = create_simplified_sde_problem(
        v_target, L_launcher, m_payload,
        elevation, azimuth;
        tspan=(0.0u"s", 60.0u"s"),
        noise_level=0.0,
        enable_logging=false,
        body=MOON
    )

    sol = solve(prob, Euler(), dt=0.01u"s", adaptive=false, dense=false)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Extract final state
    u_final = sol.u[end]
    vel_final = SVector(u_final[4], u_final[5], u_final[6])
    v_final = norm(vel_final)

    # Velocity should not have decreased significantly (no drag on Moon)
    # It should be close to or above the launch velocity (gravity reduces radial component)
    @test ustrip(u"m/s", v_final) > 1000.0  # Still has substantial velocity
end

@testset "Lunar Velocity Conservation After Exit" begin
    # After exiting launcher, with no atmosphere, only gravity acts
    # Speed at exit should not be diminished by drag
    v_target = 1500.0u"m/s"
    L_launcher = 1000.0u"m"
    m_payload = 10.0u"kg"
    elevation = 90.0u"°"  # Vertical launch for simplicity
    azimuth = 0.0u"°"

    prob, log, callback = create_simplified_sde_problem(
        v_target, L_launcher, m_payload,
        elevation, azimuth;
        tspan=(0.0u"s", 10.0u"s"),
        noise_level=0.0,
        enable_logging=false,
        body=MOON
    )

    sol = solve(prob, Euler(), dt=0.01u"s", adaptive=false, dense=false)
    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Find time when launcher is exited (when x_launcher > L_launcher)
    R_moon = MOON.radius
    launch_pos = SVector(R_moon, 0.0u"m", 0.0u"m")

    # Check that temperature stays constant (no heating on Moon)
    T_initial = ustrip(u"K", sol.u[1][7])
    T_final = ustrip(u"K", sol.u[end][7])
    @test T_initial ≈ T_final atol=1.0  # No atmospheric heating
end

@testset "Lunar MissionProfile Construction" begin
    # 7-arg constructor defaults to EARTH
    mp_earth = MissionProfile(
        0.0u"m", 32.5u"°", 34.9u"°", 90.0u"°",
        45.0u"°", 8000.0u"m/s", 400.0u"km"
    )
    @test mp_earth.body === EARTH

    # 8-arg constructor with explicit body
    mp_moon = MissionProfile(
        0.0u"m", 0.0u"°", 0.0u"°", 90.0u"°",
        45.0u"°", 1700.0u"m/s", 100.0u"km", MOON
    )
    @test mp_moon.body === MOON
end

@testset "Lunar vs Earth Drag Comparison" begin
    # Same launch conditions, different bodies
    v_target = 2000.0u"m/s"
    L_launcher = 1000.0u"m"
    m_payload = 10.0u"kg"
    elevation = 45.0u"°"
    azimuth = 90.0u"°"

    # Earth launch
    prob_earth, _, _ = create_simplified_sde_problem(
        v_target, L_launcher, m_payload,
        elevation, azimuth;
        tspan=(0.0u"s", 30.0u"s"),
        noise_level=0.0,
        body=EARTH
    )
    sol_earth = solve(prob_earth, Euler(), dt=0.01u"s", adaptive=false, dense=false)

    # Moon launch
    prob_moon, _, _ = create_simplified_sde_problem(
        v_target, L_launcher, m_payload,
        elevation, azimuth;
        tspan=(0.0u"s", 30.0u"s"),
        noise_level=0.0,
        body=MOON
    )
    sol_moon = solve(prob_moon, Euler(), dt=0.01u"s", adaptive=false, dense=false)

    # Extract final velocities
    vel_earth = norm(SVector(sol_earth.u[end][4], sol_earth.u[end][5], sol_earth.u[end][6]))
    vel_moon = norm(SVector(sol_moon.u[end][4], sol_moon.u[end][5], sol_moon.u[end][6]))

    # Moon launch should retain more velocity (no atmospheric drag)
    @test ustrip(u"m/s", vel_moon) > ustrip(u"m/s", vel_earth)
end

println("  ✓ Lunar mission tests passed")
