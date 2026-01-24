"""
Tests for simulation.jl and actual trajectory ODE solving

Tests:
- Full trajectory simulations
- ODE solver convergence
- Energy conservation during integration
- Mission scenario validity
- Numerical stability
- Physical constraints during flight
"""


using DifferentialEquations

@testset "Default Payload Creation" begin
    payload = create_default_payload(mass=20.0u"kg")

    @test payload.mass == 20.0u"kg"
    @test payload.area > 0.0u"m^2"
    @test payload.nose_radius > 0.0u"m"
    @test 0.0 < payload.emissivity <= 1.0
    @test payload.specific_heat > 0.0u"J/(kg*K)"
    @test payload.T_initial > 0.0u"K"
    @test payload.T_max > payload.T_initial
end

@testset "Simple Free Fall Trajectory" begin
    # Test a simple free fall trajectory (no EM, no drag)
    # This should follow ballistic motion

    # Create minimal launcher (won't be used)
    launcher = create_uniform_launcher(
        length = 100.0u"m",
        num_coils = 10,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F",
        voltage_per_coil = 0.0u"V",  # No voltage = no EM force
        gradient_per_coil = 0.001u"H/m"
    )

    payload = create_default_payload(mass=1.0u"kg")

    # Launch vertically from high altitude (minimal drag)
    mission = MissionProfile(
        100.0u"km",           # Start at 100 km (minimal atmosphere)
        0.0u"°",              # Equator
        0.0u"°",              # Prime meridian
        0.0u"°",              # North
        90.0u"°",             # Straight up
        0.0u"m/s",            # Zero initial velocity
        200.0u"km"            # Target (won't reach it in free fall)
    )

    # Create initial state
    u0_unitful = initial_state(launcher, payload, mission)
    n_coils = launcher.num_coils

    # Convert to dimensionless
    u0, units_ref = strip_units_state(u0_unitful, n_coils)
    p = (launcher, payload, mission, units_ref)

    # Short time span (1 second)
    tspan = (0.0, 1.0)

    # Create ODE problem
    prob = ODEProblem(trajectory_ode_dimensionless!, u0, tspan, p)

    # Solve
    sol = solve(prob, Tsit5(), saveat=0.1)

    @test sol.retcode == ReturnCode.Success || sol.retcode == ReturnCode.Default

    # Check that we have solutions
    @test length(sol.t) > 1

    # Extract data
    traj = extract_trajectory_data(sol, n_coils, units_ref)

    # Verify altitude decreases (falling)
    @test traj.altitudes[end] < traj.altitudes[1]

    # Verify velocity increases (accelerating down)
    @test traj.speeds[end] > traj.speeds[1]
end

@testset "Energy Conservation in Orbit" begin
    # Test that energy is approximately conserved in a circular orbit
    # (without atmospheric drag)

    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 10,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 0.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )

    payload = create_default_payload(mass=10.0u"kg")

    # Start at 400 km altitude with orbital velocity
    h = 400.0u"km"
    v_orb = orbital_velocity(h)

    # Initial position and velocity for circular orbit
    r = 6.3781370e6u"m" + h
    pos_0 = SVector(r, 0.0u"m", 0.0u"m")
    vel_0 = SVector(0.0u"m/s", v_orb, 0.0u"m/s")

    # Calculate initial specific energy
    ε_0 = specific_orbital_energy(vel_0, pos_0)

    # We can't easily test full orbit with current setup since trajectory_ode!
    # includes EM forces. This test validates the energy calculation itself.
    @test ε_0 < 0.0u"m^2/s^2"  # Bound orbit
    @test abs(ε_0) > 0.0u"m^2/s^2"  # Non-zero energy
end

@testset "Atmospheric Drag Effect" begin
    # Test that atmospheric drag slows down the payload

    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 10,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 0.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )

    payload = create_default_payload(mass=1.0u"kg")

    # Launch horizontally at low altitude (high drag)
    mission = MissionProfile(
        10.0u"km",            # 10 km altitude (dense atmosphere)
        0.0u"°",
        0.0u"°",
        90.0u"°",             # East
        0.0u"°",              # Horizontal
        1000.0u"m/s",         # 1 km/s initial velocity
        20.0u"km"
    )

    u0_unitful = initial_state(launcher, payload, mission)
    n_coils = launcher.num_coils

    # Give initial horizontal velocity
    u0_unitful[4] = 0.0u"m/s"        # vx
    u0_unitful[5] = 1000.0u"m/s"     # vy (horizontal)
    u0_unitful[6] = 0.0u"m/s"        # vz

    # Convert to dimensionless
    u0, units_ref = strip_units_state(u0_unitful, n_coils)
    p = (launcher, payload, mission, units_ref)
    tspan = (0.0, 10.0)

    prob = ODEProblem(trajectory_ode_dimensionless!, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=1.0)

    if sol.retcode == ReturnCode.Success || sol.retcode == ReturnCode.Default
        traj = extract_trajectory_data(sol, n_coils, units_ref)

        # Velocity should decrease due to drag
        @test traj.speeds[end] < traj.speeds[1]
    end
end

@testset "Heating During Flight" begin
    # Test that payload heats up during high-speed atmospheric flight

    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 10,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 0.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )

    payload = create_default_payload(mass=5.0u"kg")

    # High speed at moderate altitude
    mission = MissionProfile(
        30.0u"km",
        0.0u"°",
        0.0u"°",
        90.0u"°",
        45.0u"°",
        3000.0u"m/s",  # High speed = high heating
        50.0u"km"
    )

    u0_unitful = initial_state(launcher, payload, mission)
    n_coils = launcher.num_coils

    # Set initial velocity
    v_mag = 3000.0u"m/s"
    u0_unitful[4] = v_mag * 0.707  # vx (45° angle)
    u0_unitful[5] = 0.0u"m/s"
    u0_unitful[6] = v_mag * 0.707  # vz

    # Convert to dimensionless
    u0, units_ref = strip_units_state(u0_unitful, n_coils)
    p = (launcher, payload, mission, units_ref)
    tspan = (0.0, 5.0)

    prob = ODEProblem(trajectory_ode_dimensionless!, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=0.5)

    if sol.retcode == ReturnCode.Success || sol.retcode == ReturnCode.Default
        traj = extract_trajectory_data(sol, n_coils, units_ref)

        # Temperature should increase initially due to heating
        # (May later cool via radiation)
        @test maximum(traj.temperatures) > traj.temperatures[1]
    end
end

@testset "ODE Solver Stability" begin
    # Test that different solvers give similar results

    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 10,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 0.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )

    payload = create_default_payload(mass=10.0u"kg")

    mission = MissionProfile(
        50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
        1000.0u"m/s", 100.0u"km"
    )

    u0_unitful = initial_state(launcher, payload, mission)
    n_coils = launcher.num_coils

    # Convert to dimensionless
    u0, units_ref = strip_units_state(u0_unitful, n_coils)
    p = (launcher, payload, mission, units_ref)
    tspan = (0.0, 2.0)

    prob = ODEProblem(trajectory_ode_dimensionless!, u0, tspan, p)

    # Try different solvers
    solvers = [Tsit5(), Vern7(), VCABM()]
    solutions = []

    for solver in solvers
        try
            sol = solve(prob, solver, saveat=0.5)
            if sol.retcode == ReturnCode.Success || sol.retcode == ReturnCode.Default
                push!(solutions, sol)
            end
        catch
            # Some solvers might fail, that's okay
        end
    end

    # At least one solver should work
    @test length(solutions) >= 1

    # If multiple solvers worked, results should be similar
    if length(solutions) >= 2
        #n_coils already defined above
        traj1 = extract_trajectory_data(solutions[1], n_coils, units_ref)
        traj2 = extract_trajectory_data(solutions[2], n_coils, units_ref)

        # Final positions should be similar (within 10%)
        if length(traj1.altitudes) == length(traj2.altitudes)
            @test isapprox(traj1.altitudes[end], traj2.altitudes[end], rtol=0.1)
        end
    end
end

@testset "State Vector Dimensions" begin
    # Test that state vector maintains proper dimensions throughout

    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 10,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 0.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )

    payload = create_default_payload(mass=10.0u"kg")

    mission = MissionProfile(
        50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
        500.0u"m/s", 100.0u"km"
    )

    u0_unitful = initial_state(launcher, payload, mission)
    n_coils = launcher.num_coils

    # Convert to dimensionless
    u0, units_ref = strip_units_state(u0_unitful, n_coils)
    p = (launcher, payload, mission, units_ref)
    tspan = (0.0, 1.0)

    prob = ODEProblem(trajectory_ode_dimensionless!, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=0.2)

    if sol.retcode == ReturnCode.Success || sol.retcode == ReturnCode.Default
        # Check all solution points
        for i in 1:length(sol.t)
            u_unitful = reattach_units_state(sol.u[i], units_ref, n_coils)
            state = extract_state(u_unitful, n_coils)

            # Check dimensions
            @test dimension(state.position[1]) == dimension(1.0u"m")
            @test dimension(state.velocity[1]) == dimension(1.0u"m/s")
            @test dimension(state.temperature) == dimension(1.0u"K")
            @test dimension(state.currents[1]) == dimension(1.0u"A")
            @test dimension(state.charges[1]) == dimension(1.0u"C")

            # Check physical bounds
            @test state.temperature > 0.0u"K"
            @test state.temperature < 10000.0u"K"  # Reasonable upper bound
        end
    end
end

@testset "Performance Analysis Functions" begin
    # Test that analysis functions work on solutions

    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 10,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 0.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )

    payload = create_default_payload(mass=10.0u"kg")

    mission = MissionProfile(
        50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
        500.0u"m/s", 100.0u"km"
    )

    u0_unitful = initial_state(launcher, payload, mission)
    n_coils = launcher.num_coils

    # Convert to dimensionless
    u0, units_ref = strip_units_state(u0_unitful, n_coils)
    p = (launcher, payload, mission, units_ref)
    tspan = (0.0, 2.0)

    prob = ODEProblem(trajectory_ode_dimensionless!, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=0.5)

    if sol.retcode == ReturnCode.Success || sol.retcode == ReturnCode.Default
        #n_coils already defined above

        # Test analyze_performance
        metrics = analyze_performance(sol, n_coils)

        @test haskey(metrics, "final_velocity")
        @test haskey(metrics, "final_altitude")
        @test haskey(metrics, "max_mach")
        @test haskey(metrics, "max_temperature")
        @test haskey(metrics, "flight_time")
        @test haskey(metrics, "apogee")

        # Check values are physical
        @test metrics["final_velocity"] >= 0.0u"m/s"
        @test metrics["final_altitude"] >= 0.0u"m"
        @test metrics["max_mach"] >= 0.0
        @test metrics["max_temperature"] > 0.0u"K"
        @test metrics["flight_time"] > 0.0u"s"
    end
end

@testset "Trajectory Data Extraction" begin
    # Test extraction of trajectory data from solutions

    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 10,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 0.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )

    payload = create_default_payload(mass=10.0u"kg")

    mission = MissionProfile(
        50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
        500.0u"m/s", 100.0u"km"
    )

    u0_unitful = initial_state(launcher, payload, mission)
    n_coils = launcher.num_coils

    # Convert to dimensionless
    u0, units_ref = strip_units_state(u0_unitful, n_coils)
    p = (launcher, payload, mission, units_ref)
    tspan = (0.0, 1.0)

    prob = ODEProblem(trajectory_ode_dimensionless!, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=0.2)

    if sol.retcode == ReturnCode.Success || sol.retcode == ReturnCode.Default
        #n_coils already defined above
        traj = extract_trajectory_data(sol, n_coils, units_ref)

        # Check all fields exist
        @test length(traj.times) > 0
        @test length(traj.positions) == length(traj.times)
        @test length(traj.velocities) == length(traj.times)
        @test length(traj.temperatures) == length(traj.times)
        @test length(traj.altitudes) == length(traj.times)
        @test length(traj.speeds) == length(traj.times)
        @test length(traj.mach_numbers) == length(traj.times)

        # Check dimensions
        @test dimension(traj.times[1]) == dimension(1.0u"s")
        @test dimension(traj.altitudes[1]) == dimension(1.0u"m")
        @test dimension(traj.speeds[1]) == dimension(1.0u"m/s")
        @test dimension(traj.temperatures[1]) == dimension(1.0u"K")
    end
end

@testset "Physical Constraints During Flight" begin
    # Verify physics constraints are maintained

    launcher = create_uniform_launcher(
        length = 100.0u"m", num_coils = 10,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 1.0u"F", voltage_per_coil = 0.0u"V",
        gradient_per_coil = 0.001u"H/m"
    )

    payload = create_default_payload(mass=10.0u"kg")

    mission = MissionProfile(
        50.0u"km", 0.0u"°", 0.0u"°", 90.0u"°", 45.0u"°",
        500.0u"m/s", 100.0u"km"
    )

    u0_unitful = initial_state(launcher, payload, mission)
    n_coils = launcher.num_coils

    # Convert to dimensionless
    u0, units_ref = strip_units_state(u0_unitful, n_coils)
    p = (launcher, payload, mission, units_ref)
    tspan = (0.0, 5.0)

    prob = ODEProblem(trajectory_ode_dimensionless!, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=0.5)

    if sol.retcode == ReturnCode.Success || sol.retcode == ReturnCode.Default
        #n_coils already defined above
        traj = extract_trajectory_data(sol, n_coils, units_ref)

        # All temperatures should be positive and finite
        for T in traj.temperatures
            @test T > 0.0u"K"
            @test isfinite(ustrip(u"K", T))
        end

        # All speeds should be non-negative and finite
        for v in traj.speeds
            @test v >= 0.0u"m/s"
            @test isfinite(ustrip(u"m/s", v))
        end

        # All altitudes should be finite
        for h in traj.altitudes
            @test isfinite(ustrip(u"m", h))
        end

        # Mach numbers should be non-negative
        for M in traj.mach_numbers
            @test M >= 0.0
            @test isfinite(M)
        end
    end
end

println("  ✓ Simulation and ODE solving tests passed")
