"""
# Simulation Module

High-level simulation interface for mission scenarios.

Provides convenient functions for common mission types:
- LEO cargo delivery
- Hypersonic testing
- Suborbital missions
- Custom trajectories
"""

using Unitful
using DifferentialEquations

"""
    create_default_payload(;mass=20.0u"kg")

Create payload configuration with standard assumptions.

# Arguments
- `mass`: Payload mass [kg]

# Returns
- PayloadConfig with default parameters
"""
function create_default_payload(;mass=20.0u"kg")
    return PayloadConfig(
        mass,                           # mass
        π * (0.15u"m")^2,              # area (0.3m diameter)
        0.1u"m",                        # nose_radius
        0.85,                           # emissivity
        900.0u"J/(kg*K)",              # specific_heat (aluminum-like)
        288.15u"K",                     # T_initial (room temp)
        1500.0u"K"                      # T_max
    )
end

"""
    mission_leo_cargo(;
        payload_mass=20.0u"kg",
        target_altitude=400.0u"km",
        launcher_length=2000.0u"m",
        n_coils=200
    )

Simulate LEO cargo delivery mission.

# Arguments
- `payload_mass`: Mass of payload [kg]
- `target_altitude`: Target orbital altitude [m]
- `launcher_length`: Length of launcher tube [m]
- `n_coils`: Number of acceleration coils

# Returns
- Solution object from trajectory simulation
"""
function mission_leo_cargo(;
    payload_mass=20.0u"kg",
    target_altitude=400.0u"km",
    launcher_length=2000.0u"m",
    n_coils=200
)
    # Create launcher configuration
    launcher = create_uniform_launcher(
        length = launcher_length,
        num_coils = n_coils,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 10.0u"F",
        voltage_per_coil = 5000.0u"V",
        gradient_per_coil = 0.002u"H/m"
    )

    # Create payload
    payload = create_default_payload(mass=payload_mass)

    # Create mission profile (Caesarea, Israel location)
    mission = MissionProfile(
        0.0u"m",                        # launch_altitude (sea level)
        32.5u"°",                       # launch_latitude (Caesarea)
        34.9u"°",                       # launch_longitude
        90.0u"°",                       # launch_azimuth (East for max velocity boost)
        45.0u"°",                       # launch_elevation (45° standard)
        8000.0u"m/s",                  # target_velocity (orbital)
        target_altitude                 # target_altitude
    )

    # Simulate
    sol = simulate_trajectory(
        launcher,
        payload,
        mission,
        tspan = (0.0u"s", 300.0u"s"),  # 5 minutes
        saveat = 0.1u"s"                # Save every 0.1 seconds
    )

    return sol
end

"""
    mission_hypersonic_test(;
        payload_mass=10.0u"kg",
        target_mach=6.0,
        launcher_length=1000.0u"m",
        n_coils=100
    )

Simulate hypersonic test flight.

# Arguments
- `payload_mass`: Mass of test article [kg]
- `target_mach`: Target Mach number
- `launcher_length`: Length of launcher tube [m]
- `n_coils`: Number of acceleration coils

# Returns
- Solution object from trajectory simulation
"""
function mission_hypersonic_test(;
    payload_mass=10.0u"kg",
    target_mach=6.0,
    launcher_length=1000.0u"m",
    n_coils=100
)
    # Target velocity (Mach 6 at 25km altitude)
    # Speed of sound at 25km ≈ 300 m/s
    target_velocity = target_mach * 300.0u"m/s"

    # Create launcher (smaller for hypersonic testing)
    launcher = create_uniform_launcher(
        length = launcher_length,
        num_coils = n_coils,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 5.0u"F",
        voltage_per_coil = 3000.0u"V",
        gradient_per_coil = 0.002u"H/m"
    )

    # Create payload
    payload = create_default_payload(mass=payload_mass)

    # Create mission profile (shallower angle for sustained flight)
    mission = MissionProfile(
        0.0u"m",                        # launch_altitude
        32.5u"°",                       # launch_latitude (Caesarea)
        34.9u"°",                       # launch_longitude
        90.0u"°",                       # launch_azimuth
        25.0u"°",                       # launch_elevation (shallower)
        target_velocity,                # target_velocity
        25.0u"km"                       # target_altitude (hypersonic corridor)
    )

    # Simulate
    sol = simulate_trajectory(
        launcher,
        payload,
        mission,
        tspan = (0.0u"s", 120.0u"s"),  # 2 minutes
        saveat = 0.1u"s"
    )

    return sol
end

"""
    mission_suborbital(;
        payload_mass=50.0u"kg",
        range=1000.0u"km",
        launcher_length=2000.0u"m",
        n_coils=200
    )

Simulate suborbital point-to-point delivery.

# Arguments
- `payload_mass`: Mass of payload [kg]
- `range`: Ground range to target [km]
- `launcher_length`: Length of launcher tube [m]
- `n_coils`: Number of acceleration coils

# Returns
- Solution object from trajectory simulation
"""
function mission_suborbital(;
    payload_mass=50.0u"kg",
    range=1000.0u"km",
    launcher_length=2000.0u"m",
    n_coils=200
)
    # Calculate required velocity and angle for range
    # (Simplified ballistic trajectory)
    # v² sin(2θ) / g = range
    # Assume 45° for maximum range: v = sqrt(range * g)

    g = 9.81u"m/s^2"
    target_velocity = sqrt(range * g)

    # Create launcher
    launcher = create_uniform_launcher(
        length = launcher_length,
        num_coils = n_coils,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 10.0u"F",
        voltage_per_coil = 4000.0u"V",
        gradient_per_coil = 0.002u"H/m"
    )

    # Create payload
    payload = create_default_payload(mass=payload_mass)

    # Create mission profile
    mission = MissionProfile(
        0.0u"m",                        # launch_altitude
        32.5u"°",                       # launch_latitude
        34.9u"°",                       # launch_longitude
        90.0u"°",                       # launch_azimuth
        45.0u"°",                       # launch_elevation (optimal for range)
        target_velocity,                # target_velocity
        100.0u"km"                      # max altitude (approximate)
    )

    # Simulate
    sol = simulate_trajectory(
        launcher,
        payload,
        mission,
        tspan = (0.0u"s", 600.0u"s"),  # 10 minutes
        saveat = 1.0u"s"
    )

    return sol
end

"""
    analyze_performance(sol, n_coils)

Analyze simulation results and extract performance metrics.

# Arguments
- `sol`: ODE solution from trajectory simulation
- `n_coils`: Number of coils in launcher

# Returns
- Dict with performance metrics
"""
function analyze_performance(sol, n_coils)
    # Extract trajectory data
    traj = extract_trajectory_data(sol, n_coils)

    # Calculate metrics
    metrics = Dict(
        "final_velocity" => traj.speeds[end],
        "final_altitude" => traj.altitudes[end],
        "max_mach" => maximum(traj.mach_numbers),
        "max_temperature" => maximum(traj.temperatures),
        "max_gload" => 0.0,  # Would need to calculate from accelerations
        "flight_time" => traj.times[end],
        "apogee" => maximum(traj.altitudes)
    )

    return metrics
end

"""
    print_mission_summary(sol, n_coils, mission_name="Mission")

Print human-readable mission summary.

# Arguments
- `sol`: Solution from trajectory simulation
- `n_coils`: Number of coils
- `mission_name`: Name of mission for display
"""
function print_mission_summary(sol, n_coils, mission_name="Mission")
    metrics = analyze_performance(sol, n_coils)
    traj = extract_trajectory_data(sol, n_coils)

    println("=" ^ 70)
    println("$mission_name SUMMARY")
    println("=" ^ 70)
    println()
    println("FINAL STATE:")
    println("  Velocity:              $(metrics["final_velocity"])")
    println("  Altitude:              $(metrics["final_altitude"])")
    println()
    println("MISSION PERFORMANCE:")
    println("  Maximum Mach number:   $(round(metrics["max_mach"], digits=2))")
    println("  Maximum temperature:   $(metrics["max_temperature"])")
    println("  Apogee altitude:       $(metrics["apogee"])")
    println("  Flight time:           $(metrics["flight_time"])")
    println()
    println("INITIAL CONDITIONS:")
    println("  Launch velocity:       $(traj.speeds[1])")
    println("  Launch altitude:       $(traj.altitudes[1])")
    println()
    println("=" ^ 70)
end

export create_default_payload
export mission_leo_cargo, mission_hypersonic_test, mission_suborbital
export analyze_performance, print_mission_summary
