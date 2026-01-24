"""
Tests for trajectory.jl

Tests:
- Payload configuration
- Mission profile setup
- Initial state creation
- State extraction
- Launch direction calculations
"""


@testset "Payload Configuration" begin
    payload = PayloadConfig(
        20.0u"kg",           # mass
        π * (0.15u"m")^2,   # area
        0.1u"m",             # nose_radius
        0.85,                # emissivity
        900.0u"J/(kg*K)",   # specific_heat
        288.15u"K",          # T_initial
        1500.0u"K"           # T_max
    )

    @test payload.mass == 20.0u"kg"
    @test payload.area ≈ π * (0.15u"m")^2
    @test payload.nose_radius == 0.1u"m"
    @test payload.emissivity == 0.85
    @test payload.specific_heat == 900.0u"J/(kg*K)"
    @test payload.T_initial == 288.15u"K"
    @test payload.T_max == 1500.0u"K"

    # Units check
    @test dimension(payload.mass) == dimension(1.0u"kg")
    @test dimension(payload.area) == dimension(1.0u"m^2")
    @test dimension(payload.nose_radius) == dimension(1.0u"m")
    @test dimension(payload.specific_heat) == dimension(1.0u"J/(kg*K)")
    @test dimension(payload.T_initial) == dimension(1.0u"K")
end

@testset "Mission Profile" begin
    mission = MissionProfile(
        0.0u"m",             # launch_altitude
        32.5u"°",            # launch_latitude (Caesarea)
        34.9u"°",            # launch_longitude
        90.0u"°",            # launch_azimuth (East)
        45.0u"°",            # launch_elevation
        8000.0u"m/s",       # target_velocity
        400.0u"km"           # target_altitude
    )

    @test mission.launch_altitude == 0.0u"m"
    @test mission.launch_latitude == 32.5u"°"
    @test mission.launch_longitude == 34.9u"°"
    @test mission.launch_azimuth == 90.0u"°"
    @test mission.launch_elevation == 45.0u"°"
    @test mission.target_velocity == 8000.0u"m/s"
    @test mission.target_altitude == 400.0u"km"

    # Units check
    @test dimension(mission.launch_altitude) == dimension(1.0u"m")
    @test dimension(mission.launch_latitude) == dimension(1.0u"°")
    @test dimension(mission.target_velocity) == dimension(1.0u"m/s")
end

@testset "Initial State Creation" begin
    # Create simple configurations
    launcher = create_uniform_launcher(
        length = 1000.0u"m",
        num_coils = 10,
        inductance_per_coil = 0.001u"H",
        resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 10.0u"F",
        voltage_per_coil = 5000.0u"V",
        gradient_per_coil = 0.002u"H/m"
    )

    payload = PayloadConfig(
        20.0u"kg", π * (0.15u"m")^2, 0.1u"m", 0.85,
        900.0u"J/(kg*K)", 288.15u"K", 1500.0u"K"
    )

    mission = MissionProfile(
        0.0u"m", 32.5u"°", 34.9u"°", 90.0u"°", 45.0u"°",
        8000.0u"m/s", 400.0u"km"
    )

    # Create initial state
    u0 = initial_state(launcher, payload, mission)

    # Check state vector length
    # [x, y, z, vx, vy, vz, T, I1...I10, Q1...Q10]
    # 3 + 3 + 1 + 10 + 10 = 27
    @test length(u0) == 27

    # Check initial velocity is zero
    @test u0[4] == 0.0u"m/s"
    @test u0[5] == 0.0u"m/s"
    @test u0[6] == 0.0u"m/s"

    # Check initial temperature
    @test u0[7] == 288.15u"K"

    # Check initial currents are zero
    for i in 8:17
        @test u0[i] == 0.0u"A"
    end

    # Check capacitors are charged
    for i in 18:27
        @test u0[i] > 0.0u"C"
    end
end

@testset "State Extraction" begin
    # Create mock state vector
    n_coils = 5
    u = vcat(
        [1000.0u"km", 0.0u"km", 0.0u"km"],          # position
        [1.0u"km/s", 2.0u"km/s", 0.5u"km/s"],      # velocity
        [500.0u"K"],                                 # temperature
        fill(1000.0u"A", n_coils),                   # currents
        fill(50000.0u"C", n_coils)                   # charges
    )

    state = extract_state(u, n_coils)

    # Check extracted components
    @test state.position[1] == 1000.0u"km"
    @test state.position[2] == 0.0u"km"
    @test state.position[3] == 0.0u"km"

    @test state.velocity[1] == 1.0u"km/s"
    @test state.velocity[2] == 2.0u"km/s"
    @test state.velocity[3] == 0.5u"km/s"

    @test state.temperature == 500.0u"K"

    @test length(state.currents) == n_coils
    @test state.currents[1] == 1000.0u"A"

    @test length(state.charges) == n_coils
    @test state.charges[1] == 50000.0u"C"
end

@testset "Launch Direction" begin
    # Launch position at Caesarea
    h = 0.0u"m"
    lat = 32.5u"°"
    lon = 34.9u"°"

    pos = position_from_altitude_angle(h, lat, lon)

    # Vertical launch (90° elevation)
    azimuth = 0.0  # radians (North)
    elevation = π/2  # radians (vertical)

    dir = launch_direction(azimuth, elevation, pos)

    # Should point radially outward
    @test norm(dir) ≈ 1.0 rtol=1e-6  # Unit vector
    @test dot(dir, pos/norm(pos)) > 0.9  # Mostly radial

    # Horizontal launch (0° elevation, East)
    azimuth_east = π/2
    elevation_horiz = 0.0

    dir_horiz = launch_direction(azimuth_east, elevation_horiz, pos)

    @test norm(dir_horiz) ≈ 1.0 rtol=1e-6  # Unit vector
    @test abs(dot(dir_horiz, pos/norm(pos))) < 0.1  # Perpendicular to radial

    # 45° elevation
    elevation_45 = π/4
    dir_45 = launch_direction(azimuth, elevation_45, pos)

    @test norm(dir_45) ≈ 1.0 rtol=1e-6  # Unit vector
end

@testset "Physical Reasonableness" begin
    # Payload mass should be positive
    payload = PayloadConfig(
        20.0u"kg", π * (0.15u"m")^2, 0.1u"m", 0.85,
        900.0u"J/(kg*K)", 288.15u"K", 1500.0u"K"
    )
    @test payload.mass > 0.0u"kg"
    @test payload.area > 0.0u"m^2"
    @test payload.nose_radius > 0.0u"m"
    @test 0.0 < payload.emissivity <= 1.0
    @test payload.T_initial > 0.0u"K"
    @test payload.T_max > payload.T_initial

    # Mission parameters should be in valid ranges
    mission = MissionProfile(
        0.0u"m", 32.5u"°", 34.9u"°", 90.0u"°", 45.0u"°",
        8000.0u"m/s", 400.0u"km"
    )
    @test mission.launch_altitude >= 0.0u"m"
    @test -90.0u"°" <= mission.launch_latitude <= 90.0u"°"
    @test 0.0u"°" <= mission.launch_azimuth <= 360.0u"°"
    @test 0.0u"°" <= mission.launch_elevation <= 90.0u"°"
    @test mission.target_velocity > 0.0u"m/s"
    @test mission.target_altitude > 0.0u"km"
end

@testset "Unit Consistency in State Vector" begin
    launcher = create_uniform_launcher(
        length = 1000.0u"m", num_coils = 10,
        inductance_per_coil = 0.001u"H", resistance_per_coil = 0.01u"Ω",
        capacitance_per_coil = 10.0u"F", voltage_per_coil = 5000.0u"V",
        gradient_per_coil = 0.002u"H/m"
    )

    payload = PayloadConfig(
        20.0u"kg", π * (0.15u"m")^2, 0.1u"m", 0.85,
        900.0u"J/(kg*K)", 288.15u"K", 1500.0u"K"
    )

    mission = MissionProfile(
        0.0u"m", 32.5u"°", 34.9u"°", 90.0u"°", 45.0u"°",
        8000.0u"m/s", 400.0u"km"
    )

    u0 = initial_state(launcher, payload, mission)

    # Position units
    @test dimension(u0[1]) == dimension(1.0u"m")
    @test dimension(u0[2]) == dimension(1.0u"m")
    @test dimension(u0[3]) == dimension(1.0u"m")

    # Velocity units
    @test dimension(u0[4]) == dimension(1.0u"m/s")
    @test dimension(u0[5]) == dimension(1.0u"m/s")
    @test dimension(u0[6]) == dimension(1.0u"m/s")

    # Temperature units
    @test dimension(u0[7]) == dimension(1.0u"K")

    # Current units
    for i in 8:17
        @test dimension(u0[i]) == dimension(1.0u"A")
    end

    # Charge units
    for i in 18:27
        @test dimension(u0[i]) == dimension(1.0u"C")
    end
end

println("  ✓ Trajectory system tests passed")
