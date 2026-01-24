"""
# Trajectory Module

Complete system of ODEs for electromagnetic launch trajectory simulation.

Combines all physics:
- Electromagnetic acceleration
- Atmospheric drag
- Aerodynamic heating
- Gravitational forces (including J2)

## State Vector

The complete state is:
```
u = [x, y, z, vₓ, vᵧ, vᵤ, T, I₁, I₂, ..., Iₙ, Q₁, Q₂, ..., Qₙ]
```

where:
- `[x, y, z]` = position [m]
- `[vₓ, vᵧ, vᵤ]` = velocity [m/s]
- `T` = payload temperature [K]
- `Iᵢ` = current in coil i [A]
- `Qᵢ` = charge in capacitor i [C]

## Equations of Motion

Position derivatives:
```math
\\frac{d\\vec{r}}{dt} = \\vec{v}
```

Velocity derivatives (Newton's 2nd law):
```math
\\frac{d\\vec{v}}{dt} = \\frac{\\vec{F}_{em} + \\vec{F}_{drag} + \\vec{F}_{grav}}{m}
```

Temperature derivative (energy balance):
```math
\\frac{dT}{dt} = \\frac{\\dot{q} A - \\epsilon \\sigma A (T^4 - T_{amb}^4)}{m c_p}
```

Electrical circuit equations (for each coil):
```math
\\begin{aligned}
\\frac{dI_i}{dt} &= \\frac{V_i(t) - Q_i/C_i - R_i I_i}{L_i} \\\\
\\frac{dQ_i}{dt} &= I_i
\\end{aligned}
```
"""

using Unitful
using LinearAlgebra
using StaticArrays
using DifferentialEquations

"""
    PayloadConfig

Configuration for the payload being launched.

# Fields
- `mass`: Payload mass [kg]
- `area`: Reference area for drag [m²]
- `nose_radius`: Nose radius for heating [m]
- `emissivity`: Thermal emissivity (dimensionless)
- `specific_heat`: Specific heat capacity [J/(kg·K)]
- `T_initial`: Initial temperature [K]
- `T_max`: Maximum allowable temperature [K]
"""
struct PayloadConfig
    mass
    area
    nose_radius
    emissivity::Float64  # dimensionless
    specific_heat
    T_initial
    T_max
end

"""
    MissionProfile

Defines the mission parameters and launch conditions.

# Fields
- `launch_altitude`: Launch site altitude [m]
- `launch_latitude`: Launch site latitude [degrees]
- `launch_longitude`: Launch site longitude [degrees]
- `launch_azimuth`: Launch azimuth angle [degrees] (0=North, 90=East)
- `launch_elevation`: Launch elevation angle [degrees] (0=horizontal, 90=vertical)
- `target_velocity`: Target final velocity [m/s]
- `target_altitude`: Target final altitude (if orbital) [m]
"""
struct MissionProfile
    launch_altitude
    launch_latitude
    launch_longitude
    launch_azimuth
    launch_elevation
    target_velocity
    target_altitude
end

"""
    initial_state(launcher_config, payload_config, mission_profile)

Create initial state vector for ODE integration.

# Returns
State vector: [x, y, z, vₓ, vᵧ, vᵤ, T, I₁...Iₙ, Q₁...Qₙ]
"""
function initial_state(
    launcher_config::LauncherConfig,
    payload_config::PayloadConfig,
    mission_profile::MissionProfile
)
    # Initial position (at launch site, start of tube)
    pos_0 = position_from_altitude_angle(
        mission_profile.launch_altitude,
        mission_profile.launch_latitude,
        mission_profile.launch_longitude
    )

    # Initial velocity (zero - starting from rest)
    vel_0 = SVector(0.0u"m/s", 0.0u"m/s", 0.0u"m/s")

    # Initial temperature
    T_0 = payload_config.T_initial

    # Initial currents (all zero)
    n_coils = launcher_config.num_coils
    I_0 = zeros(n_coils) * 1.0u"A"

    # Initial charges (capacitors fully charged)
    Q_0 = [coil.capacitance * coil.voltage for coil in launcher_config.coils]

    # Combine into state vector
    # Note: In Julia with Unitful, we need to carefully handle arrays
    state = vcat(
        collect(pos_0),      # [x, y, z]
        collect(vel_0),      # [vₓ, vᵧ, vᵤ]
        [T_0],               # [T]
        collect(I_0),        # [I₁...Iₙ]
        collect(Q_0)         # [Q₁...Qₙ]
    )

    return state
end

"""
    extract_state(u, n_coils)

Extract components from state vector.

# Returns
Named tuple: (position, velocity, temperature, currents, charges)
"""
function extract_state(u, n_coils)
    position = SVector(u[1], u[2], u[3])
    velocity = SVector(u[4], u[5], u[6])
    temperature = u[7]
    currents = u[8:(7+n_coils)]
    charges = u[(8+n_coils):(7+2*n_coils)]

    return (
        position = position,
        velocity = velocity,
        temperature = temperature,
        currents = currents,
        charges = charges
    )
end

"""
    strip_units_state(u_unitful, n_coils)

Convert unitful state vector to dimensionless for ODE integration.

Returns: (u_dimensionless, units_tuple)
"""
function strip_units_state(u_unitful, n_coils)
    # Define reference units for normalization
    units = (
        L = 1.0u"m",           # Length
        V = 1.0u"m/s",         # Velocity
        T = 1.0u"K",           # Temperature
        I = 1.0u"A",           # Current
        Q = 1.0u"C"            # Charge
    )

    # Strip units from state vector
    u_dimensionless = zeros(Float64, length(u_unitful))

    # Position [1:3]
    for i in 1:3
        u_dimensionless[i] = ustrip(u"m", u_unitful[i])
    end

    # Velocity [4:6]
    for i in 4:6
        u_dimensionless[i] = ustrip(u"m/s", u_unitful[i])
    end

    # Temperature [7]
    u_dimensionless[7] = ustrip(u"K", u_unitful[7])

    # Currents [8:(7+n_coils)]
    for i in 1:n_coils
        u_dimensionless[7+i] = ustrip(u"A", u_unitful[7+i])
    end

    # Charges [(8+n_coils):(7+2*n_coils)]
    for i in 1:n_coils
        u_dimensionless[7+n_coils+i] = ustrip(u"C", u_unitful[7+n_coils+i])
    end

    return u_dimensionless, units
end

"""
    reattach_units_state(u_dimensionless, units, n_coils)

Convert dimensionless state vector back to unitful.
"""
function reattach_units_state(u_dimensionless, units, n_coils)
    u_unitful = similar(u_dimensionless, Any)

    # Position [1:3]
    for i in 1:3
        u_unitful[i] = u_dimensionless[i] * units.L
    end

    # Velocity [4:6]
    for i in 4:6
        u_unitful[i] = u_dimensionless[i] * units.V
    end

    # Temperature [7]
    u_unitful[7] = u_dimensionless[7] * units.T

    # Currents [8:(7+n_coils)]
    for i in 1:n_coils
        u_unitful[7+i] = u_dimensionless[7+i] * units.I
    end

    # Charges [(8+n_coils):(7+2*n_coils)]
    for i in 1:n_coils
        u_unitful[7+n_coils+i] = u_dimensionless[7+n_coils+i] * units.Q
    end

    return u_unitful
end

"""
    launch_direction(azimuth, elevation, position)

Calculate launch direction vector in ECI frame.

# Arguments
- `azimuth`: Azimuth angle [rad] (0=North, π/2=East)
- `elevation`: Elevation angle [rad] (0=horizontal, π/2=up)
- `position`: Launch position vector [m] (for local frame)

# Returns
- Unit direction vector in ECI frame
"""
function launch_direction(azimuth, elevation, position)
    # Simplified: assume launch direction in local tangent plane
    # More sophisticated: transform from local NEU to ECI

    # Local up direction (radial from Earth center)
    r_hat = position / norm(position)

    # Local north (simplified - in xy plane perpendicular to r_hat)
    z_axis = SVector(0.0, 0.0, 1.0)
    east = cross(z_axis, r_hat)
    east = east / norm(east)
    north = cross(r_hat, east)

    # Launch direction in local frame
    d_local = cos(elevation) * (cos(azimuth) * north + sin(azimuth) * east) +
              sin(elevation) * r_hat

    return d_local / norm(d_local)
end

"""
    trajectory_ode_dimensionless!(du, u, p, t)

Dimensionless ODE system for electromagnetic launch trajectory.
Works with Float64 arrays (no units) for compatibility with DifferentialEquations.jl.

# Arguments
- `du`: Derivative vector (output, dimensionless)
- `u`: State vector [x,y,z,vₓ,vᵧ,vᵤ,T,I₁...Iₙ,Q₁...Qₙ] (dimensionless)
- `p`: Parameters (launcher_config, payload_config, mission_profile, units_tuple)
- `t`: Time [s] (dimensionless)
"""
function trajectory_ode_dimensionless!(du, u, p, t)
    launcher_config, payload_config, mission_profile, units_ref = p
    n_coils = launcher_config.num_coils

    # Reattach units to state for physics calculations
    pos = SVector(u[1], u[2], u[3]) * units_ref.L
    vel = SVector(u[4], u[5], u[6]) * units_ref.V
    T = u[7] * units_ref.T
    currents = u[8:(7+n_coils)] * units_ref.I
    charges = u[(8+n_coils):(7+2*n_coils)] * units_ref.Q

    # Calculate altitude
    h = altitude_from_position(pos)

    # =================================================================
    # ELECTROMAGNETIC FORCE
    # =================================================================
    # Project position onto launcher axis to get "x" coordinate
    # (Simplified: use distance from launch site)
    launch_pos = position_from_altitude_angle(
        mission_profile.launch_altitude,
        mission_profile.launch_latitude,
        mission_profile.launch_longitude
    )
    x_launcher = norm(pos - launch_pos)  # Position along launcher

    # Only apply EM force if inside launcher tube
    if x_launcher < launcher_config.length
        F_em_mag = total_em_force(x_launcher, currents, launcher_config)

        # Direction along launcher axis (simplified: radially outward)
        launch_dir = launch_direction(
            ustrip(u"rad", mission_profile.launch_azimuth),
            ustrip(u"rad", mission_profile.launch_elevation),
            pos
        )
        F_em = F_em_mag * launch_dir
    else
        # Outside launcher - no EM force
        F_em = SVector(0.0u"N", 0.0u"N", 0.0u"N")
    end

    # =================================================================
    # DRAG FORCE (with vacuum tube support)
    # =================================================================
    if h < 100.0u"km"  # Only calculate drag in atmosphere
        # Apply vacuum reduction if inside launcher tube
        if x_launcher < launcher_config.length
            # Inside tube - reduce drag by vacuum_pressure_ratio
            F_drag_vec = drag_force(vel, h, payload_config.area) * launcher_config.vacuum_pressure_ratio
        else
            # Outside tube - normal atmospheric drag
            F_drag_vec = drag_force(vel, h, payload_config.area)
        end
    else
        F_drag_vec = SVector(0.0u"N", 0.0u"N", 0.0u"N")
    end

    # =================================================================
    # GRAVITATIONAL FORCE
    # =================================================================
    a_grav = gravity_with_J2(pos)
    F_grav = payload_config.mass * a_grav

    # =================================================================
    # TOTAL ACCELERATION
    # =================================================================
    F_total = F_em + F_drag_vec + F_grav
    a_total = F_total / payload_config.mass

    # =================================================================
    # THERMAL HEATING (with vacuum tube support)
    # =================================================================
    if h < 100.0u"km"  # Only heat in atmosphere
        # Apply vacuum reduction if inside launcher tube
        if x_launcher < launcher_config.length
            q_dot = aerodynamic_heating(vel, h, payload_config.nose_radius) * launcher_config.vacuum_pressure_ratio
        else
            q_dot = aerodynamic_heating(vel, h, payload_config.nose_radius)
        end

        # Radiative cooling (Stefan-Boltzmann)
        σ = 5.670374419e-8u"W/(m^2*K^4)"
        T_amb = 250.0u"K"
        q_radiation = payload_config.emissivity * σ * payload_config.area *
                      (T^4 - T_amb^4)

        # Net heating rate
        q_net = q_dot * payload_config.area - q_radiation

        # Temperature rate
        dT_dt = q_net / (payload_config.mass * payload_config.specific_heat)
    else
        dT_dt = 0.0u"K/s"
    end

    # =================================================================
    # ELECTRICAL CIRCUIT DYNAMICS (RLC for each coil)
    # =================================================================
    dI_dt = zeros(n_coils) * 1.0u"A/s"
    dQ_dt = zeros(n_coils) * 1.0u"C/s"

    # Use x_launcher (computed above) for coil triggering
    # x_launcher = distance traveled from launch site (0 to 4000m)

    for i in 1:n_coils
        coil = launcher_config.coils[i]
        I = currents[i]
        Q = charges[i]

        # Applied voltage (triggered based on projectile position along launcher)
        V = coil_voltage(x_launcher, coil)

        # Circuit equation: L dI/dt = V - Q/C - RI
        dI_dt[i] = (V - Q / coil.capacitance - coil.resistance * I) / coil.inductance

        # Charge rate
        dQ_dt[i] = I
    end

    # =================================================================
    # ASSEMBLE DERIVATIVE VECTOR (strip units)
    # =================================================================
    # Position derivatives = velocity (m/s)
    for i in 1:3
        du[i] = ustrip(u"m/s", vel[i])
    end

    # Velocity derivatives = acceleration (m/s²)
    for i in 1:3
        du[3+i] = ustrip(u"m/s^2", a_total[i])
    end

    # Temperature derivative (K/s)
    du[7] = ustrip(u"K/s", dT_dt)

    # Current derivatives (A/s)
    for i in 1:n_coils
        du[7+i] = ustrip(u"A/s", dI_dt[i])
    end

    # Charge derivatives (C/s = A)
    for i in 1:n_coils
        du[7+n_coils+i] = ustrip(u"A", dQ_dt[i])
    end

    return nothing
end

"""
    simulate_trajectory(
        launcher_config,
        payload_config,
        mission_profile;
        tspan=(0.0u"s", 60.0u"s"),
        solver=Tsit5()
    )

Simulate complete launch trajectory.

# Arguments
- `launcher_config`: Launcher configuration
- `payload_config`: Payload configuration
- `mission_profile`: Mission parameters
- `tspan`: Time span for simulation [s]
- `solver`: ODE solver (default: Tsit5)

# Returns
- ODE solution object
"""
function simulate_trajectory(
    launcher_config::LauncherConfig,
    payload_config::PayloadConfig,
    mission_profile::MissionProfile;
    tspan=(0.0u"s", 60.0u"s"),
    solver=Tsit5(),
    saveat=nothing
)
    # Initial state (with units)
    u0_unitful = initial_state(launcher_config, payload_config, mission_profile)
    n_coils = launcher_config.num_coils

    # Convert to dimensionless
    u0_dimensionless, units_ref = strip_units_state(u0_unitful, n_coils)

    # Dimensionless time span
    tspan_dimensionless = (ustrip(u"s", tspan[1]), ustrip(u"s", tspan[2]))

    # Parameters (include units reference)
    p = (launcher_config, payload_config, mission_profile, units_ref)

    # Create ODE problem with dimensionless function
    prob = ODEProblem(trajectory_ode_dimensionless!, u0_dimensionless, tspan_dimensionless, p)

    # Solve
    if isnothing(saveat)
        sol = solve(prob, solver)
    else
        saveat_dimensionless = ustrip(u"s", saveat)
        sol = solve(prob, solver, saveat=saveat_dimensionless)
    end

    return sol
end

"""
    extract_trajectory_data(sol, n_coils, units_ref=nothing)

Extract human-readable trajectory data from solution.

# Arguments
- `sol`: ODE solution (dimensionless)
- `n_coils`: Number of coils
- `units_ref`: Optional units reference for reattaching units

# Returns
Named tuple with time series of key variables
"""
function extract_trajectory_data(sol, n_coils, units_ref=nothing)
    # Define default units if not provided
    if isnothing(units_ref)
        units_ref = (
            L = 1.0u"m",
            V = 1.0u"m/s",
            T = 1.0u"K",
            I = 1.0u"A",
            Q = 1.0u"C"
        )
    end

    times = sol.t * 1.0u"s"  # Reattach time units

    # Extract components at each time step (reattach units)
    positions = [SVector(sol[1,i], sol[2,i], sol[3,i]) * units_ref.L for i in 1:length(times)]
    velocities = [SVector(sol[4,i], sol[5,i], sol[6,i]) * units_ref.V for i in 1:length(times)]
    temperatures = [sol[7,i] * units_ref.T for i in 1:length(times)]

    # Derived quantities
    altitudes = [altitude_from_position(pos) for pos in positions]
    speeds = [norm(vel) for vel in velocities]
    mach_numbers = [mach_number(speeds[i], altitudes[i]) for i in 1:length(times)]

    return (
        times = times,
        positions = positions,
        velocities = velocities,
        temperatures = temperatures,
        altitudes = altitudes,
        speeds = speeds,
        mach_numbers = mach_numbers
    )
end

export PayloadConfig, MissionProfile
export initial_state, extract_state, launch_direction
export strip_units_state, reattach_units_state
export trajectory_ode_dimensionless!, simulate_trajectory, extract_trajectory_data

"""
    trajectory_ode!(du, u, p, t)

Backward-compatible wrapper for trajectory_ode_dimensionless!.
Automatically handles unitful state vectors.

Note: For direct ODE problem creation, use trajectory_ode_dimensionless! with
dimensionless state vectors. Or better yet, use simulate_trajectory() which
handles all unit conversions automatically.
"""
function trajectory_ode!(du, u, p, t)
    # Check if this is a dimensionless call (p has 4 elements including units_ref)
    if length(p) == 4
        # Already dimensionless, pass through
        return trajectory_ode_dimensionless!(du, u, p, t)
    else
        # Unitful call - needs conversion
        # This path is mainly for backward compatibility with old tests
        launcher_config, payload_config, mission_profile = p
        n_coils = launcher_config.num_coils

        # Convert to dimensionless if needed
        # Note: This assumes u and du are already properly set up
        # In practice, users should use simulate_trajectory instead
        units_ref = (
            L = 1.0u"m",
            V = 1.0u"m/s",
            T = 1.0u"K",
            I = 1.0u"A",
            Q = 1.0u"C"
        )

        p_new = (launcher_config, payload_config, mission_profile, units_ref)
        return trajectory_ode_dimensionless!(du, u, p_new, t)
    end
end

export trajectory_ode!

# ============================================================================
# SIMPLIFIED TRAJECTORY MODEL (for Monte Carlo / SDE simulations)
# ============================================================================

"""
    TrajectoryLog

Stores logged data during simulation for analysis.
Fields contain vectors of values at each logged timestep.
"""
mutable struct TrajectoryLog
    t::Vector               # Time with units
    pos::Vector             # Position with units
    vel::Vector             # Velocity with units
    altitude::Vector        # Altitude with units
    speed::Vector           # Speed with units
    F_em::Vector            # EM force with units
    F_drag::Vector          # Drag force with units
    F_grav::Vector          # Gravity force with units
    temperature::Vector     # Temperature with units
    x_launcher::Vector      # Position along launcher with units
    mach::Vector{Float64}   # Mach number (dimensionless)
end

TrajectoryLog() = TrajectoryLog(
    [], [], [], [], [], [], [], [], [], [], Float64[]
)

"""
    create_logging_callback(log::TrajectoryLog, p_params)

Create a DifferentialEquations.jl callback that logs trajectory data.
Uses PeriodicCallback for compatibility with unitful time.
"""
function create_logging_callback(log::TrajectoryLog, p_params)
    function log_state(integrator)
        # Extract parameters (handle both old 10-param and new 13-param versions)
        if length(p_params) >= 13
            F_avg, L_launcher, launch_dir, m_payload, A_payload, c_p, ε, R_nose, launch_pos, R_Earth, vacuum_pressure_ratio, drag_reduction_method, d_capsule = p_params
        else
            F_avg, L_launcher, launch_dir, m_payload, A_payload, c_p, ε, R_nose, launch_pos, R_Earth = p_params
            vacuum_pressure_ratio = 1.0
        end

        # Extract state from integrator
        u = integrator.u
        t = integrator.t

        # Handle unitful state
        pos = SVector(ustrip(u"m", u[1]), ustrip(u"m", u[2]), ustrip(u"m", u[3]))
        vel = SVector(ustrip(u"m/s", u[4]), ustrip(u"m/s", u[5]), ustrip(u"m/s", u[6]))
        T = ustrip(u"K", u[7])

        # Strip units from parameters
        R_Earth_m = ustrip(u"m", R_Earth)
        L_launcher_m = ustrip(u"m", L_launcher)
        F_avg_n = ustrip(u"N", F_avg)
        A_payload_m2 = ustrip(u"m^2", A_payload)
        m_payload_kg = ustrip(u"kg", m_payload)

        # Calculate derived quantities (all unitless now)
        r = norm(pos)
        h = r - R_Earth_m
        launch_pos_m = SVector(ustrip(u"m", launch_pos[1]), ustrip(u"m", launch_pos[2]), ustrip(u"m", launch_pos[3]))
        x_launcher = norm(pos - launch_pos_m)
        speed = norm(vel)

        # Forces (recalculate from state) - unitless
        # EM force
        if x_launcher < L_launcher_m
            F_em = F_avg_n * launch_dir
        else
            F_em = SVector(0.0, 0.0, 0.0)
        end

        # Drag
        if h < 100000.0
            ρ_0 = 1.225
            H = 8500.0
            h_clamped = max(h, 0.0)
            ρ = ρ_0 * exp(-h_clamped / H)

            # Apply vacuum reduction if inside launcher
            if x_launcher < L_launcher_m
                ρ = ρ * vacuum_pressure_ratio
            end

            C_d = 0.5
            if speed > 0.0
                F_drag = -0.5 * ρ * C_d * A_payload_m2 * speed^2 * (vel / speed)
            else
                F_drag = SVector(0.0, 0.0, 0.0)
            end
        else
            F_drag = SVector(0.0, 0.0, 0.0)
        end

        # Gravity
        μ = 3.986004418e14
        if r > 0.0
            a_grav = -μ / r^3 * pos
            F_grav = m_payload_kg * a_grav
        else
            F_grav = SVector(0.0, 0.0, 0.0)
        end

        # Calculate Mach number for logging
        γ = 1.4
        R_gas = 287.0
        T_amb = 288.15 - 6.5e-3 * h
        T_amb = max(T_amb, 216.65)
        a_sound = sqrt(γ * R_gas * T_amb)
        mach = speed / a_sound

        # Log everything with units restored
        t_val = ustrip(u"s", t)
        push!(log.t, t_val * 1.0u"s")
        push!(log.pos, pos * 1.0u"m")
        push!(log.vel, vel * 1.0u"m/s")
        push!(log.altitude, h * 1.0u"m")
        push!(log.speed, speed * 1.0u"m/s")
        push!(log.F_em, F_em * 1.0u"N")
        push!(log.F_drag, F_drag * 1.0u"N")
        push!(log.F_grav, F_grav * 1.0u"N")
        push!(log.temperature, T * 1.0u"K")
        push!(log.x_launcher, x_launcher * 1.0u"m")
        push!(log.mach, mach)

        return false  # Don't terminate
    end

    # Use PeriodicCallback instead of PresetTimeCallback for unitful time compatibility
    # Log every 0.1 seconds
    return PeriodicCallback(log_state, 0.1u"s")
end

"""
    trajectory_ode_simplified!(du, u, p, t)

Simplified trajectory ODE using constant average force model.
Much more stable for SDE/Monte Carlo simulations than full RLC model.

## State Vector
u = [x, y, z, vx, vy, vz, T] (7 variables)
- Position: [x, y, z] in meters (Earth-centered inertial frame)
- Velocity: [vx, vy, vz] in m/s
- Temperature: T in Kelvin

## Parameters
p = (F_avg, L_launcher, launch_dir, m_payload, A_payload, c_p, ε, R_nose, launch_pos, R_Earth)
- F_avg: Average EM force [N]
- L_launcher: Launcher tube length [m]
- launch_dir: Unit vector pointing in launch direction
- m_payload: Payload mass [kg]
- A_payload: Cross-sectional area [m²]
- c_p: Specific heat capacity [J/(kg·K)]
- ε: Emissivity [dimensionless, 0-1]
- R_nose: Nose radius [m]
- launch_pos: Initial position vector [m]
- R_Earth: Earth radius [m]

## Physics
Implements simplified point-mass dynamics:
1. Constant electromagnetic force (inside launcher only)
2. Exponential atmosphere drag
3. Point-mass gravity (two-body problem)
4. Aerodynamic heating (Sutton-Graves formula)

# TODO: Add unit enforcement with Unitful
# TODO: Add energy conservation checking
"""
function trajectory_ode_simplified!(du, u, p, t)
    # =================================================================
    # UNPACK PARAMETERS
    # =================================================================
    F_avg, L_launcher, launch_dir, m_payload, A_payload, c_p, ε, R_nose, launch_pos, R_Earth, vacuum_pressure_ratio, drag_reduction_method, d_capsule = p

    # =================================================================
    # EXTRACT STATE
    # =================================================================
    pos = SVector(u[1], u[2], u[3])  # Position [m]
    vel = SVector(u[4], u[5], u[6])  # Velocity [m/s]
    T = u[7]                          # Temperature [K]

    # =================================================================
    # DERIVED QUANTITIES
    # =================================================================
    r = norm(pos)                     # Distance from Earth center [m]
    h = r - R_Earth                   # Altitude above surface [m]
    x_launcher = norm(pos - launch_pos)  # Distance along launcher [m]
    v_mag = norm(vel)                 # Speed [m/s]

    # Physical assertions
    @assert r > 0.0u"m" "Position magnitude must be positive (r=$r)"
    @assert T > 0.0u"K" "Temperature must be positive (T=$T K)"
    @assert isfinite(ustrip(u"m", r)) && isfinite(ustrip(u"m/s", v_mag)) && isfinite(ustrip(u"K", T)) "State contains non-finite values"

    # =================================================================
    # ELECTROMAGNETIC FORCE
    # =================================================================
    # Physics: F_em = constant force while x < L_launcher
    # This is a simplified model; real coilgun has position-dependent force
    # See em_acceleration.jl for detailed electromagnetic theory
    if x_launcher < L_launcher
        F_em = F_avg * launch_dir  # [N]
    else
        F_em = SVector(0.0u"N", 0.0u"N", 0.0u"N")
    end

    # =================================================================
    # ATMOSPHERIC DRAG FORCE
    # =================================================================
    # Physics: Exponential atmosphere model
    #   ρ(h) = ρ₀ exp(-h/H)
    # where:
    #   ρ₀ = 1.225 kg/m³  (sea level density)
    #   H = 8500 m         (scale height)
    #
    # Drag equation:
    #   F_drag = -½ ρ C_d A v² (v̂)
    # where:
    #   C_d = drag coefficient (0.5 for streamlined projectile)
    #   A = cross-sectional area [m²]
    #   v̂ = velocity direction (unit vector)
    #
    if h < 100000.0u"m"  # Only apply drag below 100 km
        ρ_0 = 1.225u"kg/m^3"      # Sea level density
        H = 8500.0u"m"             # Atmospheric scale height
        # Clamp altitude to 0 if below sea level (don't model underground)
        h_clamped = max(h, 0.0u"m")
        ρ = ρ_0 * exp(-ustrip(u"m", h_clamped) / ustrip(u"m", H))  # Exponential atmosphere

        # Apply vacuum tube pressure reduction while inside launcher
        if x_launcher < L_launcher
            ρ = ρ * vacuum_pressure_ratio  # Reduced density in vacuum tube
        end

        @assert ρ >= 0.0u"kg/m^3" "Atmospheric density must be non-negative"
        @assert ρ <= ρ_0 "Atmospheric density should not exceed sea level value (accounting for vacuum reduction)"

        # Base drag coefficient (streamlined body)
        C_d_base = 0.5  # Dimensionless

        # Calculate Mach number for drag reduction
        # Speed of sound: a = √(γ R T) where γ=1.4, R=287 J/(kg·K)
        γ = 1.4  # Dimensionless
        R_gas = 287.0u"J/(kg*K)"  # Specific gas constant for air
        # Use ambient temperature from standard atmosphere
        T_amb = 288.15u"K" - 6.5e-3u"K/m" * h_clamped  # Linear temperature lapse
        T_amb = max(T_amb, 216.65u"K")  # Min temperature (tropopause)
        a_sound = sqrt(γ * R_gas * T_amb)
        mach = v_mag / a_sound  # Dimensionless

        # Apply passive drag reduction based on Mach number and altitude
        # Strip units for drag reduction factor calculation (it's empirical)
        C_d = effective_drag_coefficient(C_d_base, drag_reduction_method, ustrip(mach), ustrip(u"m", h_clamped))

        if v_mag > 0.0u"m/s"
            # Drag force opposes velocity: F_drag = -½ ρ C_d A v² v̂
            F_drag = -0.5 * ρ * C_d * A_payload * v_mag^2 * (vel / v_mag)  # [N]
        else
            F_drag = SVector(0.0u"N", 0.0u"N", 0.0u"N")
        end
    else
        F_drag = SVector(0.0u"N", 0.0u"N", 0.0u"N")  # No drag in space
    end

    # =================================================================
    # GRAVITATIONAL FORCE
    # =================================================================
    # Physics: Newton's law of gravitation (two-body, point mass)
    #   F_grav = -GMm/r² r̂ = -μm/r² r̂
    # where:
    #   G = gravitational constant
    #   M = Earth mass
    #   μ = GM = 3.986004418×10¹⁴ m³/s² (Earth's gravitational parameter)
    #   r̂ = r/|r| (unit vector)
    #
    # This simplifies to:
    #   a_grav = -μ/r³ r
    #
    # TODO: Add J2 perturbation for more accuracy (see gravity.jl)
    μ = 3.986004418e14u"m^3/s^2"  # Earth's gravitational parameter

    if r > 0.0u"m"
        # Gravitational acceleration: a = -μ/r³ r
        a_grav = -μ / r^3 * pos  # [m/s²]

        # Sanity check: gravity should be ~9.8 m/s² at surface
        g_mag = norm(a_grav)
        if h < 10000.0u"m"  # Within 10 km of surface
            @assert g_mag > 9.0u"m/s^2" && g_mag < 11.0u"m/s^2" "Gravity unreasonable near surface (g=$g_mag)"
        end
    else
        a_grav = SVector(0.0u"m/s^2", 0.0u"m/s^2", 0.0u"m/s^2")
    end

    F_grav = m_payload * a_grav  # [N]

    # =================================================================
    # TOTAL FORCE AND ACCELERATION
    # =================================================================
    # Newton's second law: F = ma  =>  a = F/m
    F_total = F_em + F_drag + F_grav  # [N]
    a_total = F_total / m_payload      # [m/s²]

    # Sanity check on total acceleration
    # Allow up to 20,000,000 m/s² (~2,000,000 g) for extreme launcher configurations
    # (LEO-capable launchers with very high exit velocities and short tracks need these high accelerations)
    a_mag = norm(a_total)
    @assert a_mag < 20000000.0u"m/s^2" "Acceleration unreasonably high (a=$a_mag)"

    # =================================================================
    # AERODYNAMIC HEATING
    # =================================================================
    # Physics: Sutton-Graves formula for stagnation point heat flux
    #   q̇ = K √(ρ/R_nose) v³
    # where:
    #   K = 1.83×10⁻⁴ (empirical constant)
    #   ρ = atmospheric density [kg/m³]
    #   R_nose = nose radius [m]
    #   v = velocity [m/s]
    #
    # Radiative cooling (Stefan-Boltzmann law):
    #   q_rad = ε σ A (T⁴ - T_amb⁴)
    # where:
    #   ε = emissivity
    #   σ = 5.67×10⁻⁸ W/(m²·K⁴) (Stefan-Boltzmann constant)
    #   A = surface area [m²]
    #   T_amb = ambient temperature [K]
    #
    # Energy balance:
    #   m c_p dT/dt = q_heating - q_cooling
    #
    # TODO: Add conduction to internal structure
    # TODO: Use more sophisticated heating model for subsonic/transonic
    if h < 100000.0u"m" && v_mag > 100.0u"m/s"  # Only heat in atmosphere at high speed
        ρ_0 = 1.225u"kg/m^3"
        H = 8500.0u"m"
        ρ = ρ_0 * exp(-ustrip(u"m", h) / ustrip(u"m", H))  # [kg/m³]

        # Apply vacuum tube pressure reduction while inside launcher
        if x_launcher < L_launcher
            ρ = ρ * vacuum_pressure_ratio  # Reduced heating in vacuum tube
        end

        # Sutton-Graves heating
        K = 1.83e-4u"kg^0.5/m"  # Empirical constant
        q_dot = K * sqrt(ρ / R_nose) * v_mag^3  # [W/m²] Heat flux

        # Radiative cooling
        σ = 5.670374419e-8u"W/(m^2*K^4)"  # Stefan-Boltzmann constant
        T_amb = 250.0u"K"       # Ambient temperature (high altitude)

        @assert T >= T_amb "Payload cooler than ambient?"

        q_radiation = ε * σ * A_payload * (T^4 - T_amb^4)  # [W] Radiated power

        # Net heating rate
        q_net = q_dot * A_payload - q_radiation  # [W]

        # Temperature rate of change
        dT_dt = q_net / (m_payload * c_p)  # [K/s]
    else
        dT_dt = 0.0u"K/s"  # No heating in space or at low speeds
    end

    # =================================================================
    # ASSEMBLE DERIVATIVE VECTOR
    # =================================================================
    # State equations:
    #   dr/dt = v           (position derivative is velocity)
    #   dv/dt = a = F/m     (velocity derivative is acceleration)
    #   dT/dt = q_net/(mc_p) (temperature derivative from energy balance)

    du[1] = vel[1]      # dx/dt = vx [m/s]
    du[2] = vel[2]      # dy/dt = vy [m/s]
    du[3] = vel[3]      # dz/dt = vz [m/s]
    du[4] = a_total[1]  # dvx/dt = ax [m/s²]
    du[5] = a_total[2]  # dvy/dt = ay [m/s²]
    du[6] = a_total[3]  # dvz/dt = az [m/s²]
    du[7] = dT_dt       # dT/dt [K/s]

    return nothing
end

"""
    create_simplified_sde_problem(v_target, L_launcher, m_payload, elevation, azimuth; kwargs...)

Create SDE problem using simplified constant-force model (no RLC circuits).
Much more stable and faster than full electrical model.

Returns: SDEProblem ready to solve, TrajectoryLog for recording data
"""
function create_simplified_sde_problem(
    v_target, L_launcher, m_payload, elevation, azimuth;
    tspan=(0.0u"s", 180.0u"s"), noise_level=0.05, enable_logging=false, vacuum_pressure_ratio=1.0,
    drag_reduction_method=NoDragReduction(), capsule_diameter=0.4u"m"
)
    # Keep all quantities with units!
    v = v_target
    L = L_launcher
    m_base = m_payload
    θ_elev = ustrip(u"rad", elevation)  # Convert to radians (dimensionless for trig)
    θ_azim = ustrip(u"rad", azimuth)
    d_capsule = capsule_diameter

    # Calculate added mass from drag reduction hardware
    m_added = added_mass(drag_reduction_method, capsule_diameter, m_payload)
    m = m_base + m_added

    # Calculate required force (for total mass including drag reduction hardware)
    a_avg = v^2 / (2 * L)
    F_avg = m * a_avg

    # Earth parameters
    R_E = 6.3781370e6u"m"

    # Initial position
    launch_pos = SVector(R_E, 0.0u"m", 0.0u"m")

    # Launch direction (proper local coordinates) - dimensionless unit vectors
    radial_dir = SVector(1.0, 0.0, 0.0)
    east_dir = SVector(0.0, 1.0, 0.0)
    north_dir = SVector(0.0, 0.0, 1.0)
    horizontal_component = cos(θ_elev)
    vertical_component = sin(θ_elev)
    horizontal_dir = cos(θ_azim) * north_dir + sin(θ_azim) * east_dir
    launch_dir = horizontal_component * horizontal_dir + vertical_component * radial_dir
    launch_dir = launch_dir / norm(launch_dir)

    # Initial state (with units!)
    u0 = [launch_pos[1], launch_pos[2], launch_pos[3], 0.0u"m/s", 0.0u"m/s", 0.0u"m/s", 300.0u"K"]

    # Payload parameters (with units!)
    A_payload = 0.1u"m^2"
    c_p = 900.0u"J/(kg*K)"
    ε = 0.8  # Dimensionless emissivity
    R_nose = 0.1u"m"

    # Parameters tuple (includes drag reduction method and capsule diameter)
    p = (F_avg, L, launch_dir, m, A_payload, c_p, ε, R_nose, launch_pos, R_E, vacuum_pressure_ratio, drag_reduction_method, d_capsule)

    # Noise function for SDE (handle units correctly)
    # In SDEs: du = f(u,p,t)*dt + g(u,p,t)*dW where dW ~ N(0, dt)
    #
    # **With UNITFUL time** (dt has units s):
    #   dW has units sqrt(s), so g must have units [state]/sqrt(s)
    #   e.g., position noise: m/s^0.5
    #
    # **With UNITLESS time** (dt is Float64 seconds):
    #   dW is unitless, so g must have units [state]
    #   e.g., position noise: m
    #
    # This function works for UNITLESS time (compatible with SDE solvers)
    function noise!(du, u, p, t)
        # Position noise: meters (when time is unitless)
        du[1:3] .= 1.0u"m"

        # Velocity noise scales with current velocity
        v_mag = sqrt(u[4]^2 + u[5]^2 + u[6]^2)
        v_noise = max(1.0u"m/s", v_mag * noise_level)
        du[4:6] .= v_noise

        # Temperature noise: Kelvin
        du[7] = 1.0u"K"
        return nothing
    end

    # Time span (already has units)
    t0 = tspan[1]
    tf = tspan[2]

    # Create logging if requested
    log = enable_logging ? TrajectoryLog() : nothing
    callback = enable_logging ? create_logging_callback(log, p) : nothing

    # Create appropriate problem type based on noise_level
    # IMPORTANT: For unitful state vectors with mixed dimensions (m, m/s, K):
    #
    # **ODE solvers (deterministic, noise_level=0.0)**:
    #   - Use fixed-step solvers: Euler(), RK4() with dt in seconds (e.g., dt=0.01u"s")
    #   - Do NOT use saveat (interpolation fails with mixed units)
    #   - Do NOT use callbacks with preset times (unitful time incompatibility)
    #   - Recommended: solve(prob, Euler(), dt=0.01u"s", adaptive=false, dense=false)
    #
    # **SDE solvers (stochastic, noise_level>0.0)**:
    #   - SDE solvers have limited unitful support due to sqrt(dt) operations
    #   - TIME must be unitless but STATE remains unitful
    #   - Use ustrip(u"s", tspan) for time, keep u0 unitful
    #   - Recommended: solve(prob, EM(), dt=0.01, adaptive=false, dense=false)
    #
    # For this function: Problem created with unitful time for ODE compatibility.
    # For SDE: User should convert time to unitless before solving.

    if noise_level == 0.0
        # Deterministic case: create ODE problem
        prob = ODEProblem(trajectory_ode_simplified!, u0, (t0, tf), p)
    else
        # Stochastic case: create SDE problem
        prob = SDEProblem(trajectory_ode_simplified!, noise!, u0, (t0, tf), p)
    end

    return prob, log, callback
end

"""
    NOTE: SDE solvers have fundamental limitations with unitful types.

Julia's DifferentialEquations.jl SDE solvers have limited support for:
1. Unitful time (sqrt(dt) creates unit incompatibilities)
2. Heterogeneous unitful state vectors (Vector{Quantity{Float64}} with mixed dimensions)

**Recommended approach for stochastic simulations**:
- For uncertainty quantification, use `noise_level=0.0` and solve deterministically with ODE solvers
- Run Monte Carlo simulations by varying initial conditions instead of using SDE noise
- ODE solvers (Euler, RK4) fully support unitful state with mixed dimensions

**If you need SDE functionality**:
- Convert the entire problem to unitless (strip all units)
- Add units back to results after solving
- See stochastic_trajectory.jl for unitless SDE implementations

This is a known limitation of the DifferentialEquations.jl SDE infrastructure.
"""

export TrajectoryLog, create_logging_callback
export trajectory_ode_simplified!, create_simplified_sde_problem
