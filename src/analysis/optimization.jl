"""
# Optimization Module

Parameter optimization for electromagnetic launcher system.

## Optimization Objectives

1. **Maximize final velocity** for given energy budget
2. **Minimize energy consumption** for target velocity
3. **Minimize structural loads** (g-forces)
4. **Minimize thermal loading** (peak temperature)
5. **Multi-objective optimization** balancing all factors

## Design Variables

- Coil parameters: `N_coils`, spacing, inductance, capacitance, voltage
- Coil timing: `t_fire[i]` for each coil
- Launch angle: elevation and azimuth
- Payload mass allocation

## Constraints

- Structural: `a_max < a_structural_limit`
- Thermal: `T_max < T_thermal_limit`
- Energy: `E_total < E_budget`
- Velocity: `v_final ≥ v_target`
"""

using Unitful
using Optim
using Optimization
using ForwardDiff
using LinearAlgebra
using StaticArrays
using DifferentialEquations

"""
    OptimizationObjective

Enum for optimization objectives.
"""
@enum OptimizationObjective begin
    MAXIMIZE_VELOCITY
    MINIMIZE_ENERGY
    MINIMIZE_GLOAD
    MINIMIZE_THERMAL
    MULTI_OBJECTIVE
end

"""
    objective_function(params, fixed_config, objective_type)

Generic objective function for optimization.

# Arguments
- `params`: Vector of design parameters to optimize
- `fixed_config`: Fixed configuration parameters
- `objective_type`: Type of optimization objective

# Returns
- Objective value (to be minimized)
"""
function objective_function(
    params::Vector,
    fixed_config,
    objective_type::OptimizationObjective
)
    # Unpack parameters and create launcher configuration
    # This depends on what we're optimizing

    # Run simulation
    # sol = simulate_trajectory(launcher_config, payload_config, mission_profile)

    # Calculate objective
    if objective_type == MAXIMIZE_VELOCITY
        # Return negative velocity (since we minimize)
        # return -final_velocity
    elseif objective_type == MINIMIZE_ENERGY
        # return total_energy
    elseif objective_type == MINIMIZE_GLOAD
        # return max_gload
    elseif objective_type == MINIMIZE_THERMAL
        # return max_temperature
    else  # MULTI_OBJECTIVE
        # return weighted_sum
    end

    # Placeholder
    return 0.0
end

"""
    optimize_coil_timing(launcher_config, payload_config, mission_profile)

Optimize coil firing times for maximum final velocity.

# Arguments
- `launcher_config`: Initial launcher configuration
- `payload_config`: Payload configuration
- `mission_profile`: Mission profile

# Returns
- Optimized coil timing array [s]
"""
function optimize_coil_timing(
    launcher_config::LauncherConfig,
    payload_config::PayloadConfig,
    mission_profile::MissionProfile
)
    n_coils = launcher_config.num_coils

    # Initial guess: uniform spacing in time
    t_max = 2.0  # seconds (rough estimate for launcher transit time)
    t_initial = range(0.0, t_max, length=n_coils) |> collect

    # Objective: maximize exit velocity
    function obj(t_fire)
        # Update launcher config with new timing
        # Run simulation
        # Extract final velocity
        # Return negative velocity (minimization)

        # Placeholder
        return sum(t_fire)  # Replace with actual simulation
    end

    # Bounds: each coil must fire after previous one
    lower = zeros(n_coils)
    upper = fill(t_max, n_coils)

    # Optimize using BFGS or similar
    # result = optimize(obj, lower, upper, t_initial, Fminbox(BFGS()))

    # Placeholder
    return t_initial * 1.0u"s"
end

"""
    optimize_launch_angle_min_drag(
        v_target, L_launcher, m_payload;
        elevation_range=(20.0u"°", 80.0u"°"),
        azimuth=90.0u"°",
        vacuum_pressure_ratio=1.0,
        verbose=false
    )

Optimize launch elevation angle to minimize atmospheric drag losses.

Higher angles get out of the atmosphere faster but have longer flight paths.
Lower angles have shorter paths but spend more time in dense atmosphere.

# Arguments
- `v_target`: Target exit velocity [m/s]
- `L_launcher`: Launcher length [m]
- `m_payload`: Payload mass [kg]
- `elevation_range`: Range of elevations to search [degrees]
- `azimuth`: Fixed azimuth angle [degrees] (default: 90° = East)
- `vacuum_pressure_ratio`: Vacuum level (1.0 = atmospheric, 0.001 = high vacuum)
- `verbose`: Print optimization progress

# Returns
- Named tuple with optimal elevation, final velocity, drag loss, and altitude
"""
function optimize_launch_angle_min_drag(
    v_target, L_launcher, m_payload;
    elevation_range=(20.0u"°", 80.0u"°"),
    azimuth=90.0u"°",
    vacuum_pressure_ratio=1.0,
    verbose=false
)
    # Strip units for optimization
    function strip_unit(x)
        if hasfield(typeof(x), :val)
            return Float64(x.val)
        else
            return Float64(x)
        end
    end

    elev_min = strip_unit(elevation_range[1])
    elev_max = strip_unit(elevation_range[2])

    # Objective: minimize drag loss (maximize final velocity)
    function objective(elev_deg)
        try
            # Create problem with this elevation
            # noise_level=0.0 creates ODEProblem for deterministic simulation
            prob, log, callback = create_simplified_sde_problem(
                v_target, L_launcher, m_payload,
                elev_deg[1] * 1.0u"°", azimuth;
                tspan=(0.0u"s", 120.0u"s"),
                noise_level=0.0,
                enable_logging=true,
                vacuum_pressure_ratio=vacuum_pressure_ratio
            )

            # Solve with ODE solver (Euler supports unitful time)
            # Use fixed timestep without saveat (unitful compatibility)
            sol = solve(prob, Euler(), dt=0.1u"s", adaptive=false, dense=false)

            if sol.retcode != SciMLBase.ReturnCode.Success
                return 1e10  # Penalty for failed solution
            end

            # Extract final velocity (what we want to maximize)
            u_final = sol.u[end]
            vel_final = SVector(u_final[4], u_final[5], u_final[6])
            v_final = norm(vel_final)

            # Return negative velocity (we minimize, but want to maximize velocity)
            # Strip units for optimizer (Optim.jl requires Float64)
            return -ustrip(u"m/s", v_final)

        catch e
            if verbose
                println("  Error at elevation $(elev_deg[1])°: $e")
            end
            return 1e10  # Penalty for errors
        end
    end

    # Initial guess: 45° (compromise between altitude gain and path length)
    x0 = [45.0]

    if verbose
        println("Optimizing launch angle to minimize drag...")
        println("  Elevation range: $(elev_min)° to $(elev_max)°")
        println("  Azimuth (fixed): $(strip_unit(azimuth))°")
        println("  Vacuum ratio: $vacuum_pressure_ratio")
        println()
    end

    # Optimize using bounded Nelder-Mead (derivative-free, good for noisy objectives)
    result = optimize(
        objective,
        [elev_min],
        [elev_max],
        x0,
        Fminbox(NelderMead()),
        Optim.Options(
            iterations=20,
            show_trace=verbose
        )
    )

    optimal_elevation = result.minimizer[1]
    optimal_obj = result.minimum

    # Run simulation at optimal angle to get full results
    # noise_level=0.0 creates ODEProblem for deterministic simulation
    prob, log, callback = create_simplified_sde_problem(
        v_target, L_launcher, m_payload,
        optimal_elevation * 1.0u"°", azimuth;
        tspan=(0.0u"s", 120.0u"s"),
        noise_level=0.0,
        enable_logging=true,
        vacuum_pressure_ratio=vacuum_pressure_ratio
    )

    # Solve with ODE solver (Euler supports unitful time)
    sol = solve(prob, Euler(), dt=0.1u"s", adaptive=false, dense=false)

    # Extract results
    R_Earth = 6.3781370e6u"m"
    u_final = sol.u[end]
    pos_final = SVector(u_final[1], u_final[2], u_final[3])
    vel_final = SVector(u_final[4], u_final[5], u_final[6])

    r_final = norm(pos_final)
    v_final = norm(vel_final)
    h_final = r_final - R_Earth

    # Calculate drag loss
    drag_loss = 0.0u"m/s"
    if !isnothing(log) && length(log.t) > 10
        em_on_idx = findfirst(x -> norm(x) < 1.0, log.F_em)
        if !isnothing(em_on_idx)
            v_exit = log.speed[em_on_idx]
            drag_loss = v_exit - v_final
        end
    end

    if verbose
        println()
        println("Optimization complete:")
        println("  Optimal elevation: $(round(optimal_elevation, digits=2))°")
        println("  Final velocity: $(round(ustrip(u"m/s", v_final), digits=1)) m/s")
        println("  Final altitude: $(round(ustrip(u"km", h_final), digits=2)) km")
        println("  Drag loss: $(round(ustrip(u"m/s", drag_loss), digits=1)) m/s")
    end

    return (
        elevation = optimal_elevation * 1.0u"°",
        azimuth = azimuth,
        v_final = v_final,
        h_final = h_final,
        drag_loss = drag_loss,
        converged = Optim.converged(result)
    )
end

"""
    optimize_launch_angle(launcher_config, payload_config, target_altitude)

Optimize launch elevation and azimuth for orbital insertion.

DEPRECATED: Use optimize_launch_angle_min_drag for simplified interface.

# Arguments
- `launcher_config`: Launcher configuration
- `payload_config`: Payload configuration
- `target_altitude`: Target orbital altitude [m]

# Returns
- Optimal elevation [degrees], azimuth [degrees]
"""
function optimize_launch_angle(
    launcher_config::LauncherConfig,
    payload_config::PayloadConfig,
    target_altitude
)
    # Objective: minimize Δv required for circularization
    function obj(angles)
        elevation, azimuth = angles

        # Create mission profile with these angles
        # Simulate trajectory
        # Calculate Δv for orbit circularization
        # Return Δv

        # Placeholder
        return (elevation - 45.0)^2 + (azimuth - 90.0)^2
    end

    # Initial guess
    x0 = [45.0, 90.0]  # 45° elevation, 90° azimuth (East)

    # Bounds
    lower = [20.0, 0.0]
    upper = [60.0, 180.0]

    # Optimize
    # result = optimize(obj, lower, upper, x0, Fminbox(BFGS()))

    # Placeholder
    return (45.0u"°", 90.0u"°")
end

"""
    optimize_payload_mass(launcher_config, target_velocity, energy_budget)

Find maximum payload mass achievable for target velocity and energy budget.

# Arguments
- `launcher_config`: Launcher configuration
- `target_velocity`: Target final velocity [m/s]
- `energy_budget`: Available electrical energy [J]

# Returns
- Maximum payload mass [kg]
"""
function optimize_payload_mass(
    launcher_config::LauncherConfig,
    target_velocity,
    energy_budget
)
    # Binary search or optimization to find max mass
    # that still reaches target velocity with given energy

    function can_reach_velocity(mass)
        # Create payload config with this mass
        # Run simulation
        # Check if v_final >= target_velocity

        # Placeholder
        return true
    end

    # Binary search
    m_min = 1.0u"kg"
    m_max = 1000.0u"kg"

    # Placeholder
    return 50.0u"kg"
end

"""
    design_trade_study(
        launcher_length_range,
        n_coils_range,
        payload_mass_range,
        target_velocity
    )

Perform parametric trade study across design space.

# Returns
- DataFrame with performance metrics for each configuration
"""
function design_trade_study(
    launcher_length_range,
    n_coils_range,
    payload_mass_range,
    target_velocity
)
    # Placeholder for comprehensive trade study
    # Would iterate through parameter combinations
    # Run simulations for each
    # Collect performance metrics

    results = Dict(
        "launcher_length" => Float64[],
        "n_coils" => Int[],
        "payload_mass" => Float64[],
        "final_velocity" => Float64[],
        "energy_required" => Float64[],
        "max_gload" => Float64[],
        "max_temperature" => Float64[]
    )

    return results
end

export OptimizationObjective
export MAXIMIZE_VELOCITY, MINIMIZE_ENERGY, MINIMIZE_GLOAD, MINIMIZE_THERMAL, MULTI_OBJECTIVE
export objective_function
export optimize_coil_timing, optimize_launch_angle, optimize_launch_angle_min_drag, optimize_payload_mass
export design_trade_study
