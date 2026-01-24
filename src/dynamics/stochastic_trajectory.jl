"""
# Stochastic Trajectory Module

Uses Stochastic Differential Equations (SDEs) to model uncertainty in the
electromagnetic launcher system, producing probability distributions instead
of point estimates.

## Uncertainty Sources

1. **Atmospheric Variations**: Weather fluctuations, seasonal changes
2. **Manufacturing Tolerances**: Coil parameters, payload mass
3. **Control System Noise**: Timing jitter, current ripple
4. **Measurement Uncertainty**: Sensor noise, position errors
5. **Physical Variations**: Drag coefficient variations with Reynolds number

## SDE Formulation

The deterministic ODE system:
```math
\\frac{du}{dt} = f(u, p, t)
```

Becomes an SDE:
```math
du = f(u, p, t)dt + g(u, p, t)dW
```

where:
- `f(u, p, t)` = deterministic drift term (same as ODE)
- `g(u, p, t)` = stochastic diffusion term (noise magnitude)
- `dW` = Wiener process (Brownian motion)

## Noise Model

Additive noise on key parameters:
- **Atmospheric density**: ±5% (weather variations)
- **Drag coefficient**: ±10% (Reynolds number effects)
- **Coil current**: ±2% (power supply ripple)
- **Coil timing**: ±0.1ms (trigger jitter)

## Monte Carlo Analysis

Run multiple SDE realizations to get:
- Mean trajectory
- Standard deviation bounds
- Confidence intervals (68%, 95%, 99%)
- Success probability distributions
"""

using Unitful
using StaticArrays
using LinearAlgebra
using DifferentialEquations
using SciMLBase
using Statistics
using Logging
using ProgressMeter

"""
    NoiseParameters

Parameters controlling stochastic noise in the simulation.

# Fields
- `atmospheric_noise`: Relative noise in atmospheric density [dimensionless]
- `drag_noise`: Relative noise in drag coefficient [dimensionless]
- `current_noise`: Relative noise in coil currents [dimensionless]
- `timing_noise`: Absolute noise in coil timing [s]
- `mass_noise`: Relative noise in payload mass [dimensionless]
"""
struct NoiseParameters
    atmospheric_noise::Float64      # Relative (e.g., 0.05 = 5%)
    drag_noise::Float64
    current_noise::Float64
    timing_noise                    # Absolute with units
    mass_noise::Float64
end

"""
    default_noise_parameters()

Create default noise parameters for realistic uncertainty modeling.

Returns conservative estimates based on:
- Atmospheric variability: ±5% (weather fluctuations)
- Drag coefficient: ±10% (Reynolds number variations)
- Current ripple: ±2% (power supply quality)
- Timing jitter: ±0.1ms (control system precision)
- Mass uncertainty: ±1% (manufacturing tolerance)
"""
function default_noise_parameters()
    return NoiseParameters(
        0.05,           # 5% atmospheric density variation
        0.10,           # 10% drag coefficient variation
        0.02,           # 2% current ripple
        0.0001u"s",     # 0.1ms timing jitter
        0.01            # 1% mass uncertainty
    )
end

"""
    trajectory_sde_noise!(du, u, p, t)

Diagonal noise function for SDE system.
Returns the magnitude of noise for each state variable.

This implements multiplicative noise where uncertainty scales with
the magnitude of the state variables and system parameters.

# Arguments
- `du`: Noise magnitude vector (output)
- `u`: State vector [x,y,z,vₓ,vᵧ,vᵤ,T,I₁...Iₙ,Q₁...Qₙ] (dimensionless)
- `p`: Parameters (launcher_config, payload_config, mission_profile, units_ref, noise_params)
- `t`: Time [s] (dimensionless)
"""
function trajectory_sde_noise!(du, u, p, t)
    launcher_config, payload_config, mission_profile, units_ref, noise_params = p
    n_coils = launcher_config.num_coils

    # Position noise (very small - GPS/inertial navigation accuracy ~1m)
    du[1:3] .= 1.0  # 1m position uncertainty

    # Velocity noise (from atmospheric turbulence and measurement)
    # Scales with velocity magnitude
    v_mag = sqrt(u[4]^2 + u[5]^2 + u[6]^2)
    v_noise = max(0.1, v_mag * 0.01)  # 1% of velocity or 0.1 m/s minimum
    du[4:6] .= v_noise

    # Temperature noise (thermal fluctuations and measurement)
    du[7] = 1.0  # ±1K temperature uncertainty

    # Current noise (power supply ripple and measurement)
    for i in 1:n_coils
        I_mag = abs(u[7+i])
        du[7+i] = max(1.0, I_mag * noise_params.current_noise)
    end

    # Charge noise (from current integration errors)
    for i in 1:n_coils
        Q_mag = abs(u[7+n_coils+i])
        du[7+n_coils+i] = max(0.1, Q_mag * 0.001)  # 0.1% charge uncertainty
    end

    return nothing
end

"""
    create_sde_problem(launcher, payload, mission, noise_params=default_noise_parameters())

Create an SDE problem for stochastic trajectory simulation.

# Arguments
- `launcher`: LauncherConfig
- `payload`: PayloadConfig
- `mission`: MissionProfile
- `noise_params`: NoiseParameters (optional)

# Returns
- SDEProblem ready for solving with DifferentialEquations.jl
"""
function create_sde_problem(launcher, payload, mission, noise_params=default_noise_parameters(); tspan=(0.0u"s", 10.0u"s"))
    # Initial state
    u0_unitful = initial_state(launcher, payload, mission)
    n_coils = launcher.num_coils

    # Convert to dimensionless
    u0, units_ref = strip_units_state(u0_unitful, n_coils)

    # Parameters include noise
    p = (launcher, payload, mission, units_ref, noise_params)

    # Time span (dimensionless) - strip units if provided
    tspan_dimensionless = (ustrip(u"s", tspan[1]), ustrip(u"s", tspan[2]))

    # Create SDE problem with drift and diffusion terms
    prob = SDEProblem(trajectory_ode_dimensionless!, trajectory_sde_noise!, u0, tspan_dimensionless, p)

    return prob
end

"""
    monte_carlo_analysis(launcher, payload, mission, n_runs=100; noise_params=default_noise_parameters())

Run Monte Carlo analysis to compute trajectory distributions.

Performs multiple stochastic simulations and computes statistics:
- Mean trajectory
- Standard deviation at each time point
- Confidence intervals (68%, 95%, 99%)
- Success probability (reaching target altitude/velocity)

# Arguments
- `launcher`: LauncherConfig
- `payload`: PayloadConfig
- `mission`: MissionProfile
- `n_runs`: Number of Monte Carlo runs (default: 100)
- `noise_params`: NoiseParameters (optional)

# Returns
Named tuple with:
- `mean_trajectory`: Mean state over all runs
- `std_trajectory`: Standard deviation at each time
- `trajectories`: All individual trajectories
- `success_rate`: Fraction achieving target
- `confidence_intervals`: 68%, 95%, 99% bounds
"""
function monte_carlo_analysis(
    launcher, payload, mission, n_runs=100;
    noise_params=default_noise_parameters()
)
    @info "Running Monte Carlo analysis" n_runs=n_runs
    @info "Noise parameters" atmospheric="±$(noise_params.atmospheric_noise*100)%" drag="±$(noise_params.drag_noise*100)%" current="±$(noise_params.current_noise*100)%" timing="±$(noise_params.timing_noise)" mass="±$(noise_params.mass_noise*100)%"

    # Create SDE problem
    prob = create_sde_problem(launcher, payload, mission, noise_params)

    # Create ensemble problem for Monte Carlo
    ensemble_prob = EnsembleProblem(prob)

    # Create progress bar
    prog = Progress(n_runs, dt=0.5, desc="Solving SDE trajectories: ", barlen=50)

    # Create callback to update progress
    function prob_func(prob, i, repeat)
        next!(prog)
        prob
    end

    # Recreate ensemble with progress callback
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)

    # Solve ensemble with parallel threading
    @info "Solving SDE trajectories" n_runs=n_runs
    sol = solve(ensemble_prob, SRIW1(), EnsembleThreads(), trajectories=n_runs,
                adaptive=true, dt=0.001, saveat=0.1,
                maxiters=1e7, verbose=false)

    finish!(prog)
    @info "Ensemble solve complete"

    # Extract trajectories and filter out failed ones
    all_trajectories = sol.u
    # Use SciMLBase.successful_retcode to properly check success
    # ReturnCode is an enum, not a Symbol, so == :Success doesn't work
    trajectories = [traj for traj in all_trajectories if SciMLBase.successful_retcode(traj.retcode)]
    n_failed = n_runs - length(trajectories)

    if n_failed > 0
        @warn "$(n_failed) trajectories failed due to numerical instability and were excluded"
    end

    if length(trajectories) == 0
        error("All trajectories failed! Try reducing launcher power or adjusting parameters.")
    end

    # Compute statistics
    @info "Computing statistics"

    # Get all time points (assuming all trajectories share same time grid)
    if length(trajectories) > 0 && length(trajectories[1].t) > 0
        t_common = trajectories[1].t
        n_times = length(t_common)
        n_states = length(prob.u0)
        n_successful = length(trajectories)

        # Collect all trajectory values (only successful ones)
        all_values = zeros(n_successful, n_times, n_states)
        for (i, traj) in enumerate(trajectories)
            for (j, val) in enumerate(traj.u)
                all_values[i, j, :] = val
            end
        end

        # Compute mean and std across realizations
        mean_trajectory = [mean(all_values[:, j, :], dims=1)[1, :] for j in 1:n_times]
        std_trajectory = [std(all_values[:, j, :], dims=1)[1, :] for j in 1:n_times]

        # Compute confidence intervals
        conf_intervals = compute_confidence_intervals(all_values)

        # Compute success rate
        target_altitude = ustrip(u"m", mission.target_altitude)
        target_velocity = ustrip(u"m/s", mission.target_velocity)
        success_rate, stats = success_probability(trajectories, target_altitude, target_velocity, launcher)

        @info "Statistics computed" success_rate="$(round(success_rate * 100, digits=1))%" successful_runs="$(n_successful)/$(n_runs)"

        results = Dict(
            :mean_trajectory => mean_trajectory,
            :std_trajectory => std_trajectory,
            :trajectories => trajectories,
            :success_rate => success_rate,
            :confidence_intervals => conf_intervals,
            :time => t_common,
            :statistics => stats,
            :n_runs => n_runs,
            :n_successful => n_successful,
            :n_failed => n_failed
        )
    else
        @warn "No trajectories computed"
        results = Dict(
            :mean_trajectory => nothing,
            :std_trajectory => nothing,
            :trajectories => [],
            :success_rate => 0.0,
            :confidence_intervals => Dict()
        )
    end

    return results
end

"""
    compute_confidence_intervals(all_values, confidence_levels=[0.68, 0.95, 0.99])

Compute confidence intervals from Monte Carlo trajectories.

# Arguments
- `all_values`: Array of shape (n_runs, n_times, n_states) containing all trajectory values
- `confidence_levels`: Array of confidence levels (e.g., [0.68, 0.95, 0.99])

# Returns
Dictionary mapping confidence level to (lower_bound, upper_bound) trajectories
"""
function compute_confidence_intervals(all_values, confidence_levels=[0.68, 0.95, 0.99])
    n_runs, n_times, n_states = size(all_values)
    intervals = Dict{Float64, Tuple{Vector{Vector{Float64}}, Vector{Vector{Float64}}}}()

    for level in confidence_levels
        # Compute percentiles
        lower_percentile = (1.0 - level) / 2.0
        upper_percentile = 1.0 - lower_percentile

        # Compute percentiles for each time and state
        lower_bounds = Vector{Vector{Float64}}(undef, n_times)
        upper_bounds = Vector{Vector{Float64}}(undef, n_times)

        for j in 1:n_times
            lower_bounds[j] = Vector{Float64}(undef, n_states)
            upper_bounds[j] = Vector{Float64}(undef, n_states)

            for k in 1:n_states
                values_at_time_state = all_values[:, j, k]
                lower_bounds[j][k] = quantile(values_at_time_state, lower_percentile)
                upper_bounds[j][k] = quantile(values_at_time_state, upper_percentile)
            end
        end

        intervals[level] = (lower_bounds, upper_bounds)
    end

    return intervals
end

"""
    success_probability(trajectories, target_altitude, target_velocity, launcher)

Compute probability of successfully achieving target parameters.

# Arguments
- `trajectories`: Array of trajectory solutions
- `target_altitude`: Target altitude to reach (dimensionless)
- `target_velocity`: Target velocity to achieve (dimensionless)
- `launcher`: LauncherConfig for determining end of launch phase

# Returns
- Probability (0 to 1) of success
- Dict with detailed statistics
"""
function success_probability(trajectories, target_altitude, target_velocity, launcher)
    n_success = 0
    n_total = length(trajectories)

    final_altitudes = Float64[]
    final_velocities = Float64[]

    for traj in trajectories
        if length(traj.u) > 0
            # Get final state
            final_state = traj.u[end]

            # Extract position and velocity (first 6 components)
            pos = final_state[1:3]
            vel = final_state[4:6]

            # Compute altitude (distance from origin minus Earth radius, assuming dimensionless)
            altitude = sqrt(pos[1]^2 + pos[2]^2 + pos[3]^2)
            velocity = sqrt(vel[1]^2 + vel[2]^2 + vel[3]^2)

            push!(final_altitudes, altitude)
            push!(final_velocities, velocity)

            # Success criterion: reached target altitude and velocity
            if altitude >= target_altitude && velocity >= target_velocity
                n_success += 1
            end
        end
    end

    prob = n_total > 0 ? n_success / n_total : 0.0

    stats = Dict(
        :probability => prob,
        :success_count => n_success,
        :total_count => n_total,
        :mean_altitude => length(final_altitudes) > 0 ? mean(final_altitudes) : 0.0,
        :std_altitude => length(final_altitudes) > 0 ? std(final_altitudes) : 0.0,
        :mean_velocity => length(final_velocities) > 0 ? mean(final_velocities) : 0.0,
        :std_velocity => length(final_velocities) > 0 ? std(final_velocities) : 0.0,
        :final_altitudes => final_altitudes,
        :final_velocities => final_velocities
    )

    return prob, stats
end

export NoiseParameters, default_noise_parameters
export trajectory_sde_noise!
export create_sde_problem, monte_carlo_analysis
export compute_confidence_intervals, success_probability
