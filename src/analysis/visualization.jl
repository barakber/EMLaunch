"""
# Visualization Module

Plotting and visualization functions for trajectory and performance analysis.
Supports multiple celestial bodies via the `body` keyword argument.

Provides plots for:
- 3D trajectory
- Altitude vs time
- Velocity vs time
- Mach number profile
- Temperature evolution
- G-loading
- Forces breakdown
"""

using Unitful
using Plots
using Printf

"""
    plot_trajectory_3d(sol, n_coils; body=EARTH)

Create 3D plot of launch trajectory.

# Arguments
- `sol`: ODE solution
- `n_coils`: Number of coils
- `body`: CelestialBody (default: EARTH)

# Returns
- Plots object
"""
function plot_trajectory_3d(sol, n_coils; body=EARTH)
    traj = extract_trajectory_data(sol, n_coils; body=body)

    x = [ustrip(u"km", pos[1]) for pos in traj.positions]
    y = [ustrip(u"km", pos[2]) for pos in traj.positions]
    z = [ustrip(u"km", pos[3]) for pos in traj.positions]

    p = plot(x, y, z,
        xlabel = "X [km]",
        ylabel = "Y [km]",
        zlabel = "Z [km]",
        title = "Launch Trajectory (3D) - $(body.name)",
        label = "Trajectory",
        linewidth = 2,
        color = :blue,
        legend = :topright
    )

    scatter!(p, [x[1]], [y[1]], [z[1]],
        label = "Launch",
        markersize = 8,
        color = :green
    )

    scatter!(p, [x[end]], [y[end]], [z[end]],
        label = "Final",
        markersize = 8,
        color = :red
    )

    return p
end

"""
    plot_altitude_profile(sol, n_coils; body=EARTH)

Plot altitude vs time.

# Arguments
- `sol`: ODE solution
- `n_coils`: Number of coils
- `body`: CelestialBody (default: EARTH)

# Returns
- Plots object
"""
function plot_altitude_profile(sol, n_coils; body=EARTH)
    traj = extract_trajectory_data(sol, n_coils; body=body)

    times_s = ustrip.(u"s", traj.times)
    altitudes_km = ustrip.(u"km", traj.altitudes)

    p = plot(times_s, altitudes_km,
        xlabel = "Time [s]",
        ylabel = "Altitude [km]",
        title = "Altitude Profile - $(body.name)",
        label = "Altitude",
        linewidth = 2,
        color = :blue,
        legend = :best,
        grid = true
    )

    # Only show Karman line for bodies with atmosphere
    if has_atmosphere(body)
        hline!(p, [100.0], label="Karman line (100 km)", linestyle=:dash, color=:gray)
    end

    return p
end

"""
    plot_velocity_profile(sol, n_coils; body=EARTH)

Plot velocity magnitude vs time.

# Arguments
- `sol`: ODE solution
- `n_coils`: Number of coils
- `body`: CelestialBody (default: EARTH)

# Returns
- Plots object
"""
function plot_velocity_profile(sol, n_coils; body=EARTH)
    traj = extract_trajectory_data(sol, n_coils; body=body)

    times_s = ustrip.(u"s", traj.times)
    speeds_km_s = ustrip.(u"km/s", traj.speeds)

    p = plot(times_s, speeds_km_s,
        xlabel = "Time [s]",
        ylabel = "Velocity [km/s]",
        title = "Velocity Profile - $(body.name)",
        label = "Speed",
        linewidth = 2,
        color = :red,
        legend = :best,
        grid = true
    )

    # Compute orbital velocity for the body
    v_orb_km_s = ustrip(u"km/s", orbital_velocity(0.0u"km"; body=body))
    hline!(p, [v_orb_km_s], label="Orbital velocity ($(round(v_orb_km_s, digits=1)) km/s)", linestyle=:dash, color=:green)

    return p
end

"""
    plot_mach_profile(sol, n_coils; body=EARTH)

Plot Mach number vs time or altitude.

# Arguments
- `sol`: ODE solution
- `n_coils`: Number of coils
- `body`: CelestialBody (default: EARTH)

# Returns
- Plots object
"""
function plot_mach_profile(sol, n_coils; body=EARTH)
    traj = extract_trajectory_data(sol, n_coils; body=body)

    times_s = ustrip.(u"s", traj.times)
    mach = traj.mach_numbers

    p = plot(times_s, mach,
        xlabel = "Time [s]",
        ylabel = "Mach Number",
        title = "Mach Number Profile - $(body.name)",
        label = "Mach",
        linewidth = 2,
        color = :purple,
        legend = :best,
        grid = true
    )

    if has_atmosphere(body)
        hline!(p, [1.0], label="Mach 1 (Sonic)", linestyle=:dash, color=:gray, alpha=0.5)
        hline!(p, [5.0], label="Mach 5 (Hypersonic)", linestyle=:dash, color=:orange, alpha=0.5)
    end

    return p
end

"""
    plot_temperature_profile(sol, n_coils; body=EARTH)

Plot payload temperature vs time.

# Arguments
- `sol`: ODE solution
- `n_coils`: Number of coils
- `body`: CelestialBody (default: EARTH)

# Returns
- Plots object
"""
function plot_temperature_profile(sol, n_coils; body=EARTH)
    traj = extract_trajectory_data(sol, n_coils; body=body)

    times_s = ustrip.(u"s", traj.times)
    temps_K = ustrip.(u"K", traj.temperatures)

    p = plot(times_s, temps_K,
        xlabel = "Time [s]",
        ylabel = "Temperature [K]",
        title = "Payload Temperature - $(body.name)",
        label = "Temperature",
        linewidth = 2,
        color = :orange,
        legend = :best,
        grid = true
    )

    hline!(p, [1500.0], label="Thermal limit (1500 K)", linestyle=:dash, color=:red)

    return p
end

"""
    plot_mission_overview(sol, n_coils; save_path=nothing, body=EARTH)

Create comprehensive multi-panel overview plot.

# Arguments
- `sol`: ODE solution
- `n_coils`: Number of coils
- `save_path`: Optional path to save figure
- `body`: CelestialBody (default: EARTH)

# Returns
- Plots object with multiple subplots
"""
function plot_mission_overview(sol, n_coils; save_path=nothing, body=EARTH)
    p1 = plot_altitude_profile(sol, n_coils; body=body)
    p2 = plot_velocity_profile(sol, n_coils; body=body)
    p3 = plot_mach_profile(sol, n_coils; body=body)
    p4 = plot_temperature_profile(sol, n_coils; body=body)

    p = plot(p1, p2, p3, p4,
        layout = (2, 2),
        size = (1200, 900),
        plot_title = "Mission Overview - $(body.name)",
        margin = 5Plots.mm
    )

    if !isnothing(save_path)
        savefig(p, save_path)
        println("Saved plot to: $save_path")
    end

    return p
end

"""
    plot_phase_diagram(sol, n_coils; body=EARTH)

Create phase diagram (altitude vs velocity).

# Arguments
- `sol`: ODE solution
- `n_coils`: Number of coils
- `body`: CelestialBody (default: EARTH)

# Returns
- Plots object
"""
function plot_phase_diagram(sol, n_coils; body=EARTH)
    traj = extract_trajectory_data(sol, n_coils; body=body)

    altitudes_km = ustrip.(u"km", traj.altitudes)
    speeds_km_s = ustrip.(u"km/s", traj.speeds)

    p = plot(speeds_km_s, altitudes_km,
        xlabel = "Velocity [km/s]",
        ylabel = "Altitude [km]",
        title = "Phase Diagram - $(body.name)",
        label = "Trajectory",
        linewidth = 2,
        color = :blue,
        legend = :best,
        grid = true
    )

    scatter!(p, [speeds_km_s[1]], [altitudes_km[1]],
        label = "Launch",
        markersize = 8,
        color = :green
    )

    scatter!(p, [speeds_km_s[end]], [altitudes_km[end]],
        label = "Final",
        markersize = 8,
        color = :red
    )

    return p
end

"""
    plot_energy_diagram(sol, n_coils, payload_mass; body=EARTH)

Plot kinetic and potential energy evolution.

# Arguments
- `sol`: ODE solution
- `n_coils`: Number of coils
- `payload_mass`: Payload mass [kg]
- `body`: CelestialBody (default: EARTH)

# Returns
- Plots object
"""
function plot_energy_diagram(sol, n_coils, payload_mass; body=EARTH)
    traj = extract_trajectory_data(sol, n_coils; body=body)

    times_s = ustrip.(u"s", traj.times)

    μ = ustrip(u"m^3/s^2", body.mu)
    R = ustrip(u"m", body.radius)

    KE = [0.5 * ustrip(u"kg", payload_mass) * ustrip(u"m/s", v)^2 / 1e6
          for v in traj.speeds]  # MJ

    PE = [-ustrip(u"kg", payload_mass) * μ /
          (R + ustrip(u"m", h)) / 1e6
          for h in traj.altitudes]  # MJ

    TE = KE .+ PE

    p = plot(times_s, [KE PE TE],
        xlabel = "Time [s]",
        ylabel = "Energy [MJ]",
        title = "Energy Evolution - $(body.name)",
        label = ["Kinetic" "Potential" "Total"],
        linewidth = 2,
        legend = :best,
        grid = true
    )

    return p
end

"""
    animate_trajectory(sol, n_coils; fps=30, duration=10, body=EARTH)

Create animation of trajectory (requires Plots animation support).

# Arguments
- `sol`: ODE solution
- `n_coils`: Number of coils
- `fps`: Frames per second
- `duration`: Animation duration [s]
- `body`: CelestialBody (default: EARTH)

# Returns
- Animation object
"""
function animate_trajectory(sol, n_coils; fps=30, duration=10, body=EARTH)
    traj = extract_trajectory_data(sol, n_coils; body=body)

    n_frames = fps * duration
    indices = round.(Int, range(1, length(traj.times), length=n_frames))

    anim = @animate for i in indices
        alt = ustrip(u"km", traj.altitudes[i])
        vel = ustrip(u"km/s", traj.speeds[i])
        t = ustrip(u"s", traj.times[i])

        altitudes_km = ustrip.(u"km", traj.altitudes[1:i])
        speeds_km_s = ustrip.(u"km/s", traj.speeds[1:i])

        plot(speeds_km_s, altitudes_km,
            xlabel = "Velocity [km/s]",
            ylabel = "Altitude [km]",
            title = @sprintf("t = %.1f s | v = %.2f km/s | h = %.1f km", t, vel, alt),
            label = "Trajectory",
            linewidth = 2,
            color = :blue,
            xlims = (0, maximum(ustrip.(u"km/s", traj.speeds)) * 1.1),
            ylims = (0, maximum(ustrip.(u"km", traj.altitudes)) * 1.1),
            legend = :bottomright,
            grid = true
        )

        scatter!([vel], [alt],
            label = "Current",
            markersize = 10,
            color = :red
        )
    end

    return anim
end

# ============================================================================
# Monte Carlo Visualization Functions (PlotlyJS)
# ============================================================================

using PlotlyJS

"""
    extract_altitude(state::Vector; body=EARTH)

Extract altitude (in km) from a state vector.

# Arguments
- `state`: State vector [x, y, z, vx, vy, vz, T]
- `body`: CelestialBody (default: EARTH)

# Returns
- Altitude in kilometers above body's surface
"""
function extract_altitude(state::Vector; body=EARTH)
    pos = [state[1], state[2], state[3]]
    R = ustrip(u"m", body.radius)
    return (sqrt(sum(pos.^2)) - R) / 1000  # km
end

"""
    extract_velocity(state::Vector)

Extract velocity magnitude (in m/s) from a state vector.

# Arguments
- `state`: State vector [x, y, z, vx, vy, vz, T]

# Returns
- Velocity magnitude in m/s
"""
function extract_velocity(state::Vector)
    vel = [state[4], state[5], state[6]]
    return sqrt(sum(vel.^2))  # m/s
end

"""
    extract_temperature(state::Vector)

Extract temperature (in K) from a state vector.

# Arguments
- `state`: State vector [x, y, z, vx, vy, vz, T]

# Returns
- Temperature in Kelvin
"""
function extract_temperature(state::Vector)
    return state[7]  # K
end

"""
    compute_percentiles(trajectories, time_points, extractor_fn, percentiles=[0.5, 0.05, 0.95])

Compute specified percentiles from Monte Carlo trajectories.

# Arguments
- `trajectories`: Vector of trajectory solutions
- `time_points`: Time points for the trajectories
- `extractor_fn`: Function to extract desired quantity from state vector
- `percentiles`: Vector of percentiles to compute (e.g., [0.5, 0.25, 0.75] for median, 25th, 75th)

# Returns
- Tuple of vectors, one for each requested percentile in the same order
"""
function compute_percentiles(trajectories, time_points, extractor_fn, percentiles=[0.5, 0.05, 0.95])
    n_times = length(time_points)
    n_percentiles = length(percentiles)

    result_vals = [zeros(n_times) for _ in 1:n_percentiles]

    for t_idx in 1:n_times
        values = Float64[]
        for traj in trajectories
            if length(traj.u) >= t_idx
                push!(values, extractor_fn(traj.u[t_idx]))
            end
        end

        if !isempty(values)
            sorted_vals = sort(values)
            n = length(sorted_vals)

            for (i, p) in enumerate(percentiles)
                idx = max(1, min(n, round(Int, n * p)))
                result_vals[i][t_idx] = sorted_vals[idx]
            end
        end
    end

    return tuple(result_vals...)
end

"""
    create_trajectory_dashboard(results; height=900, body=EARTH)

Create an interactive dashboard with altitude, velocity, and temperature plots.

# Arguments
- `results`: Results dictionary from `monte_carlo_analysis`
- `height`: Plot height in pixels (default: 900)
- `body`: CelestialBody (default: EARTH)

# Returns
- PlotlyJS plot object with three subplots showing trajectory statistics
"""
function create_trajectory_dashboard(results; height=900, body=EARTH)
    time_points = results[:time]
    trajectories = results[:trajectories]

    # Create body-aware altitude extractor
    alt_fn = state -> extract_altitude(state; body=body)

    alt_median, alt_p5, alt_p95 = compute_percentiles(trajectories, time_points, alt_fn)
    vel_median, vel_p5, vel_p95 = compute_percentiles(trajectories, time_points, extract_velocity)
    temp_median, temp_p5, temp_p95 = compute_percentiles(trajectories, time_points, extract_temperature)

    # Altitude traces
    trace_alt_band = PlotlyJS.scatter(
        x=vcat(time_points, reverse(time_points)),
        y=vcat(alt_p95, reverse(alt_p5)),
        fill="toself",
        fillcolor="rgba(0,100,255,0.2)",
        line=attr(color="rgba(255,255,255,0)"),
        name="90% CI",
        showlegend=true,
        hoverinfo="skip"
    )

    trace_alt_median = PlotlyJS.scatter(
        x=time_points,
        y=alt_median,
        mode="lines",
        name="Median Altitude",
        line=attr(color="rgb(0,100,255)", width=3),
        hovertemplate="Time: %{x:.2f}s<br>Altitude: %{y:.2f} km<extra></extra>"
    )

    # Conditional Karman line (only for atmospheric bodies)
    alt_traces = [trace_alt_band, trace_alt_median]
    if has_atmosphere(body)
        trace_alt_karman = PlotlyJS.scatter(
            x=time_points,
            y=fill(100.0, length(time_points)),
            mode="lines",
            name="Karman Line (100 km)",
            line=attr(color="red", width=2, dash="dash"),
            hoverinfo="skip"
        )
        push!(alt_traces, trace_alt_karman)
    end

    # Velocity traces
    trace_vel_band = PlotlyJS.scatter(
        x=vcat(time_points, reverse(time_points)),
        y=vcat(vel_p95, reverse(vel_p5)),
        fill="toself",
        fillcolor="rgba(0,200,100,0.2)",
        line=attr(color="rgba(255,255,255,0)"),
        name="90% CI",
        showlegend=true,
        hoverinfo="skip"
    )

    trace_vel_median = PlotlyJS.scatter(
        x=time_points,
        y=vel_median,
        mode="lines",
        name="Median Velocity",
        line=attr(color="rgb(0,200,100)", width=3),
        hovertemplate="Time: %{x:.2f}s<br>Velocity: %{y:.2f} m/s<br>Mach: %{customdata:.2f}<extra></extra>",
        customdata=vel_median ./ 340.0
    )

    v_orb = ustrip(u"m/s", orbital_velocity(0.0u"km"; body=body))
    trace_vel_orbital = PlotlyJS.scatter(
        x=time_points,
        y=fill(v_orb, length(time_points)),
        mode="lines",
        name="Orbital Velocity ($(round(Int, v_orb)) m/s)",
        line=attr(color="red", width=2, dash="dash"),
        hoverinfo="skip"
    )

    # Temperature traces
    trace_temp_band = PlotlyJS.scatter(
        x=vcat(time_points, reverse(time_points)),
        y=vcat(temp_p95, reverse(temp_p5)),
        fill="toself",
        fillcolor="rgba(255,140,0,0.2)",
        line=attr(color="rgba(255,255,255,0)"),
        name="90% CI",
        showlegend=true,
        hoverinfo="skip"
    )

    trace_temp_median = PlotlyJS.scatter(
        x=time_points,
        y=temp_median,
        mode="lines",
        name="Median Temperature",
        line=attr(color="rgb(255,140,0)", width=3),
        hovertemplate="Time: %{x:.2f}s<br>Temperature: %{y:.2f} K<extra></extra>"
    )

    trace_temp_max = PlotlyJS.scatter(
        x=time_points,
        y=fill(1800.0, length(time_points)),
        mode="lines",
        name="Max Temp (1800 K)",
        line=attr(color="red", width=2, dash="dash"),
        hoverinfo="skip"
    )

    # Create subplots
    p = PlotlyJS.make_subplots(
        rows=3, cols=1,
        subplot_titles=["Altitude vs Time" "Velocity vs Time" "Temperature vs Time"],
        vertical_spacing=0.1
    )

    for trace in alt_traces
        PlotlyJS.add_trace!(p, trace, row=1, col=1)
    end

    PlotlyJS.add_trace!(p, trace_vel_band, row=2, col=1)
    PlotlyJS.add_trace!(p, trace_vel_median, row=2, col=1)
    PlotlyJS.add_trace!(p, trace_vel_orbital, row=2, col=1)

    PlotlyJS.add_trace!(p, trace_temp_band, row=3, col=1)
    PlotlyJS.add_trace!(p, trace_temp_median, row=3, col=1)
    PlotlyJS.add_trace!(p, trace_temp_max, row=3, col=1)

    PlotlyJS.relayout!(p,
        height=height,
        title_text="EMLaunch Trajectory Analysis - $(body.name)",
        showlegend=true,
        hovermode="x unified",
        template="plotly_white"
    )

    PlotlyJS.relayout!(p,
        xaxis_title="Time (s)",
        yaxis_title="Altitude (km)",
        xaxis2_title="Time (s)",
        yaxis2_title="Velocity (m/s)",
        xaxis3_title="Time (s)",
        yaxis3_title="Temperature (K)"
    )

    return p
end

"""
    plot_monte_carlo_results(results; size=(1200, 900), body=EARTH)

Create static plots for Monte Carlo results using Plots.jl.

# Arguments
- `results`: Results dictionary from `monte_carlo_analysis`
- `size`: Plot size (width, height) in pixels
- `body`: CelestialBody (default: EARTH)

# Returns
- Plots object with three subplots showing trajectory statistics
"""
function plot_monte_carlo_results(results; size=(1200, 900), body=EARTH)
    time_points = results[:time]
    trajectories = results[:trajectories]

    alt_fn = state -> extract_altitude(state; body=body)

    # 50% CI
    alt_median, alt_p25, alt_p75 = compute_percentiles(trajectories, time_points, alt_fn, [0.5, 0.25, 0.75])
    vel_median, vel_p25, vel_p75 = compute_percentiles(trajectories, time_points, extract_velocity, [0.5, 0.25, 0.75])
    temp_median, temp_p25, temp_p75 = compute_percentiles(trajectories, time_points, extract_temperature, [0.5, 0.25, 0.75])

    # 90% CI
    _, alt_p5, alt_p95 = compute_percentiles(trajectories, time_points, alt_fn, [0.5, 0.05, 0.95])
    _, vel_p5, vel_p95 = compute_percentiles(trajectories, time_points, extract_velocity, [0.5, 0.05, 0.95])
    _, temp_p5, temp_p95 = compute_percentiles(trajectories, time_points, extract_temperature, [0.5, 0.05, 0.95])

    # 98% CI
    _, alt_p1, alt_p99 = compute_percentiles(trajectories, time_points, alt_fn, [0.5, 0.01, 0.99])
    _, vel_p1, vel_p99 = compute_percentiles(trajectories, time_points, extract_velocity, [0.5, 0.01, 0.99])
    _, temp_p1, temp_p99 = compute_percentiles(trajectories, time_points, extract_temperature, [0.5, 0.01, 0.99])

    # Altitude plot
    p1 = Plots.plot(legend=:topleft, grid=true, xlabel="Time [s]", ylabel="Altitude [km]", title="Altitude vs Time - $(body.name)")

    Plots.plot!(p1, time_points, alt_median,
        ribbon=(alt_median .- alt_p1, alt_p99 .- alt_median),
        fillalpha=0.15, linewidth=0, color=:blue, label="98% CI")

    Plots.plot!(p1, time_points, alt_median,
        ribbon=(alt_median .- alt_p5, alt_p95 .- alt_median),
        fillalpha=0.25, linewidth=0, color=:blue, label="90% CI")

    Plots.plot!(p1, time_points, alt_median,
        ribbon=(alt_median .- alt_p25, alt_p75 .- alt_median),
        fillalpha=0.40, linewidth=0, color=:blue, label="50% CI")

    Plots.plot!(p1, time_points, alt_median,
        linewidth=2.5, color=:blue, label="Median")

    if has_atmosphere(body)
        Plots.hline!(p1, [100.0], label="Karman line", linestyle=:dash, color=:red, linewidth=2)
    end

    # Velocity plot
    p2 = Plots.plot(legend=:topright, grid=true, xlabel="Time [s]", ylabel="Velocity [m/s]", title="Velocity vs Time - $(body.name)")

    Plots.plot!(p2, time_points, vel_median,
        ribbon=(vel_median .- vel_p1, vel_p99 .- vel_median),
        fillalpha=0.15, linewidth=0, color=:green, label="98% CI")

    Plots.plot!(p2, time_points, vel_median,
        ribbon=(vel_median .- vel_p5, vel_p95 .- vel_median),
        fillalpha=0.25, linewidth=0, color=:green, label="90% CI")

    Plots.plot!(p2, time_points, vel_median,
        ribbon=(vel_median .- vel_p25, vel_p75 .- vel_median),
        fillalpha=0.40, linewidth=0, color=:green, label="50% CI")

    Plots.plot!(p2, time_points, vel_median,
        linewidth=2.5, color=:green, label="Median")

    v_orb = ustrip(u"m/s", orbital_velocity(0.0u"km"; body=body))
    Plots.hline!(p2, [v_orb], label="Orbital velocity", linestyle=:dash, color=:red, linewidth=2)

    # Temperature plot
    p3 = Plots.plot(legend=:topleft, grid=true, xlabel="Time [s]", ylabel="Temperature [K]", title="Temperature vs Time - $(body.name)")

    Plots.plot!(p3, time_points, temp_median,
        ribbon=(temp_median .- temp_p1, temp_p99 .- temp_median),
        fillalpha=0.15, linewidth=0, color=:orange, label="98% CI")

    Plots.plot!(p3, time_points, temp_median,
        ribbon=(temp_median .- temp_p5, temp_p95 .- temp_median),
        fillalpha=0.25, linewidth=0, color=:orange, label="90% CI")

    Plots.plot!(p3, time_points, temp_median,
        ribbon=(temp_median .- temp_p25, temp_p75 .- temp_median),
        fillalpha=0.40, linewidth=0, color=:orange, label="50% CI")

    Plots.plot!(p3, time_points, temp_median,
        linewidth=2.5, color=:orange, label="Median")

    Plots.hline!(p3, [2500.0], label="Max temperature", linestyle=:dash, color=:red, linewidth=2)

    p = Plots.plot(p1, p2, p3,
        layout=(3, 1),
        size=size,
        plot_title="Monte Carlo Trajectory Analysis - $(body.name)",
        margin=5Plots.mm
    )

    return p
end

"""
    print_summary_statistics(results; body=EARTH)

Print a summary of Monte Carlo simulation results.

# Arguments
- `results`: Results dictionary from `monte_carlo_analysis`
- `body`: CelestialBody (default: EARTH)
"""
function print_summary_statistics(results; body=EARTH)
    mean_traj = results[:mean_trajectory]
    success_rate = results[:success_rate]

    final_state = mean_traj[end]

    idx_x, idx_y, idx_z = 1, 2, 3
    idx_vx, idx_vy, idx_vz = 4, 5, 6
    idx_T = 7

    position = [final_state[idx_x], final_state[idx_y], final_state[idx_z]]
    velocity = [final_state[idx_vx], final_state[idx_vy], final_state[idx_vz]]
    R = ustrip(u"m", body.radius)
    altitude = sqrt(sum(position.^2)) - R
    speed = sqrt(sum(velocity.^2))
    temperature = final_state[idx_T]

    println("\n" * "="^70)
    println("MONTE CARLO SIMULATION RESULTS ($(body.name))")
    println("="^70)

    n_runs = get(results, :n_runs, 0)
    n_successful = get(results, :n_successful, 0)
    n_failed = get(results, :n_failed, 0)

    if n_runs > 0
        println("\nSimulation Statistics:")
        println("  Total runs:       $(n_runs)")
        println("  Successful runs:  $(n_successful)")
        if n_failed > 0
            println("  Failed runs:      $(n_failed) (numerical instability)")
        end
    end

    println("\nFinal State (Mean Trajectory):")
    println("  Altitude:     $(round(altitude/1000, digits=2)) km")
    if has_atmosphere(body)
        println("  Speed:        $(round(speed, digits=2)) m/s (Mach $(round(speed/340, digits=2)))")
    else
        println("  Speed:        $(round(speed, digits=2)) m/s")
    end
    println("  Temperature:  $(round(temperature, digits=2)) K")
    println("\nMission Performance:")
    println("  Success rate: $(round(success_rate*100, digits=1))%")

    if has_atmosphere(body)
        if altitude > 100000
            println("\n  SUCCESS! Crossed the Karman line (100 km)")
            println("  Altitude above Karman line: $(round((altitude-100000)/1000, digits=2)) km")
        else
            println("\n  Did not reach space (100 km Karman line)")
            println("  Altitude deficit: $(round((100000-altitude)/1000, digits=2)) km")
        end
    else
        v_esc = ustrip(u"m/s", escape_velocity(0.0u"km"; body=body))
        if speed > v_esc
            println("\n  Achieved escape velocity!")
        end
    end
    println("="^70)
end

export plot_trajectory_3d, plot_altitude_profile, plot_velocity_profile
export plot_mach_profile, plot_temperature_profile
export plot_mission_overview, plot_phase_diagram, plot_energy_diagram
export animate_trajectory
export extract_altitude, extract_velocity, extract_temperature
export compute_percentiles, create_trajectory_dashboard, plot_monte_carlo_results, print_summary_statistics
