"""
# EMLaunch.jl

Electromagnetic launch system modeling package.
Models the physics of electromagnetic acceleration, atmospheric flight,
orbital insertion, and uncertainty quantification with type-safe units throughout.

## Modules
- `PhysicalParameters`: Physical constants and material properties
- `EMAcceleration`: Electromagnetic coilgun dynamics
- `Atmosphere`: Atmospheric models (drag, heating, density)
- `DragReduction`: Passive drag reduction methods (aerospikes, waveriders, etc.)
- `Thermal`: Ablative heat shields and thermal protection systems
- `Gravity`: Gravitational models including perturbations
- `Trajectory`: ODE system for full trajectory simulation
- `StochasticTrajectory`: SDE system for uncertainty quantification
- `Optimization`: Parameter optimization for mission profiles
- `Simulation`: High-level simulation interface

Author: Barak Bercovitz (barakber@gmail.com)
"""
module EMLaunch

using Unitful
using UnitfulAstro
using PhysicalConstants.CODATA2018: c_0, μ_0, ε_0, k_B, σ
using StaticArrays
using LinearAlgebra

# Export main types and functions
export LauncherConfig, PayloadConfig, MissionProfile
export simulate_launch, optimize_parameters, optimize_launch_angle_min_drag
export plot_trajectory, plot_forces, plot_thermal

# Export launcher/payload creation
export create_uniform_launcher, create_default_payload

# Export stochastic analysis functions
export NoiseParameters, default_noise_parameters
export trajectory_sde_noise!
export create_sde_problem, create_simplified_sde_problem, monte_carlo_analysis
export compute_confidence_intervals, success_probability

# Export thermal/ablation functions
export AblativeMaterial, HeatShield, ABLATIVE_MATERIALS
export ablation_rate, ablation_thickness_rate
export required_heat_shield_mass, required_heat_shield_thickness
export update_heat_shield!, heat_shield_intact, heat_shield_margin
export thermal_survival_check
export print_heat_shield_summary, print_ablative_materials

# Export drag reduction functions
export DragReductionMethod, NoDragReduction
export Aerospike, Aerodisk, BoatTailing, Waverider, PorousSurface
export CompositeDragReduction
export drag_reduction_factor, added_mass, effective_drag_coefficient
export print_drag_reduction_summary

# Include submodules
# Core physics modules
include("core/physical_parameters.jl")
include("core/celestial_body.jl")
include("core/atmosphere.jl")
include("core/gravity.jl")

# Propulsion systems
include("propulsion/em_acceleration.jl")

# Flight dynamics and thermal protection
include("dynamics/trajectory.jl")
include("dynamics/stochastic_trajectory.jl")
include("dynamics/drag_reduction.jl")
include("dynamics/thermal.jl")

# Analysis and optimization
include("analysis/optimization.jl")
include("analysis/simulation.jl")
include("analysis/visualization.jl")

end # module
