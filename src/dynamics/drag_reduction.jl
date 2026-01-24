"""
# Passive Drag Reduction Module

Models passive (non-powered) methods to reduce atmospheric drag on hypersonic projectiles.

These methods don't require onboard energy, real-time control, or active systems.
They rely on geometry, materials, and static features to manipulate airflow, shock waves,
and boundary layers.

## Passive Drag Reduction Methods

### 1. Geometric Methods
- **Aerospike**: Forward-facing spike creates detached shock wave (20-50% reduction)
- **Aerodisk**: Disk at nose creates multiple shock interactions (up to 30% reduction)
- **Boat Tailing**: Tapered rear reduces base drag (10-20% reduction)
- **Waverider**: Lifting body rides own shock wave (15-25% reduction)

### 2. Surface Methods
- **Porous Surface**: Micro-perforations allow passive transpiration (10-20% reduction)
- **Surface Features**: Strategic roughness/dimples modify boundary layer (5-15% reduction)

### 3. Material Methods
- **Ablative Coating**: Vaporization creates boundary layer (5-10% indirect reduction)

## Effectiveness vs Mach Number

Most methods are Mach-dependent:
- Peak effectiveness: Mach 5-15 (hypersonic regime)
- Reduced effectiveness: Mach 1-5 (transonic/supersonic)
- Minimal effect: Mach < 1 (subsonic)

## Limitations

- Typical reductions: 10-40% of base drag
- At extreme velocities (>7 km/s) near sea level, insufficient alone
- Best combined with high-altitude launch (e.g., StarTram at 20 km)
- Mass penalties: 5-15% of payload mass for hardware

## References

1. Kalimuthu, R. et al. (2008). "Experimental Investigation on Spiked Body in Hypersonic Flow." *Aeronautical Journal*.
2. Ahmed, M. Y. M. & Qin, N. (2010). "Drag Reduction Using Aerodisks for Hypersonic Hemispherical Bodies." *AIAA Journal of Spacecraft*.
3. Reding, J. P. & Jecmen, D. M. (1981). "Base Region Flow Field of Spiked Bodies." *AIAA Journal*.
4. Cui, K. et al. (2015). "Hypersonic Drag and Heat Reduction Using a Forward-Facing Jet." *Aerospace Science and Technology*.
"""

using Unitful
using Printf
using StaticArrays

# =============================================================================
# DRAG REDUCTION METHOD TYPES
# =============================================================================

"""
Base type for passive drag reduction methods.
"""
abstract type DragReductionMethod end

"""
No drag reduction (baseline).
"""
struct NoDragReduction <: DragReductionMethod
end

"""
    Aerospike(length_ratio, diameter, material_density)

Forward-facing spike creates detached shock wave ahead of the body.

# Physics
- Creates low-pressure recirculation zone
- Pushes bow shock away from body
- Reduces wave drag by 20-50% at Mach 5-10

# Parameters
- `length_ratio`: Spike length / capsule diameter (typically 1.0-2.0)
- `diameter`: Spike diameter [m]
- `material_density`: Material density [kg/m³] (e.g., 2700 for aluminum, 8000 for tungsten)

# Mass Penalty
Spike mass ≈ π/4 × diameter² × length × density
Typically adds 5-10% to payload mass

# Limitations
- Extreme heating on spike tip (requires UHTC: ZrB₂, HfC)
- Structural loads during acceleration
- Effectiveness peaks at Mach 5-15
"""
struct Aerospike <: DragReductionMethod
    length_ratio::Float64  # Length / capsule diameter
    diameter::typeof(1.0u"m")
    material_density::typeof(1.0u"kg/m^3")
end

"""
    Aerodisk(disk_diameter_ratio, num_disks, material_density)

Disk(s) at nose create multiple shock interactions.

# Physics
- Multiple shock reflections rarify airflow
- Reduces drag by up to 30% in hypersonic flow
- Multi-row disk (MRD) improves stability

# Parameters
- `disk_diameter_ratio`: Disk diameter / capsule diameter (typically 0.5-0.8)
- `num_disks`: Number of stacked disks (1-3)
- `material_density`: Material density [kg/m³]

# Mass Penalty
Lighter than aerospike (disk vs cone)
Typically adds 3-7% to payload mass
"""
struct Aerodisk <: DragReductionMethod
    disk_diameter_ratio::Float64
    num_disks::Int
    material_density::typeof(1.0u"kg/m^3")
end

"""
    BoatTailing(taper_angle, length_ratio)

Tapered rear section reduces base drag.

# Physics
- Narrows wake recirculation zone
- Base drag can be 40% of total drag
- Reduces base drag by 10-20% even in hypersonic

# Parameters
- `taper_angle`: Taper angle [degrees] (typically 5-10°)
- `length_ratio`: Taper length / capsule diameter (typically 0.5-1.0)

# Mass Penalty
Minimal (reshaping, no added mass)
May actually reduce mass vs cylindrical base
"""
struct BoatTailing <: DragReductionMethod
    taper_angle::typeof(1.0u"°")
    length_ratio::Float64
end

"""
    Waverider(lift_to_drag_ratio, slenderness)

Lifting body that "rides" its own shock wave.

# Physics
- Caret or conical waverider shape
- Channels airflow to minimize pressure drag
- Reduces drag coefficient by 15-25% vs blunt shapes

# Parameters
- `lift_to_drag_ratio`: L/D ratio (typically 3-10)
- `slenderness`: Length/diameter ratio (typically 10-20)

# Tradeoffs
- Lower drag but higher stagnation heating
- Requires precise aerodynamic design
- Not suitable for all payload geometries
"""
struct Waverider <: DragReductionMethod
    lift_to_drag_ratio::Float64
    slenderness::Float64
end

"""
    PorousSurface(porosity, pore_diameter)

Micro-perforated surface allows passive transpiration cooling.

# Physics
- Air seeps through, cooling boundary layer
- Delays shock-induced separation
- Reduces skin friction drag by 10-20%

# Parameters
- `porosity`: Fraction of surface area that is pores (0.0-1.0)
- `pore_diameter`: Average pore size [m] (typically 10-100 μm)

# Mass Penalty
Porous ceramics/sintered metals slightly denser
Adds 2-5% mass
"""
struct PorousSurface <: DragReductionMethod
    porosity::Float64  # 0.0 to 1.0
    pore_diameter::typeof(1.0u"m")
end

"""
    CompositeDragReduction(methods)

Combine multiple passive drag reduction methods.

# Example
```julia
composite = CompositeDragReduction([
    Aerospike(1.5, 0.05u"m", 8000.0u"kg/m^3"),  # Tungsten spike
    BoatTailing(7.0u"°", 0.8),
    PorousSurface(0.15, 50.0e-6u"m")
])
```

Reductions are compounded multiplicatively (conservative):
C_d_effective = C_d_base × (1 - r₁) × (1 - r₂) × ...
"""
struct CompositeDragReduction <: DragReductionMethod
    methods::Vector{DragReductionMethod}
end

# =============================================================================
# DRAG REDUCTION FACTOR CALCULATIONS
# =============================================================================

"""
    drag_reduction_factor(method, mach, altitude)

Calculate drag reduction factor (0.0 to 1.0) for a given method.

Returns fraction by which drag is reduced:
- 0.0 = no reduction
- 0.3 = 30% reduction
- 0.5 = 50% reduction

# Arguments
- `method`: Drag reduction method
- `mach`: Mach number (dimensionless, can be unitful or unitless)
- `altitude`: Altitude (can be [m] with units or unitless)

# Returns
- Reduction factor (0.0 to 1.0, unitless)

# Note
This function handles both unitful and unitless inputs.
For ODE integration (performance), pass unitless values.
For analysis, you can pass units which will be stripped internally.
"""
function drag_reduction_factor(method::NoDragReduction, mach, altitude)
    return 0.0  # No reduction
end

function drag_reduction_factor(spike::Aerospike, mach, altitude)
    # Aerospike effectiveness vs Mach number
    # Peak at Mach 5-15, less effective outside this range

    # Base reduction: 20-50% depending on spike length
    base_reduction = 0.20 + 0.30 * min(1.0, spike.length_ratio / 2.0)

    # Mach effectiveness curve (empirical from literature)
    if mach < 1.0
        mach_factor = 0.1  # Minimal effect subsonic
    elseif mach < 3.0
        mach_factor = 0.3 + 0.4 * (mach - 1.0) / 2.0  # Ramp up
    elseif mach < 15.0
        mach_factor = 1.0  # Peak effectiveness
    else
        mach_factor = 1.0 - 0.05 * (mach - 15.0)  # Gradual decline
        mach_factor = max(0.5, mach_factor)
    end

    return base_reduction * mach_factor
end

function drag_reduction_factor(disk::Aerodisk, mach, altitude)
    # Aerodisk effectiveness
    # Up to 30% reduction at hypersonic speeds

    # Base reduction: scales with disk size and number
    base_reduction = 0.15 + 0.10 * (disk.disk_diameter_ratio - 0.5) / 0.3
    base_reduction += 0.05 * (disk.num_disks - 1)  # MRD bonus
    base_reduction = min(0.30, base_reduction)

    # Mach effectiveness
    if mach < 2.0
        mach_factor = 0.2
    elseif mach < 10.0
        mach_factor = 0.5 + 0.5 * (mach - 2.0) / 8.0
    else
        mach_factor = 1.0
    end

    return base_reduction * mach_factor
end

function drag_reduction_factor(tailing::BoatTailing, mach, altitude)
    # Boat tailing reduces base drag
    # 10-20% reduction across speed ranges

    # Base reduction: scales with taper angle
    angle_deg = ustrip(u"°", tailing.taper_angle)
    base_reduction = 0.10 + 0.10 * min(1.0, angle_deg / 10.0)

    # Effective across all Mach (primarily base drag)
    # Slightly less effective at extreme hypersonic
    if mach > 20.0
        mach_factor = 0.8
    else
        mach_factor = 1.0
    end

    return base_reduction * mach_factor
end

function drag_reduction_factor(rider::Waverider, mach, altitude)
    # Waverider effectiveness
    # 15-25% reduction in drag coefficient

    # Base reduction: scales with L/D and slenderness
    base_reduction = 0.15 + 0.10 * min(1.0, rider.lift_to_drag_ratio / 10.0)

    # Most effective in hypersonic regime
    if mach < 3.0
        mach_factor = 0.3
    elseif mach < 8.0
        mach_factor = 0.5 + 0.5 * (mach - 3.0) / 5.0
    else
        mach_factor = 1.0
    end

    return base_reduction * mach_factor
end

function drag_reduction_factor(porous::PorousSurface, mach, altitude)
    # Porous surface transpiration cooling
    # 10-20% skin friction reduction

    # Base reduction: scales with porosity
    base_reduction = 0.10 + 0.10 * porous.porosity / 0.2  # Assumes ~20% optimal
    base_reduction = min(0.20, base_reduction)

    # Effective across Mach, but best in hypersonic
    if mach < 2.0
        mach_factor = 0.4
    elseif mach < 8.0
        mach_factor = 0.6 + 0.4 * (mach - 2.0) / 6.0
    else
        mach_factor = 1.0
    end

    return base_reduction * mach_factor
end

function drag_reduction_factor(composite::CompositeDragReduction, mach, altitude)
    # Compound reductions multiplicatively (conservative)
    # C_d_eff = C_d × (1 - r₁) × (1 - r₂) × ...

    total_reduction = 0.0
    remaining_drag = 1.0

    for method in composite.methods
        r = drag_reduction_factor(method, mach, altitude)
        remaining_drag *= (1.0 - r)
    end

    total_reduction = 1.0 - remaining_drag
    return total_reduction
end

# =============================================================================
# MASS PENALTY CALCULATIONS
# =============================================================================

"""
    added_mass(method, capsule_diameter, capsule_mass)

Calculate added mass [kg] for drag reduction hardware.

# Arguments
- `method`: Drag reduction method
- `capsule_diameter`: Capsule diameter [m]
- `capsule_mass`: Baseline capsule mass [kg]

# Returns
- Added mass [kg]
"""
function added_mass(method::NoDragReduction, capsule_diameter, capsule_mass)
    return 0.0u"kg"
end

function added_mass(spike::Aerospike, capsule_diameter, capsule_mass)
    # Spike mass = π/4 × d² × L × ρ
    spike_length = spike.length_ratio * capsule_diameter
    volume = π/4 * spike.diameter^2 * spike_length
    mass = volume * spike.material_density

    return mass
end

function added_mass(disk::Aerodisk, capsule_diameter, capsule_mass)
    # Disk mass ≈ π × (D²/4) × thickness × ρ × num_disks
    # Assume thickness = 0.01 × capsule_diameter
    disk_diameter = disk.disk_diameter_ratio * capsule_diameter
    thickness = 0.01 * capsule_diameter

    volume_per_disk = π/4 * disk_diameter^2 * thickness
    mass = volume_per_disk * disk.num_disks * disk.material_density

    return mass
end

function added_mass(tailing::BoatTailing, capsule_diameter, capsule_mass)
    # Boat tailing is reshaping, not added mass
    # May actually reduce mass vs cylindrical
    # Return small negative or zero
    return 0.0u"kg"
end

function added_mass(rider::Waverider, capsule_diameter, capsule_mass)
    # Waverider is a different geometry, not added hardware
    # Mass change depends on design
    # Assume neutral for now
    return 0.0u"kg"
end

function added_mass(porous::PorousSurface, capsule_diameter, capsule_mass)
    # Porous ceramics/sintered metals slightly denser
    # Assume 2-5% increase proportional to porosity
    mass_fraction = 0.02 + 0.03 * porous.porosity / 0.2
    return mass_fraction * capsule_mass
end

function added_mass(composite::CompositeDragReduction, capsule_diameter, capsule_mass)
    # Sum all added masses
    total_mass = 0.0u"kg"
    for method in composite.methods
        total_mass += added_mass(method, capsule_diameter, capsule_mass)
    end
    return total_mass
end

# =============================================================================
# EFFECTIVE DRAG COEFFICIENT
# =============================================================================

"""
    effective_drag_coefficient(C_d_base, method, mach, altitude)

Calculate effective drag coefficient with drag reduction applied.

# Arguments
- `C_d_base`: Baseline drag coefficient (dimensionless)
- `method`: Drag reduction method
- `mach`: Mach number
- `altitude`: Altitude [m]

# Returns
- Effective drag coefficient (dimensionless)
"""
function effective_drag_coefficient(C_d_base, method::DragReductionMethod, mach, altitude)
    reduction = drag_reduction_factor(method, mach, altitude)
    return C_d_base * (1.0 - reduction)
end

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

"""
    print_drag_reduction_summary(method, capsule_diameter, capsule_mass, mach, altitude)

Print summary of drag reduction method performance.
"""
function print_drag_reduction_summary(method::DragReductionMethod, capsule_diameter, capsule_mass, mach, altitude)
    println("=" ^ 70)
    println("DRAG REDUCTION SUMMARY")
    println("=" ^ 70)
    println()
    println("Method: $(typeof(method))")
    println()

    # Reduction factor
    reduction = drag_reduction_factor(method, mach, altitude)
    @printf("Drag reduction: %.1f%%\n", reduction * 100)

    # Mass penalty
    mass_added = added_mass(method, capsule_diameter, capsule_mass)
    mass_pct = 100.0 * ustrip(u"kg", mass_added) / ustrip(u"kg", capsule_mass)
    @printf("Added mass: %.2f kg (%.1f%%)\n", ustrip(u"kg", mass_added), mass_pct)

    # Effective C_d
    C_d_base = 0.50  # Typical streamlined body
    C_d_eff = effective_drag_coefficient(C_d_base, method, mach, altitude)
    @printf("C_d: %.3f → %.3f\n", C_d_base, C_d_eff)

    println()
    println("=" ^ 70)
end

# Export types
export DragReductionMethod, NoDragReduction
export Aerospike, Aerodisk, BoatTailing, Waverider, PorousSurface
export CompositeDragReduction

# Export functions
export drag_reduction_factor, added_mass, effective_drag_coefficient
export print_drag_reduction_summary
