"""
# Thermal Protection & Ablation Module

Models ablative heat shields, thermal protection systems (TPS), and thermal survival.

## Ablation Physics

When heat flux exceeds material limits, ablative materials vaporize/char, carrying away heat:

```math
\\dot{m}_{abl} = \\frac{\\dot{q}_{net}}{H_{abl}}
```

where:
- `ṁ_abl` = ablation rate [kg/(m²·s)]
- `q̇_net` = net heat flux [W/m²]
- `H_abl` = effective heat of ablation [J/kg]

## Heat Shield Mass

Required heat shield mass depends on total heat load:

```math
m_{shield} = \\frac{Q_{total}}{H_{abl}} \\cdot A
```

where:
- `Q_total` = ∫ q̇ dt = total heat load [J/m²]
- `A` = protected area [m²]

## Common Ablative Materials

| Material | H_abl [MJ/kg] | ρ [kg/m³] | Max Temp [K] | Notes |
|----------|---------------|-----------|--------------|-------|
| Carbon Phenolic | 8-12 | 1400 | 3800 | Apollo heat shield |
| PICA | 10-15 | 250 | 3000 | Mars entry, Stardust |
| AVCOAT | 8-10 | 500 | 2800 | Orion capsule |
| Cork | 3-5 | 240 | 1800 | Early reentry vehicles |

## Thermal Survival Criteria

1. **Heat Shield Intact**: Remaining thickness > 0
2. **Temperature Limits**: Structure temp < T_max
3. **Heat Flux Limits**: q̇ < q̇_max (for radiative cooling)
4. **Total Heat Load**: Integrated heat < shield capacity
"""

using Unitful
using Printf
using StaticArrays

# =============================================================================
# ABLATIVE MATERIAL PROPERTIES
# =============================================================================

"""
Ablative material properties for thermal protection.

# Fields
- `name`: Material name
- `H_abl`: Effective heat of ablation [J/kg]
- `density`: Material density [kg/m³]
- `k_thermal`: Thermal conductivity [W/(m·K)]
- `T_max`: Maximum substrate temperature [K]
- `emissivity`: Thermal emissivity for radiation (dimensionless)
- `char_fraction`: Fraction that forms protective char layer (0-1)
"""
struct AblativeMaterial
    name::String
    H_abl::typeof(1.0u"J/kg")           # Heat of ablation
    density::typeof(1.0u"kg/m^3")       # Material density
    k_thermal::typeof(1.0u"W/(m*K)")    # Thermal conductivity
    T_max::typeof(1.0u"K")              # Max substrate temperature
    emissivity::Float64                  # Thermal emissivity
    char_fraction::Float64               # Char formation fraction
end

"""
Common ablative materials library.
"""
const ABLATIVE_MATERIALS = Dict(
    :carbon_phenolic => AblativeMaterial(
        "Carbon Phenolic",
        10.0e6u"J/kg",      # 10 MJ/kg
        1400.0u"kg/m^3",
        2.0u"W/(m*K)",
        3800.0u"K",
        0.85,
        0.60                # 60% forms char
    ),
    :pica => AblativeMaterial(
        "PICA (Phenolic Impregnated Carbon Ablator)",
        12.5e6u"J/kg",      # 12.5 MJ/kg (highly efficient)
        250.0u"kg/m^3",     # Very lightweight
        0.5u"W/(m*K)",
        3000.0u"K",
        0.90,
        0.70
    ),
    :avcoat => AblativeMaterial(
        "AVCOAT",
        9.0e6u"J/kg",       # 9 MJ/kg
        500.0u"kg/m^3",
        0.8u"W/(m*K)",
        2800.0u"K",
        0.85,
        0.50
    ),
    :cork => AblativeMaterial(
        "Cork/Epoxy",
        4.0e6u"J/kg",       # 4 MJ/kg (lower performance)
        240.0u"kg/m^3",
        0.15u"W/(m*K)",
        1800.0u"K",
        0.75,
        0.30
    )
)

# =============================================================================
# HEAT SHIELD CONFIGURATION
# =============================================================================

"""
Heat shield configuration and state.

# Fields
- `material`: Ablative material properties
- `initial_thickness`: Initial heat shield thickness [m]
- `area`: Protected surface area [m²]
- `thickness`: Current remaining thickness [m] (state variable)
- `mass`: Current heat shield mass [kg] (state variable)
- `total_heat_absorbed`: Cumulative heat absorbed [J] (state variable)
"""
mutable struct HeatShield
    material::AblativeMaterial
    initial_thickness::typeof(1.0u"m")
    area::typeof(1.0u"m^2")
    # State variables (updated during simulation)
    thickness::typeof(1.0u"m")
    mass::typeof(1.0u"kg")
    total_heat_absorbed::typeof(1.0u"J")
end

"""
    HeatShield(material, thickness, area)

Create a heat shield with given material, thickness, and area.

# Arguments
- `material`: Ablative material (symbol from ABLATIVE_MATERIALS or AblativeMaterial)
- `thickness`: Initial thickness [m]
- `area`: Protected surface area [m²]

# Returns
- HeatShield object
"""
function HeatShield(material, thickness, area)
    mat = isa(material, Symbol) ? ABLATIVE_MATERIALS[material] : material

    mass = mat.density * thickness * area

    return HeatShield(
        mat,
        thickness,
        area,
        thickness,  # Initial thickness
        mass,       # Initial mass
        0.0u"J"     # No heat absorbed yet
    )
end

# =============================================================================
# ABLATION CALCULATIONS
# =============================================================================

"""
    ablation_rate(q_dot, material)

Calculate ablation rate for given heat flux and material.

```math
\\dot{m}_{abl} = \\frac{\\dot{q}_{net}}{H_{abl}}
```

# Arguments
- `q_dot`: Net heat flux [W/m²]
- `material`: Ablative material

# Returns
- Ablation rate [kg/(m²·s)]
"""
function ablation_rate(q_dot, material::AblativeMaterial)
    # Only ablate for positive heat flux
    if q_dot <= 0.0u"W/m^2"
        return 0.0u"kg/(m^2*s)"
    end

    # Ablation rate = heat flux / heat of ablation
    m_dot = q_dot / material.H_abl

    return m_dot
end

"""
    ablation_thickness_rate(q_dot, material)

Calculate rate of thickness loss due to ablation.

# Arguments
- `q_dot`: Net heat flux [W/m²]
- `material`: Ablative material

# Returns
- Thickness loss rate [m/s]
"""
function ablation_thickness_rate(q_dot, material::AblativeMaterial)
    m_dot = ablation_rate(q_dot, material)

    # Convert mass rate to thickness rate
    # ṁ = ρ · Ȧ, where A = area · thickness
    # For constant area: ṁ = ρ · area · (dh/dt)
    # Therefore: dh/dt = ṁ / (ρ · area) = ṁ / ρ  (per unit area)

    h_dot = m_dot / material.density

    return h_dot
end

"""
    required_heat_shield_mass(total_heat_load, area, material)

Calculate required heat shield mass for a given total heat load.

# Arguments
- `total_heat_load`: Total integrated heat load [J/m²]
- `area`: Protected surface area [m²]
- `material`: Ablative material (symbol or AblativeMaterial)

# Returns
- Required heat shield mass [kg]
"""
function required_heat_shield_mass(total_heat_load, area, material)
    mat = isa(material, Symbol) ? ABLATIVE_MATERIALS[material] : material

    # Mass ablated = Q / H_abl
    m_ablated = total_heat_load * area / mat.H_abl

    # Add safety factor (typically 1.5-2.0x)
    safety_factor = 1.5
    m_required = m_ablated * safety_factor

    return m_required
end

"""
    required_heat_shield_thickness(total_heat_load, material)

Calculate required heat shield thickness for a given total heat load.

# Arguments
- `total_heat_load`: Total integrated heat load [J/m²]
- `material`: Ablative material (symbol or AblativeMaterial)

# Returns
- Required thickness [m]
"""
function required_heat_shield_thickness(total_heat_load, material)
    mat = isa(material, Symbol) ? ABLATIVE_MATERIALS[material] : material

    # Thickness ablated = Q / (H_abl · ρ)
    h_ablated = total_heat_load / (mat.H_abl * mat.density)

    # Add safety factor
    safety_factor = 1.5
    h_required = h_ablated * safety_factor

    return h_required
end

# =============================================================================
# HEAT SHIELD DYNAMICS (for integration into trajectory)
# =============================================================================

"""
    update_heat_shield!(shield, q_dot, dt)

Update heat shield state due to ablation over time step dt.

# Arguments
- `shield`: HeatShield object (modified in place)
- `q_dot`: Heat flux [W/m²]
- `dt`: Time step [s]

# Returns
- `(mass_lost, heat_absorbed)`: Tuple of mass lost [kg] and heat absorbed [J]
"""
function update_heat_shield!(shield::HeatShield, q_dot, dt)
    # Calculate ablation rate
    h_dot = ablation_thickness_rate(q_dot, shield.material)

    # Update thickness
    Δh = h_dot * dt
    shield.thickness = max(0.0u"m", shield.thickness - Δh)

    # Calculate mass lost
    Δm = shield.material.density * Δh * shield.area
    shield.mass = max(0.0u"kg", shield.mass - Δm)

    # Calculate heat absorbed
    ΔQ = q_dot * shield.area * dt
    shield.total_heat_absorbed += ΔQ

    return (mass_lost=Δm, heat_absorbed=ΔQ)
end

"""
    heat_shield_intact(shield)

Check if heat shield is still intact (thickness > 0).

# Arguments
- `shield`: HeatShield object

# Returns
- `true` if intact, `false` if burned through
"""
function heat_shield_intact(shield::HeatShield)
    return shield.thickness > 0.0u"m"
end

"""
    heat_shield_margin(shield)

Calculate heat shield margin (remaining thickness / initial thickness).

# Arguments
- `shield`: HeatShield object

# Returns
- Margin ratio (0.0 = completely ablated, 1.0 = pristine)
"""
function heat_shield_margin(shield::HeatShield)
    return ustrip(shield.thickness / shield.initial_thickness)
end

# =============================================================================
# THERMAL SURVIVAL ANALYSIS
# =============================================================================

"""
    thermal_survival_check(q_dot, T_structure, shield, dt)

Check multiple thermal survival criteria.

# Arguments
- `q_dot`: Current heat flux [W/m²]
- `T_structure`: Current structure temperature [K]
- `shield`: HeatShield object (can be nothing for no shield)
- `dt`: Time step [s]

# Returns
Named tuple with:
- `survived`: Overall survival status (boolean)
- `shield_intact`: Heat shield status (boolean)
- `temp_ok`: Temperature within limits (boolean)
- `margin`: Heat shield margin (0-1)
- `failure_mode`: Reason for failure (if any)
"""
function thermal_survival_check(q_dot, T_structure, shield, dt)
    # No heat shield case
    if isnothing(shield)
        # Check if structure can survive via radiative cooling alone
        # Most materials fail above 1500-2000 K
        T_limit = 1500.0u"K"
        temp_ok = T_structure < T_limit

        return (
            survived = temp_ok,
            shield_intact = false,
            temp_ok = temp_ok,
            margin = 0.0,
            failure_mode = temp_ok ? nothing : "Structure temperature exceeded"
        )
    end

    # With heat shield
    shield_ok = heat_shield_intact(shield)
    temp_ok = T_structure < shield.material.T_max
    margin = heat_shield_margin(shield)

    survived = shield_ok && temp_ok

    failure_mode = nothing
    if !shield_ok
        failure_mode = "Heat shield burned through"
    elseif !temp_ok
        failure_mode = "Structure temperature exceeded $(shield.material.T_max)"
    end

    return (
        survived = survived,
        shield_intact = shield_ok,
        temp_ok = temp_ok,
        margin = margin,
        failure_mode = failure_mode
    )
end

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

"""
    print_heat_shield_summary(shield)

Print summary of heat shield status.
"""
function print_heat_shield_summary(shield::HeatShield)
    println("=" ^ 70)
    println("HEAT SHIELD SUMMARY: $(shield.material.name)")
    println("=" ^ 70)
    println()
    @printf("Initial thickness:    %.4f m (%.2f cm)\n",
            ustrip(u"m", shield.initial_thickness),
            ustrip(u"cm", shield.initial_thickness))
    @printf("Remaining thickness:  %.4f m (%.2f cm)\n",
            ustrip(u"m", shield.thickness),
            ustrip(u"cm", shield.thickness))
    @printf("Thickness ablated:    %.4f m (%.2f cm)\n",
            ustrip(u"m", shield.initial_thickness - shield.thickness),
            ustrip(u"cm", shield.initial_thickness - shield.thickness))
    @printf("Margin:               %.1f%%\n", heat_shield_margin(shield) * 100)
    println()
    @printf("Initial mass:         %.3f kg\n", ustrip(u"kg", shield.material.density * shield.initial_thickness * shield.area))
    @printf("Remaining mass:       %.3f kg\n", ustrip(u"kg", shield.mass))
    @printf("Mass ablated:         %.3f kg\n",
            ustrip(u"kg", shield.material.density * shield.initial_thickness * shield.area - shield.mass))
    println()
    @printf("Total heat absorbed:  %.2f MJ\n", ustrip(u"MJ", shield.total_heat_absorbed))
    @printf("Heat load per area:   %.2f MJ/m²\n",
            ustrip(u"MJ/m^2", shield.total_heat_absorbed / shield.area))
    println()

    if heat_shield_intact(shield)
        println("Status: ✓ INTACT")
    else
        println("Status: ✗ BURNED THROUGH")
    end
    println("=" ^ 70)
end

"""
    print_ablative_materials()

Print table of available ablative materials and their properties.
"""
function print_ablative_materials()
    println("=" ^ 80)
    println("ABLATIVE MATERIALS LIBRARY")
    println("=" ^ 80)
    println()
    @printf("%-30s %-15s %-15s %-15s\n",
            "Material", "H_abl [MJ/kg]", "Density [kg/m³]", "T_max [K]")
    println("-" ^ 80)

    for (key, mat) in ABLATIVE_MATERIALS
        @printf("%-30s %-15.1f %-15.0f %-15.0f\n",
                mat.name,
                ustrip(u"MJ/kg", mat.H_abl),
                ustrip(u"kg/m^3", mat.density),
                ustrip(u"K", mat.T_max))
    end
    println("=" ^ 80)
end

# Export types
export AblativeMaterial, HeatShield, ABLATIVE_MATERIALS

# Export functions
export ablation_rate, ablation_thickness_rate
export required_heat_shield_mass, required_heat_shield_thickness
export update_heat_shield!, heat_shield_intact, heat_shield_margin
export thermal_survival_check
export print_heat_shield_summary, print_ablative_materials
