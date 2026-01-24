"""
# Atmosphere Module

Models atmospheric properties, drag forces, and aerodynamic heating
using the US Standard Atmosphere 1976 model with type-safe units.

## Atmospheric Density

Exponential model (simplified):
```math
\\rho(h) = \\rho_0 \\exp\\left(-\\frac{h}{H}\\right)
```

Layered model (US Standard Atmosphere 1976):
```math
\\rho(h) = \\rho_b \\left(\\frac{T(h)}{T_b}\\right)^{\\left(1 + \\frac{g_0 M}{R T_b}\\right)}
```

## Drag Force

```math
F_d = \\frac{1}{2} \\rho(h) v^2 C_d A
```

where:
- `ρ(h)` = atmospheric density at altitude h [kg/m³]
- `v` = velocity [m/s]
- `C_d` = drag coefficient (Mach-dependent)
- `A` = reference area [m²]

## Aerodynamic Heating

Detra-Kemp-Riddell correlation for stagnation point heating:
```math
\\dot{q} = K \\sqrt{\\frac{\\rho}{R_n}} v^3
```

where:
- `K ≈ 1.83×10⁻⁴` [kg^0.5/m] for Earth
- `R_n` = nose radius [m]
- `q̇` = heat flux [W/m²]
"""

using Unitful
using Printf
using StaticArrays
using NaNMath

# =============================================================================
# ATMOSPHERIC LAYERS (US Standard Atmosphere 1976)
# =============================================================================

"""
US Standard Atmosphere 1976 layer definitions.

Each layer has:
- Base altitude [m]
- Base temperature [K]
- Base pressure [Pa]
- Temperature lapse rate [K/m]
"""
const ATMOSPHERE_LAYERS = [
    # h_base,    T_base,      P_base,         lapse_rate
    (0.0u"m",    288.15u"K",  101325.0u"Pa",  -0.0065u"K/m"),  # Troposphere
    (11000.0u"m", 216.65u"K", 22632.1u"Pa",   0.0u"K/m"),      # Tropopause
    (20000.0u"m", 216.65u"K", 5474.89u"Pa",   0.001u"K/m"),    # Stratosphere 1
    (32000.0u"m", 228.65u"K", 868.019u"Pa",   0.0028u"K/m"),   # Stratosphere 2
    (47000.0u"m", 270.65u"K", 110.906u"Pa",   0.0u"K/m"),      # Stratopause
    (51000.0u"m", 270.65u"K", 66.9389u"Pa",   -0.0028u"K/m"),  # Mesosphere 1
    (71000.0u"m", 214.65u"K", 3.95642u"Pa",   -0.002u"K/m"),   # Mesosphere 2
]

"""
    get_layer_index(altitude)

Get the atmospheric layer index for a given altitude.

# Arguments
- `altitude`: Altitude above sea level [m]

# Returns
- Layer index (1 to 7)
"""
function get_layer_index(altitude)
    for i in length(ATMOSPHERE_LAYERS):-1:1
        if altitude >= ATMOSPHERE_LAYERS[i][1]
            return i
        end
    end
    return 1  # Below lowest layer
end

"""
    temperature(altitude)

Calculate atmospheric temperature at given altitude using US Standard Atmosphere 1976.

# Arguments
- `altitude`: Altitude above sea level [m]

# Returns
- Temperature [K]
"""
function temperature(altitude)
    idx = get_layer_index(altitude)
    h_base, T_base, _, lapse_rate = ATMOSPHERE_LAYERS[idx]

    # Temperature at altitude
    T = T_base + lapse_rate * (altitude - h_base)

    # For extreme altitudes above 86 km, use constant temperature model
    # Use the temperature at 86 km as the constant value to avoid discontinuity
    if altitude > 86000.0u"m"
        # Calculate T at 86 km using the last layer
        T_86km = 214.65u"K" + (-0.002u"K/m") * (86000.0u"m" - 71000.0u"m")
        return max(T_86km, 50.0u"K")  # Minimum 50 K
    end

    # Safety check: never return negative temperature
    if T < 50.0u"K"
        return 50.0u"K"
    end

    return T
end

"""
    pressure(altitude)

Calculate atmospheric pressure at given altitude using US Standard Atmosphere 1976.

# Arguments
- `altitude`: Altitude above sea level [m]

# Returns
- Pressure [Pa]
"""
function pressure(altitude)
    # Constants
    g = 9.80665u"m/s^2"
    M = 0.0289644u"kg/mol"  # Molar mass of air
    R = 8.31432u"J/(mol*K)"  # Universal gas constant

    # For extreme altitudes above 86 km, use exponential decay
    # First calculate reference pressure at 86 km for continuity
    if altitude > 86000.0u"m"
        # Calculate pressure at 86 km using the layer model
        h_86km = 86000.0u"m"
        idx_86 = get_layer_index(h_86km)
        h_base, T_base, P_base, lapse_rate = ATMOSPHERE_LAYERS[idx_86]

        # Calculate P_86km
        T_86 = temperature(h_86km)
        if abs(ustrip(u"K/m", lapse_rate)) < 1e-10
            P_86km = P_base * exp(-g * M * (h_86km - h_base) / (R * T_base))
        else
            P_86km = P_base * (T_86 / T_base)^(-g * M / (R * lapse_rate))
        end

        # Exponential decay above 86 km
        H_scale = 8500.0u"m"
        return P_86km * exp(-(altitude - 86000.0u"m") / H_scale)
    end

    idx = get_layer_index(altitude)
    h_base, T_base, P_base, lapse_rate = ATMOSPHERE_LAYERS[idx]

    if abs(ustrip(u"K/m", lapse_rate)) < 1e-10  # Isothermal layer
        # Exponential model
        P = P_base * exp(-g * M * (altitude - h_base) / (R * T_base))
    else  # Temperature gradient
        # Power law
        T = temperature(altitude)  # Use the safe temperature function
        # Ensure positive temperature ratio
        if T > 0.0u"K" && T_base > 0.0u"K"
            P = P_base * (T / T_base)^(-g * M / (R * lapse_rate))
        else
            # Fallback to exponential if temperature is problematic
            H_scale = 8500.0u"m"
            P = P_base * exp(-(altitude - h_base) / H_scale)
        end
    end

    return P
end

"""
    density(altitude)

Calculate atmospheric density at given altitude.

Uses ideal gas law: ρ = P/(R_specific * T)

# Arguments
- `altitude`: Altitude above sea level [m]

# Returns
- Density [kg/m³]
"""
function density(altitude)
    P = pressure(altitude)
    T = temperature(altitude)
    R_specific = 287.05u"J/(kg*K)"  # Specific gas constant for air

    ρ = P / (R_specific * T)
    # Force unit simplification to kg/m³
    return uconvert(u"kg/m^3", ρ)
end

"""
    density_exponential(altitude, H_scale=8500.0u"m")

Simplified exponential atmospheric density model.

```math
\\rho(h) = \\rho_0 \\exp(-h/H)
```

Faster but less accurate than US Standard Atmosphere.

# Arguments
- `altitude`: Altitude above sea level [m]
- `H_scale`: Scale height [m] (default: 8500 m)

# Returns
- Density [kg/m³]
"""
function density_exponential(altitude, H_scale=8500.0u"m")
    ρ_0 = 1.225u"kg/m^3"
    return ρ_0 * exp(-altitude / H_scale)
end

"""
    speed_of_sound(altitude)

Calculate speed of sound at given altitude.

```math
a = \\sqrt{\\gamma R T}
```

# Arguments
- `altitude`: Altitude above sea level [m]

# Returns
- Speed of sound [m/s]
"""
function speed_of_sound(altitude)
    T = temperature(altitude)
    γ = 1.4  # Specific heat ratio for air
    R = 287.05u"J/(kg*K)"

    # Use NaNMath.sqrt to handle invalid values gracefully
    # Returns NaN for negative arguments, allowing SDE solver to reject the step
    arg = γ * R * T
    a = NaNMath.sqrt(ustrip(u"m^2/s^2", arg)) * u"m/s"
    return a
end

"""
    mach_number(velocity, altitude)

Calculate Mach number at given velocity and altitude.

# Arguments
- `velocity`: Velocity [m/s]
- `altitude`: Altitude [m]

# Returns
- Mach number (dimensionless)
"""
function mach_number(velocity, altitude)
    a = speed_of_sound(altitude)
    return ustrip(velocity / a)
end

"""
    drag_coefficient(mach, regime=:hypersonic)

Calculate drag coefficient based on Mach number.

# Arguments
- `mach`: Mach number (dimensionless)
- `regime`: Flow regime (:subsonic, :transonic, :supersonic, :hypersonic)

# Returns
- Drag coefficient C_d (dimensionless)
"""
function drag_coefficient(mach; regime=:auto)
    # Auto-detect regime
    if regime == :auto
        if mach < 0.8
            regime = :subsonic
        elseif mach < 1.2
            regime = :transonic
        elseif mach < 5.0
            regime = :supersonic
        else
            regime = :hypersonic
        end
    end

    # Drag coefficient correlations
    if regime == :subsonic
        return 0.40  # Streamlined body
    elseif regime == :transonic
        # Drag rise in transonic regime
        return 0.40 + 0.45 * (mach - 0.8) / 0.4  # Linear interpolation
    elseif regime == :supersonic
        # Modified Newtonian
        return 0.85
    else  # hypersonic
        return 0.90  # Slight increase at hypersonic speeds
    end
end

"""
    drag_force(velocity, altitude, area, C_d=nothing)

Calculate aerodynamic drag force.

```math
F_d = \\frac{1}{2} \\rho v^2 C_d A
```

# Arguments
- `velocity`: Velocity vector [m/s] (can be scalar or 3D vector)
- `altitude`: Altitude [m]
- `area`: Reference area [m²]
- `C_d`: Drag coefficient (if nothing, auto-calculated from Mach)

# Returns
- Drag force magnitude [N] (or vector if velocity is vector)
"""
function drag_force(velocity, altitude, area, C_d=nothing)
    ρ = density(altitude)

    # Calculate drag coefficient if not provided
    if isnothing(C_d)
        v_mag = isa(velocity, Number) ? velocity : norm(velocity)
        M = mach_number(v_mag, altitude)
        C_d = drag_coefficient(M)
    end

    # Drag force (always opposes velocity)
    if isa(velocity, Number)
        # Scalar velocity
        F_d = 0.5 * ρ * velocity^2 * C_d * area
        return F_d
    else
        # Vector velocity
        v_mag = norm(velocity)
        # Handle zero velocity case to avoid NaN
        if v_mag < 1e-10u"m/s"
            return SVector(0.0u"N", 0.0u"N", 0.0u"N")
        end
        F_d_mag = 0.5 * ρ * v_mag^2 * C_d * area
        # Return vector opposing velocity
        return -F_d_mag * velocity / v_mag
    end
end

"""
    aerodynamic_heating(velocity, altitude, nose_radius)

Calculate stagnation point heating rate using Detra-Kemp-Riddell correlation.

```math
\\dot{q} = K \\sqrt{\\frac{\\rho}{R_n}} v^3
```

# Arguments
- `velocity`: Velocity [m/s]
- `altitude`: Altitude [m]
- `nose_radius`: Nose radius [m]

# Returns
- Heat flux [W/m²]
"""
function aerodynamic_heating(velocity, altitude, nose_radius)
    K = 1.83e-4u"kg^0.5/m"
    ρ = density(altitude)

    # Handle vector velocity
    v = isa(velocity, Number) ? velocity : norm(velocity)

    q_dot = K * sqrt(ρ / nose_radius) * v^3

    return q_dot
end

"""
    reynolds_number(velocity, altitude, length_scale)

Calculate Reynolds number.

```math
Re = \\frac{\\rho v L}{\\mu}
```

# Arguments
- `velocity`: Velocity [m/s]
- `altitude`: Altitude [m]
- `length_scale`: Characteristic length [m]

# Returns
- Reynolds number (dimensionless)
"""
function reynolds_number(velocity, altitude, length_scale)
    ρ = density(altitude)
    T = temperature(altitude)

    # Sutherland's formula for dynamic viscosity
    μ_ref = 1.789e-5u"Pa*s"
    T_ref = 288.15u"K"
    S = 110.4u"K"  # Sutherland constant

    μ = μ_ref * (T / T_ref)^1.5 * (T_ref + S) / (T + S)

    # Handle vector velocity
    v = isa(velocity, Number) ? velocity : norm(velocity)

    Re = ρ * v * length_scale / μ

    return ustrip(Re)  # Dimensionless
end

"""
    print_atmosphere_profile(altitudes)

Print atmospheric properties at various altitudes.

# Arguments
- `altitudes`: Array of altitudes [m]
"""
function print_atmosphere_profile(altitudes)
    println("=" ^ 80)
    println("ATMOSPHERIC PROFILE (US Standard Atmosphere 1976)")
    println("=" ^ 80)
    println()
    @printf("%-10s %-12s %-12s %-12s %-12s\n",
            "Alt [km]", "T [K]", "P [Pa]", "ρ [kg/m³]", "a [m/s]")
    println("-" ^ 80)

    for h in altitudes
        T = temperature(h)
        P = pressure(h)
        ρ = density(h)
        a = speed_of_sound(h)

        @printf("%-10.1f %-12.2f %-12.2f %-12.6f %-12.2f\n",
                ustrip(u"km", h),
                ustrip(u"K", T),
                ustrip(u"Pa", P),
                ustrip(u"kg/m^3", ρ),
                ustrip(u"m/s", a))
    end
    println("=" ^ 80)
end

export temperature, pressure, density, density_exponential, speed_of_sound
export mach_number, drag_coefficient, drag_force
export aerodynamic_heating, reynolds_number
export print_atmosphere_profile
