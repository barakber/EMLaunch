"""
# Electromagnetic Acceleration Module

Models the electromagnetic coilgun dynamics for the EM launcher.

## Physics

The launcher uses a series of electromagnetic coils that generate traveling
magnetic waves to accelerate a conducting payload capsule.

### Lorentz Force on Conducting Projectile

For a coilgun with a ferromagnetic or conducting projectile:

```math
F_{em} = \\frac{1}{2} \\frac{dL}{dx} I^2
```

where:
- `F_em` = electromagnetic force [N]
- `L(x)` = position-dependent inductance [H]
- `I` = coil current [A]
- `dL/dx` = inductance gradient [H/m]

### RLC Circuit Dynamics

Each coil is powered by a capacitor bank forming an RLC circuit:

```math
L \\frac{dI}{dt} + RI + \\frac{Q}{C} = V_0(t)
```

where:
- `L` = inductance [H]
- `R` = resistance [Ω]
- `C` = capacitance [F]
- `Q` = charge [C]
- `V_0(t)` = applied voltage [V]

### Coupled System

The full system couples the mechanical motion with electrical dynamics:

```math
\\begin{aligned}
\\frac{dx}{dt} &= v \\\\
\\frac{dv}{dt} &= \\frac{F_{em}(x,I)}{m} \\\\
\\frac{dI}{dt} &= \\frac{V_0(t) - Q/C - RI}{L(x)} \\\\
\\frac{dQ}{dt} &= I
\\end{aligned}
```
"""

using Unitful
using StaticArrays

"""
    CoilConfig

Configuration for a single electromagnetic coil in the launcher.

# Fields
- `position`: Position along launcher axis [m]
- `inductance`: Self-inductance of coil [H]
- `resistance`: Coil resistance [Ω]
- `capacitance`: Capacitor bank capacitance [F]
- `voltage`: Initial capacitor voltage [V]
- `trigger_time`: Time to fire this coil [s]
- `gradient`: Inductance gradient dL/dx [H/m]
"""
struct CoilConfig
    position
    inductance
    resistance
    capacitance
    voltage
    trigger_time
    gradient
end

"""
    LauncherConfig

Complete configuration of the electromagnetic launcher system.

# Fields
- `length`: Total launcher tube length [m]
- `coils`: Array of coil configurations
- `num_coils`: Number of coils (convenience)
- `vacuum_pressure_ratio`: Vacuum level in tube (1.0 = atmospheric, 0.001 = high vacuum, 0.0001 = ultra-high vacuum)
"""
struct LauncherConfig
    length
    coils::Vector{CoilConfig}
    num_coils::Int
    vacuum_pressure_ratio::Float64
end

"""
    create_uniform_launcher(;
        length,
        num_coils,
        inductance_per_coil,
        resistance_per_coil,
        capacitance_per_coil,
        voltage_per_coil,
        gradient_per_coil,
        vacuum_pressure_ratio=1.0
    )

Create a launcher with uniformly spaced coils having identical properties.

# Arguments
- `length`: Total launcher length [m]
- `num_coils`: Number of coils
- `inductance_per_coil`: Inductance for each coil [H]
- `resistance_per_coil`: Resistance for each coil [Ω]
- `capacitance_per_coil`: Capacitance for each coil [F]
- `voltage_per_coil`: Initial voltage for each coil [V]
- `gradient_per_coil`: Inductance gradient for each coil [H/m]
- `vacuum_pressure_ratio`: Vacuum level (1.0 = atmospheric, 0.001 = 99.9% vacuum, 0.0001 = 99.99% vacuum)

# Returns
- `LauncherConfig`
"""
function create_uniform_launcher(;
    length,
    num_coils::Integer,
    inductance_per_coil,
    resistance_per_coil,
    capacitance_per_coil,
    voltage_per_coil,
    gradient_per_coil,
    vacuum_pressure_ratio=1.0
)
    spacing = length / num_coils

    coils = [
        CoilConfig(
            i * spacing,                    # position along launcher axis
            inductance_per_coil,            # inductance
            resistance_per_coil,            # resistance
            capacitance_per_coil,           # capacitance
            voltage_per_coil,               # voltage
            0.0u"s",                        # trigger_time (not used with position-based triggering)
            gradient_per_coil               # gradient
        )
        for i in 0:(num_coils-1)
    ]

    return LauncherConfig(length, coils, num_coils, vacuum_pressure_ratio)
end

"""
    inductance_gradient(coil_inductance, coil_length)

Estimate the inductance gradient for a coil.

The inductance varies as the projectile enters/exits the coil.
This is a simplified model using a triangular profile.

# Arguments
- `coil_inductance`: Peak inductance [H]
- `coil_length`: Physical length of coil [m]

# Returns
- Inductance gradient [H/m]
"""
function inductance_gradient(coil_inductance, coil_length)
    # Peak gradient occurs at coil edges
    # Simplified model: gradient ~ L / (L_coil/2)
    return 2.0 * coil_inductance / coil_length
end

"""
    inductance_profile(x, coil_config)

Calculate the inductance as a function of projectile position.

Uses a triangular profile:
- Increases linearly as projectile enters coil
- Decreases linearly as projectile exits coil
- Zero when projectile is outside coil region

# Arguments
- `x`: Projectile position [m]
- `coil_config`: Configuration of the coil

# Returns
- Inductance at position x [H]
"""
function inductance_profile(x, coil_config::CoilConfig)
    # Assume coil has some spatial extent
    coil_length = 0.2u"m"  # Typical coil length
    x_coil = coil_config.position
    L_max = coil_config.inductance

    # Distance from coil center
    Δx = abs(x - x_coil)

    # Triangular profile
    if Δx < coil_length / 2
        # Inside coil region
        return L_max * (1.0 - 2.0 * Δx / coil_length)
    else
        # Outside coil region
        return 0.0u"H"
    end
end

"""
    em_force(x, I, coil_config)

Calculate electromagnetic force on projectile from a single coil.

```math
F_{em} = \\frac{1}{2} \\frac{dL}{dx} I^2
```

# Arguments
- `x`: Projectile position [m]
- `I`: Coil current [A]
- `coil_config`: Coil configuration

# Returns
- Electromagnetic force [N]
"""
function em_force(x, I, coil_config::CoilConfig)
    # Use the configured inductance gradient
    dL_dx = coil_config.gradient

    # Asymmetric range: longer when approaching, very short after passing
    # This matches real coilgun behavior where coils are turned off quickly after passing
    range_before = 6.0u"m"   # Attractive range before coil
    range_after = 0.5u"m"    # Very short range after coil center

    # Distance from coil center
    Δx = x - coil_config.position

    # Sign convention: positive force accelerates projectile forward
    if Δx < -range_before
        # Too far before coil
        return 0.0u"N"
    elseif Δx < 0.0u"m"
        # Approaching coil - attractive force (positive)
        # Force increases linearly as projectile gets closer
        return 0.5 * dL_dx * I^2 * (1.0 + ustrip(u"m", Δx) / ustrip(u"m", range_before))
    elseif Δx < range_after
        # Just past coil center - force drops quickly
        # Force decreases linearly to zero
        return 0.5 * dL_dx * I^2 * (1.0 - ustrip(u"m", Δx) / ustrip(u"m", range_after))
    else
        # Past the effective range
        return 0.0u"N"
    end
end

"""
    total_em_force(x, currents, launcher_config)

Calculate total electromagnetic force from all coils.

# Arguments
- `x`: Projectile position [m]
- `currents`: Array of current in each coil [A]
- `launcher_config`: Complete launcher configuration

# Returns
- Total electromagnetic force [N]
"""
function total_em_force(x, currents, launcher_config::LauncherConfig)
    F_total = 0.0u"N"

    for (coil, I) in zip(launcher_config.coils, currents)
        F_total += em_force(x, I, coil)
    end

    return F_total
end

"""
    coil_voltage(x_or_t, coil_config)

Determine applied voltage to coil based on projectile position or time.

# Position-based triggering:
Voltage is applied when projectile approaches within trigger distance of the coil.

# Time-based triggering:
Voltage is applied at or after the trigger_time specified in coil_config.

# Arguments
- `x_or_t`: Projectile position [m] or time [s]
- `coil_config`: Coil configuration

# Returns
- Applied voltage [V]
"""
function coil_voltage(x_or_t, coil_config::CoilConfig)
    # The capacitors start fully charged (Q = C×V). The "voltage" here represents
    # an external charging source that maintains the capacitor at its rated voltage.
    # When triggered, the source is disconnected (V=0) allowing the cap to discharge
    # through the coil's RLC circuit, generating current and EM force.
    #   NOT triggered: V = V0 (maintain capacitor charge, switch open)
    #   TRIGGERED:     V = 0  (disconnect supply, allow RLC discharge)

    # Check if this is time-based (has time dimension) or position-based
    if dimension(x_or_t) == dimension(1.0u"s")
        # Time-based triggering
        if x_or_t >= coil_config.trigger_time
            return 0.0u"V"  # Triggered: allow capacitor discharge
        else
            return coil_config.voltage  # Not yet triggered: maintain charge
        end
    else
        # Position-based triggering
        # Trigger when projectile is within trigger_distance of the coil
        # This gives the coil time to build up current before projectile arrives
        trigger_distance = 100.0u"m"  # Lead distance for coil activation

        # Handle both dimensionless and unitful position
        if typeof(x_or_t) <: AbstractFloat
            # Dimensionless position - strip units from coil position
            coil_pos = ustrip(u"m", coil_config.position)
            trigger_dist = ustrip(u"m", trigger_distance)

            if x_or_t >= (coil_pos - trigger_dist) && x_or_t <= (coil_pos + trigger_dist)
                return 0.0u"V"  # Triggered: allow capacitor discharge
            else
                return coil_config.voltage  # Not triggered: maintain charge
            end
        else
            # Unitful position - compare directly
            if x_or_t >= (coil_config.position - trigger_distance) && x_or_t <= (coil_config.position + trigger_distance)
                return 0.0u"V"  # Triggered: allow capacitor discharge
            else
                return coil_config.voltage  # Not triggered: maintain charge
            end
        end
    end
end

"""
    energy_efficiency(v_final, m_payload, launcher_config)

Calculate energy efficiency of the launcher system.

```math
\\eta = \\frac{KE_{payload}}{E_{electrical}}
```

# Arguments
- `v_final`: Final velocity of payload [m/s]
- `m_payload`: Payload mass [kg]
- `launcher_config`: Launcher configuration

# Returns
- Efficiency (dimensionless, 0 to 1)
"""
function energy_efficiency(v_final, m_payload, launcher_config::LauncherConfig)
    # Kinetic energy of payload
    KE = 0.5 * m_payload * v_final^2

    # Electrical energy stored in all capacitors
    E_elec = 0.0u"J"
    for coil in launcher_config.coils
        E_elec += 0.5 * coil.capacitance * coil.voltage^2
    end

    # Efficiency
    return ustrip(KE / E_elec)
end

"""
    required_capacitor_energy(v_final, m_payload, efficiency)

Calculate required total capacitor energy to reach target velocity.

# Arguments
- `v_final`: Target final velocity [m/s]
- `m_payload`: Payload mass [kg]
- `efficiency`: Expected system efficiency (dimensionless)

# Returns
- Required electrical energy [J]
"""
function required_capacitor_energy(v_final, m_payload, efficiency)
    KE = 0.5 * m_payload * v_final^2
    return KE / efficiency
end

"""
    print_launcher_config(config::LauncherConfig)

Print human-readable summary of launcher configuration.
"""
function print_launcher_config(config::LauncherConfig)
    println("ELECTROMAGNETIC LAUNCHER CONFIGURATION")
    println("=" ^ 60)
    println("Total length:        $(config.length)")
    println("Number of coils:     $(config.num_coils)")
    println("Coil spacing:        $(config.length / config.num_coils)")
    println()
    println("TYPICAL COIL PARAMETERS:")
    if length(config.coils) > 0
        coil = config.coils[1]
        println("  Inductance:        $(coil.inductance)")
        println("  Resistance:        $(coil.resistance)")
        println("  Capacitance:       $(coil.capacitance)")
        println("  Voltage:           $(coil.voltage)")
        println("  Inductance grad:   $(coil.gradient)")

        E_per_coil = 0.5 * coil.capacitance * coil.voltage^2
        E_total = config.num_coils * E_per_coil

        println()
        println("ENERGY:")
        println("  Energy per coil:   $(E_per_coil)")
        println("  Total energy:      $(E_total)")
        println("  Total energy (MJ): $(ustrip(u"MJ", E_total)) MJ")
    end
    println("=" ^ 60)
end

export CoilConfig, LauncherConfig
export create_uniform_launcher
export inductance_gradient, inductance_profile
export em_force, total_em_force, coil_voltage
export energy_efficiency, required_capacitor_energy
export print_launcher_config
