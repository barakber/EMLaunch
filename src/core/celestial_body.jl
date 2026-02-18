"""
# Celestial Body Module

Defines types and constants for multi-body support (Earth, Moon, etc.).
All body-specific physical parameters are encapsulated in `CelestialBody`,
allowing the same simulation engine to work for any celestial body.

## Type Hierarchy

```
AtmosphereModel (abstract)
├── VacuumAtmosphere          — no atmosphere (Moon, asteroids)
└── StandardAtmosphere1976    — layered Earth atmosphere model
```

## Usage

All physics functions accept a `body` keyword argument defaulting to `EARTH`:
```julia
gravity_point_mass(position; body=EARTH)
gravity_point_mass(position; body=MOON)
```
"""

using Unitful

# =============================================================================
# ATMOSPHERE TYPES
# =============================================================================

"""
    AtmosphereModel

Abstract type for atmosphere models. Dispatch on subtypes to select
vacuum vs. atmospheric physics.
"""
abstract type AtmosphereModel end

"""
    VacuumAtmosphere <: AtmosphereModel

Represents the absence of atmosphere (e.g., Moon, asteroids).
All atmospheric drag and heating functions return zero.
"""
struct VacuumAtmosphere <: AtmosphereModel end

"""
    StandardAtmosphere1976 <: AtmosphereModel

US Standard Atmosphere 1976 layered model. Stores all atmospheric
constants needed by the atmosphere module.

# Fields
- `layers`: Vector of (base_altitude, base_temperature, base_pressure, lapse_rate) tuples
- `sea_level_density`: [kg/m³]
- `sea_level_pressure`: [Pa]
- `sea_level_temperature`: [K]
- `scale_height`: [m]
- `gas_constant_specific`: [J/(kg·K)]
- `heat_ratio`: γ (dimensionless)
- `molar_mass`: [kg/mol]
- `dynamic_viscosity_ref`: [Pa·s]
- `viscosity_ref_temp`: [K]
- `sutherland_constant`: [K]
"""
struct StandardAtmosphere1976 <: AtmosphereModel
    layers::Vector
    sea_level_density
    sea_level_pressure
    sea_level_temperature
    scale_height
    gas_constant_specific
    heat_ratio::Float64
    molar_mass
    dynamic_viscosity_ref
    viscosity_ref_temp
    sutherland_constant
end

# =============================================================================
# CELESTIAL BODY TYPE
# =============================================================================

"""
    CelestialBody

Encapsulates all body-specific physical parameters for a celestial body.

# Fields
- `name`: Human-readable name (e.g., "Earth", "Moon")
- `mu`: Gravitational parameter GM [m³/s²]
- `radius`: Equatorial radius [m]
- `J2`: Oblateness coefficient (dimensionless Float64)
- `rotation_rate`: Rotation rate [rad/s]
- `surface_gravity`: Surface gravitational acceleration [m/s²]
- `atmosphere`: AtmosphereModel (VacuumAtmosphere or StandardAtmosphere1976)
"""
struct CelestialBody
    name::String
    mu
    radius
    J2::Float64
    rotation_rate
    surface_gravity
    atmosphere::AtmosphereModel
end

"""
    has_atmosphere(body::CelestialBody) -> Bool

Return `true` if the body has an atmosphere (not vacuum).
"""
has_atmosphere(body::CelestialBody) = !(body.atmosphere isa VacuumAtmosphere)

# =============================================================================
# EARTH
# =============================================================================

const EARTH_ATMOSPHERE = StandardAtmosphere1976(
    # Layers: (h_base, T_base, P_base, lapse_rate)
    [
        (0.0u"m",      288.15u"K", 101325.0u"Pa", -0.0065u"K/m"),  # Troposphere
        (11000.0u"m",  216.65u"K", 22632.1u"Pa",   0.0u"K/m"),     # Tropopause
        (20000.0u"m",  216.65u"K", 5474.89u"Pa",   0.001u"K/m"),   # Stratosphere 1
        (32000.0u"m",  228.65u"K", 868.019u"Pa",   0.0028u"K/m"),  # Stratosphere 2
        (47000.0u"m",  270.65u"K", 110.906u"Pa",   0.0u"K/m"),     # Stratopause
        (51000.0u"m",  270.65u"K", 66.9389u"Pa",  -0.0028u"K/m"),  # Mesosphere 1
        (71000.0u"m",  214.65u"K", 3.95642u"Pa",  -0.002u"K/m"),   # Mesosphere 2
    ],
    1.225u"kg/m^3",           # sea_level_density
    101325.0u"Pa",            # sea_level_pressure
    288.15u"K",               # sea_level_temperature
    8500.0u"m",               # scale_height
    287.05u"J/(kg*K)",        # gas_constant_specific (R for air)
    1.4,                      # heat_ratio (γ for air)
    0.0289644u"kg/mol",      # molar_mass
    1.789e-5u"Pa*s",         # dynamic_viscosity_ref
    288.15u"K",               # viscosity_ref_temp
    110.4u"K"                 # sutherland_constant
)

const EARTH = CelestialBody(
    "Earth",
    3.986004418e14u"m^3/s^2",   # mu (GM)
    6.3781370e6u"m",             # radius (equatorial)
    1.08263e-3,                  # J2
    7.2921159e-5u"rad/s",       # rotation_rate
    9.80665u"m/s^2",            # surface_gravity
    EARTH_ATMOSPHERE
)

# =============================================================================
# MOON
# =============================================================================

const MOON = CelestialBody(
    "Moon",
    4.9048695e12u"m^3/s^2",     # mu (GM)
    1.7374e6u"m",                # radius (mean)
    2.033e-4,                    # J2 (much smaller than Earth)
    2.6617e-6u"rad/s",          # rotation_rate (synchronous)
    1.625u"m/s^2",              # surface_gravity
    VacuumAtmosphere()
)

export AtmosphereModel, VacuumAtmosphere, StandardAtmosphere1976
export CelestialBody, has_atmosphere
export EARTH, EARTH_ATMOSPHERE, MOON
