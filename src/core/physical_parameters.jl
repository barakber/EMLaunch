"""
# Physical Parameters Module

All physical constants and material properties with type-safe units using Unitful.jl.
This ensures dimensional consistency throughout all calculations.
"""

using Unitful
using UnitfulAstro
using PhysicalConstants.CODATA2018
using NaNMath

# =============================================================================
# EARTH PARAMETERS
# =============================================================================

"""Gravitational parameter μ = GM for Earth"""
const μ_Earth = 3.986004418e14u"m^3/s^2"

"""Earth equatorial radius"""
const R_Earth = 6.3781370e6u"m"

"""Earth's J2 oblateness coefficient (dimensionless)"""
const J2_Earth = 1.08263e-3

"""Standard gravitational acceleration at sea level"""
const g_0 = 9.80665u"m/s^2"

"""Earth's rotation rate"""
const ω_Earth = 7.2921159e-5u"rad/s"

# =============================================================================
# ATMOSPHERIC PARAMETERS (US Standard Atmosphere 1976)
# =============================================================================

"""Sea level atmospheric density"""
const ρ_0 = 1.225u"kg/m^3"

"""Sea level pressure"""
const P_0 = 101325.0u"Pa"

"""Sea level temperature"""
const T_0 = 288.15u"K"

"""Atmospheric scale height (simplified model)"""
const H_scale = 8500.0u"m"

"""Gas constant for air"""
const R_air = 287.05u"J/(kg*K)"

"""Specific heat ratio for air"""
const γ_air = 1.4  # dimensionless

"""Dynamic viscosity at sea level"""
const μ_air_0 = 1.789e-5u"Pa*s"

"""Sutherland's constant for air"""
const T_sutherland = 110.4u"K"

# =============================================================================
# EM LAUNCHER TARGET PARAMETERS
# =============================================================================

"""Target final velocity for orbital insertion (8 km/s)"""
const v_target_orbital = 8000.0u"m/s"

"""Target velocity for hypersonic testing (Mach 6)"""
const v_target_hypersonic = 2000.0u"m/s"  # ~Mach 6 at altitude

"""Target payload fraction (45% vs 4% for rockets)"""
const payload_fraction_target = 0.45  # dimensionless

"""Estimated launcher tube length range"""
const L_launcher_min = 1000.0u"m"
const L_launcher_max = 3000.0u"m"
const L_launcher_nominal = 2000.0u"m"

"""Average acceleration for nominal launcher"""
const a_avg_nominal = v_target_orbital^2 / (2 * L_launcher_nominal)

"""Maximum structural acceleration limit (conservative estimate)"""
const a_max_structural = 50000.0 * g_0  # 50,000 g

"""Payload mass range"""
const m_payload_min = 5.0u"kg"
const m_payload_max = 100.0u"kg"
const m_payload_nominal = 20.0u"kg"

# =============================================================================
# ELECTROMAGNETIC COIL PARAMETERS
# =============================================================================

"""Number of acceleration coils (estimate)"""
const N_coils_min = 100
const N_coils_max = 500
const N_coils_nominal = 200

"""Maximum coil current (very high for mass driver)"""
const I_coil_max = 100_000.0u"A"

"""Typical coil inductance per meter"""
const L_coil_per_m = 1.0e-6u"H/m"

"""Coil resistance (estimate)"""
const R_coil = 0.01u"Ω"

"""Capacitor bank capacity range"""
const C_capacitor_min = 1.0u"F"
const C_capacitor_max = 100.0u"F"

"""Capacitor voltage range"""
const V_capacitor_max = 10_000.0u"V"

"""Coil efficiency (electrical to kinetic energy)"""
const η_coil = 0.60  # 60% efficiency (dimensionless)

# =============================================================================
# THERMAL PARAMETERS
# =============================================================================

"""Aerodynamic heating coefficient (Detra-Kemp-Riddell)"""
const K_heating = 1.83e-4u"kg^0.5/m"

"""Stefan-Boltzmann constant (from PhysicalConstants)"""
const σ_SB = 5.670374419e-8u"W/(m^2*K^4)"

"""Payload nose radius (assumed)"""
const R_nose = 0.1u"m"

"""Payload thermal protection emissivity"""
const ε_thermal = 0.85  # dimensionless

"""Maximum allowable payload temperature"""
const T_max_payload = 1500.0u"K"

"""Ambient temperature (space)"""
const T_ambient = 250.0u"K"

"""Specific heat capacity of payload (aluminum-like)"""
const c_p_payload = 900.0u"J/(kg*K)"

# =============================================================================
# DRAG PARAMETERS
# =============================================================================

"""Payload reference area (assumed cylindrical)"""
const A_ref = π * (0.15u"m")^2  # 0.3m diameter cylinder

"""Drag coefficient (supersonic/hypersonic regime)"""
const C_d_supersonic = 0.85  # dimensionless

"""Drag coefficient (subsonic regime)"""
const C_d_subsonic = 0.40  # dimensionless

# =============================================================================
# MISSION PARAMETERS
# =============================================================================

"""Low Earth Orbit altitude target"""
const h_LEO = 400.0u"km"

"""Hypersonic test altitude"""
const h_hypersonic_test = 25.0u"km"

"""Launch site latitude (Caesarea, Israel)"""
const lat_launch = 32.5u"°"

"""Launch site longitude (Caesarea, Israel)"""
const lon_launch = 34.9u"°"

"""Launch site altitude (sea level)"""
const h_launch = 0.0u"m"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

"""
    coil_spacing(N_coils, L_launcher)

Calculate spacing between coils.

# Arguments
- `N_coils`: Number of coils (dimensionless)
- `L_launcher`: Total launcher length [m]

# Returns
- Coil spacing [m]
"""
function coil_spacing(N_coils::Integer, L_launcher)
    return L_launcher / N_coils
end

"""
    average_acceleration(v_final, L_launcher)

Calculate average acceleration needed to reach final velocity.

Uses kinematic equation: v² = 2aL

# Arguments
- `v_final`: Final velocity [m/s]
- `L_launcher`: Launcher length [m]

# Returns
- Average acceleration [m/s²]
"""
function average_acceleration(v_final, L_launcher)
    return v_final^2 / (2 * L_launcher)
end

"""
    g_loading(acceleration)

Convert acceleration to g-loading (multiple of standard gravity).

# Arguments
- `acceleration`: Acceleration [m/s²]

# Returns
- G-loading (dimensionless)
"""
function g_loading(acceleration)
    return ustrip(u"m/s^2", acceleration) / ustrip(u"m/s^2", g_0)
end

# Note: mach_number function moved to atmosphere.jl to avoid dispatch confusion
# atmosphere.jl has mach_number(velocity, altitude) which computes Mach from altitude
# This version that takes temperature directly is not used in the codebase

# =============================================================================
# DISPLAY FUNCTIONS
# =============================================================================

"""
    print_launcher_parameters()

Print summary of EM launcher parameters with proper units.
"""
function print_launcher_parameters()
    println("=" ^ 70)
    println("ELECTROMAGNETIC LAUNCHER - SYSTEM PARAMETERS")
    println("=" ^ 70)
    println()

    println("TARGET PERFORMANCE:")
    println("  Orbital velocity target:     $v_target_orbital")
    println("  Hypersonic test velocity:    $v_target_hypersonic")
    println("  Payload fraction target:     $(payload_fraction_target * 100)%")
    println()

    println("LAUNCHER GEOMETRY:")
    println("  Nominal tube length:         $L_launcher_nominal")
    println("  Tube length range:           $L_launcher_min - $L_launcher_max")
    println("  Nominal coil count:          $N_coils_nominal")
    println("  Coil spacing (nominal):      $(coil_spacing(N_coils_nominal, L_launcher_nominal))")
    println()

    println("ACCELERATION PROFILE:")
    println("  Average acceleration:        $a_avg_nominal")
    println("  G-loading (average):         $(round(g_loading(a_avg_nominal), digits=1)) g")
    println("  Maximum structural limit:    $a_max_structural")
    println("  G-loading (maximum):         $(round(g_loading(a_max_structural), digits=1)) g")
    println()

    println("PAYLOAD SPECIFICATIONS:")
    println("  Mass range:                  $m_payload_min - $m_payload_max")
    println("  Nominal mass:                $m_payload_nominal")
    println("  Reference area:              $A_ref")
    println("  Max temperature:             $T_max_payload")
    println()

    println("ELECTROMAGNETIC SYSTEM:")
    println("  Maximum coil current:        $I_coil_max")
    println("  Coil resistance:             $R_coil")
    println("  Inductance per meter:        $L_coil_per_m")
    println("  Maximum capacitor voltage:   $V_capacitor_max")
    println("  System efficiency:           $(η_coil * 100)%")
    println()

    println("MISSION PROFILES:")
    println("  LEO altitude target:         $h_LEO")
    println("  Hypersonic test altitude:    $h_hypersonic_test")
    println("  Launch site (Caesarea):      $lat_launch, $lon_launch")
    println()
    println("=" ^ 70)
end

export μ_Earth, R_Earth, J2_Earth, g_0, ω_Earth
export ρ_0, P_0, T_0, H_scale, R_air, γ_air, μ_air_0, T_sutherland
export v_target_orbital, v_target_hypersonic, payload_fraction_target
export L_launcher_min, L_launcher_max, L_launcher_nominal, a_avg_nominal, a_max_structural
export m_payload_min, m_payload_max, m_payload_nominal
export N_coils_min, N_coils_max, N_coils_nominal
export I_coil_max, L_coil_per_m, R_coil, C_capacitor_min, C_capacitor_max, V_capacitor_max, η_coil
export K_heating, σ_SB, R_nose, ε_thermal, T_max_payload, T_ambient, c_p_payload
export A_ref, C_d_supersonic, C_d_subsonic
export h_LEO, h_hypersonic_test, lat_launch, lon_launch, h_launch
export coil_spacing, average_acceleration, g_loading
export print_launcher_parameters
