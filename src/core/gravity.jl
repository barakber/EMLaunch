"""
# Gravity Module

Models gravitational acceleration including oblateness (J2 perturbation)
and other effects important for launch trajectory simulation.
Supports multiple celestial bodies via the `body` keyword argument.

## Point Mass Gravity

Simple inverse-square law:
```math
\\vec{a}_g = -\\frac{\\mu}{r^3} \\vec{r}
```

where:
- `μ = GM` = gravitational parameter [m³/s²]
- `r` = distance from body center [m]
- `r̂` = unit vector pointing from body center to object

## J2 Perturbation (Oblateness)

The J2 term accounts for oblateness:

```math
\\vec{a}_{J2} = \\frac{3\\mu J_2 R^2}{2r^4} \\left[
\\left(5\\frac{z^2}{r^2} - 1\\right)\\hat{x} +
\\left(5\\frac{z^2}{r^2} - 1\\right)\\hat{y} +
\\left(5\\frac{z^2}{r^2} - 3\\right)\\hat{z}
\\right]
```

## Total Gravitational Acceleration

```math
\\vec{a}_{total} = \\vec{a}_g + \\vec{a}_{J2}
```
"""

using Unitful
using LinearAlgebra
using StaticArrays

"""
    gravity_point_mass(position; body=EARTH)

Calculate gravitational acceleration using point mass model.

```math
\\vec{a}_g = -\\frac{\\mu}{r^3} \\vec{r}
```

# Arguments
- `position`: Position vector from body center [m] (3D vector or scalar radius)
- `body`: CelestialBody (default: EARTH)

# Returns
- Gravitational acceleration [m/s²] (3D vector or scalar)
"""
function gravity_point_mass(position; body=EARTH)
    μ = body.mu

    if isa(position, Number)
        # Scalar radius - return scalar acceleration magnitude
        r = position
        a = -μ / r^2
        return a
    else
        # Vector position - return vector acceleration
        r_vec = position
        r = norm(r_vec)
        a_vec = -μ / r^3 * r_vec
        return a_vec
    end
end

"""
    gravity_with_J2(position; body=EARTH)

Calculate gravitational acceleration including J2 oblateness perturbation.

# Arguments
- `position`: Position vector [x, y, z] from body center [m]
- `body`: CelestialBody (default: EARTH)

# Returns
- Total gravitational acceleration vector [m/s²]
"""
function gravity_with_J2(position; body=EARTH)
    μ = body.mu
    J2 = body.J2
    R = body.radius

    # Position components
    x, y, z = position[1], position[2], position[3]
    r = norm(position)

    # Point mass gravity
    a_point = -μ / r^3 * position

    # J2 perturbation
    factor = (3 * μ * J2 * R^2) / (2 * r^5)  # Divide by r^5, not r^4

    z2_r2 = (z / r)^2

    a_J2_x = factor * (5 * z2_r2 - 1) * x
    a_J2_y = factor * (5 * z2_r2 - 1) * y
    a_J2_z = factor * (5 * z2_r2 - 3) * z

    a_J2 = SVector(a_J2_x, a_J2_y, a_J2_z)

    return a_point + a_J2
end

"""
    altitude_from_position(position; body=EARTH)

Calculate altitude above body's surface from position vector.

# Arguments
- `position`: Position vector from body center [m]
- `body`: CelestialBody (default: EARTH)

# Returns
- Altitude above surface [m]
"""
function altitude_from_position(position; body=EARTH)
    R = body.radius

    if isa(position, Number)
        r = position
    else
        r = norm(position)
    end

    return r - R
end

"""
    position_from_altitude_angle(altitude, latitude, longitude; body=EARTH)

Calculate ECI position vector from altitude and geographic coordinates.

Simplified model assuming spherical body and ignoring rotation effects
for initial position calculation.

# Arguments
- `altitude`: Altitude above surface [m]
- `latitude`: Geodetic latitude [degrees or radians]
- `longitude`: Geodetic longitude [degrees or radians]
- `body`: CelestialBody (default: EARTH)

# Returns
- Position vector [x, y, z] in body-centered inertial frame [m]
"""
function position_from_altitude_angle(altitude, latitude, longitude; body=EARTH)
    R = body.radius

    # Convert to radians if necessary
    if unit(latitude) == u"°"
        lat_rad = ustrip(u"rad", latitude)
    else
        lat_rad = ustrip(latitude)
    end

    if unit(longitude) == u"°"
        lon_rad = ustrip(u"rad", longitude)
    else
        lon_rad = ustrip(longitude)
    end

    # Radius from body center
    r = R + altitude

    # Convert to Cartesian (simplified spherical body)
    x = r * cos(lat_rad) * cos(lon_rad)
    y = r * cos(lat_rad) * sin(lon_rad)
    z = r * sin(lat_rad)

    return SVector(x, y, z)
end

"""
    escape_velocity(altitude; body=EARTH)

Calculate escape velocity at given altitude.

```math
v_{esc} = \\sqrt{\\frac{2\\mu}{r}}
```

# Arguments
- `altitude`: Altitude above surface [m]
- `body`: CelestialBody (default: EARTH)

# Returns
- Escape velocity [m/s]
"""
function escape_velocity(altitude; body=EARTH)
    μ = body.mu
    R = body.radius

    r = R + altitude

    v_esc = sqrt(2 * μ / r)

    return v_esc
end

"""
    orbital_velocity(altitude; body=EARTH)

Calculate circular orbital velocity at given altitude.

```math
v_{orb} = \\sqrt{\\frac{\\mu}{r}}
```

# Arguments
- `altitude`: Altitude above surface [m]
- `body`: CelestialBody (default: EARTH)

# Returns
- Circular orbital velocity [m/s]
"""
function orbital_velocity(altitude; body=EARTH)
    μ = body.mu
    R = body.radius

    r = R + altitude

    v_orb = sqrt(μ / r)

    return v_orb
end

"""
    specific_orbital_energy(velocity, position; body=EARTH)

Calculate specific orbital energy (energy per unit mass).

```math
\\epsilon = \\frac{v^2}{2} - \\frac{\\mu}{r}
```

# Arguments
- `velocity`: Velocity magnitude or vector [m/s]
- `position`: Distance from body center or position vector [m]
- `body`: CelestialBody (default: EARTH)

# Returns
- Specific orbital energy [J/kg] or [m²/s²]
"""
function specific_orbital_energy(velocity, position; body=EARTH)
    μ = body.mu

    v = isa(velocity, Number) ? velocity : norm(velocity)
    r = isa(position, Number) ? position : norm(position)

    ε = 0.5 * v^2 - μ / r

    return ε
end

"""
    orbital_elements(position, velocity; body=EARTH)

Calculate classical orbital elements from state vectors.

Returns:
- `a`: Semi-major axis [m]
- `e`: Eccentricity (dimensionless)
- `i`: Inclination [rad]
- `Ω`: Longitude of ascending node [rad]
- `ω`: Argument of periapsis [rad]
- `ν`: True anomaly [rad]

# Arguments
- `position`: Position vector [m]
- `velocity`: Velocity vector [m/s]
- `body`: CelestialBody (default: EARTH)

# Returns
- Named tuple with orbital elements
"""
function orbital_elements(position, velocity; body=EARTH)
    μ = body.mu

    # Position and velocity magnitudes
    r_vec = position
    v_vec = velocity
    r = norm(r_vec)
    v = norm(v_vec)

    # Angular momentum vector
    h_vec = cross(r_vec, v_vec)
    h = norm(h_vec)

    # Node vector
    k_hat = SVector(0.0u"m", 0.0u"m", 1.0u"m")
    n_vec = cross(k_hat / 1.0u"m", h_vec)  # Remove units from k_hat for cross product
    n = norm(n_vec)

    # Eccentricity vector
    e_vec = cross(v_vec, h_vec) / μ - r_vec / r
    e = norm(e_vec)

    # Specific orbital energy
    ε = v^2 / 2 - μ / r

    # Semi-major axis
    if abs(ustrip(u"m^2/s^2", ε)) > 1e-10  # Not parabolic
        a = -μ / (2 * ε)
    else
        a = Inf * u"m"  # Parabolic orbit
    end

    # Inclination
    i = acos(ustrip(h_vec[3] / h)) * u"rad"

    # Longitude of ascending node
    if n > 0.0u"m^2/s"
        Ω = acos(ustrip(n_vec[1] / n)) * u"rad"
        if ustrip(n_vec[2]) < 0
            Ω = 2π * u"rad" - Ω
        end
    else
        Ω = 0.0u"rad"
    end

    # Argument of periapsis
    if n > 0.0u"m^2/s" && e > 1e-10
        ω = acos(ustrip(dot(n_vec, e_vec) / (n * e))) * u"rad"
        if ustrip(e_vec[3]) < 0
            ω = 2π * u"rad" - ω
        end
    else
        ω = 0.0u"rad"
    end

    # True anomaly
    if e > 1e-10
        ν = acos(ustrip(dot(e_vec, r_vec) / (e * r))) * u"rad"
        if ustrip(dot(r_vec, v_vec)) < 0
            ν = 2π * u"rad" - ν
        end
    else
        ν = 0.0u"rad"
    end

    return (
        a = a,          # Semi-major axis
        e = e,          # Eccentricity
        i = i,          # Inclination
        Ω = Ω,          # Longitude of ascending node
        ω = ω,          # Argument of periapsis
        ν = ν,          # True anomaly
        h = h           # Angular momentum magnitude
    )
end

"""
    print_orbital_info(position, velocity; body=EARTH)

Print human-readable orbital information.

# Arguments
- `position`: Position vector [m]
- `velocity`: Velocity vector [m/s]
- `body`: CelestialBody (default: EARTH)
"""
function print_orbital_info(position, velocity; body=EARTH)
    elements = orbital_elements(position, velocity; body=body)
    r = norm(position)
    v = norm(velocity)
    h = altitude_from_position(position; body=body)
    ε = specific_orbital_energy(velocity, position; body=body)

    println("=" ^ 70)
    println("ORBITAL INFORMATION ($(body.name))")
    println("=" ^ 70)
    println()
    println("STATE VECTORS:")
    println("  Radius:                $(r)")
    println("  Altitude:              $(h)")
    println("  Velocity:              $(v)")
    println("  Specific energy:       $(ε)")
    println()
    println("ORBITAL ELEMENTS:")
    println("  Semi-major axis (a):   $(elements.a)")
    println("  Eccentricity (e):      $(elements.e)")
    println("  Inclination (i):       $(uconvert(u"°", elements.i))")
    println("  Long. asc. node (Ω):   $(uconvert(u"°", elements.Ω))")
    println("  Arg. periapsis (ω):    $(uconvert(u"°", elements.ω))")
    println("  True anomaly (ν):      $(uconvert(u"°", elements.ν))")
    println()

    # Orbit type
    R = body.radius
    if elements.e < 1e-3
        println("  Orbit type:            Circular")
    elseif elements.e < 1.0
        println("  Orbit type:            Elliptical")
        r_p = elements.a * (1 - elements.e)
        r_a = elements.a * (1 + elements.e)
        h_p = r_p - R
        h_a = r_a - R
        println("  Periapsis altitude:    $(h_p)")
        println("  Apoapsis altitude:     $(h_a)")
    elseif elements.e ≈ 1.0
        println("  Orbit type:            Parabolic (escape)")
    else
        println("  Orbit type:            Hyperbolic (escape)")
    end
    println("=" ^ 70)
end

export gravity_point_mass, gravity_with_J2
export altitude_from_position, position_from_altitude_angle
export escape_velocity, orbital_velocity, specific_orbital_energy
export orbital_elements, print_orbital_info
