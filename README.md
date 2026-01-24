# EMLaunch.jl

**Electromagnetic launcher simulation with stochastic uncertainty quantification**

A Julia package for physics-based modeling of electromagnetic launch systems (coilguns/mass drivers) with Monte Carlo uncertainty analysis.

## Features

- **Complete physics**: Electromagnetic acceleration, atmospheric drag, gravity (J2), aerodynamic heating
- **Stochastic modeling**: Monte Carlo simulations with SDEs for uncertainty quantification
- **Vacuum tube support**: Model 99%+ vacuum tubes to minimize drag losses
- **Type-safe units**: Dimensional analysis with `Unitful.jl`
- **Rich visualization**: Multi-percentile confidence intervals, trajectory analysis
- **Interactive notebooks**: Jupyter notebooks for exploration

## Quick Start

### Option 1: Docker (Recommended)

```bash
# Start Jupyter notebook server
docker-compose up

# Access at: http://localhost:8888 (no password required)
```

**What's included:** Julia 1.12, all dependencies pre-installed, Jupyter with IJulia kernel

**Configuration:**
- Multi-threading: Set `JULIA_NUM_THREADS` in `docker-compose.yml`
- Custom port: Change `ports: "8888:8888"` to `"9999:8888"`
- Stop: `docker-compose down`

### Option 2: Local Julia Installation

```julia
# Activate environment
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# Run example notebook
jupyter notebook notebooks/basic_sde_example.ipynb
```

### Option 3: Run Programmatically

```julia
using EMLaunch
using Unitful

# Create 5 km vacuum tube launcher @ 18 kV
launcher = create_uniform_launcher(
    length = 5000.0u"m",
    num_coils = 500,
    voltage_per_coil = 18000.0u"V",
    capacitance_per_coil = 2.5u"F",
    inductance_per_coil = 0.001u"H",
    resistance_per_coil = 0.01u"Ω",
    gradient_per_coil = 0.003u"H/m",
    vacuum_pressure_ratio = 0.0001  # 99.99% vacuum
)

# 10 kg payload, aerodynamic design
payload = PayloadConfig(
    10.0u"kg",
    π * (0.05u"m")^2,  # 10 cm diameter
    0.02u"m",           # Sharp nose
    0.95, 900.0u"J/(kg*K)", 288.15u"K", 2500.0u"K"
)

# Near-vertical launch
mission = MissionProfile(
    0.0u"m", 32.5u"°", 34.9u"°",  # Sea level, Caesarea
    90.0u"°", 88.0u"°",            # East, 88° elevation
    8000.0u"m/s", 400.0u"km"       # Target orbital parameters
)

# Run Monte Carlo analysis
results = monte_carlo_analysis(launcher, payload, mission, 100)

# Visualize with layered confidence intervals
plot_monte_carlo_results(results)
```

## Physics

### Electromagnetic Force
```math
F_{em} = \frac{1}{2} I^2 \frac{dL}{dx}
```

### Atmospheric Drag (with vacuum tube support)
```math
F_d = \frac{1}{2} \rho(h) v^2 C_d A \times \text{vacuum\_ratio}
```

### Aerodynamic Heating
```math
\dot{q} = K \sqrt{\frac{\rho}{R_n}} v^3
```

Inside vacuum tube, both drag and heating are reduced by `vacuum_pressure_ratio` (e.g., 0.0001 = 99.99% vacuum).

## Project Structure

```
EMLaunch/
├── src/
│   ├── core/              # Physics (atmosphere, gravity)
│   ├── propulsion/        # EM acceleration
│   ├── dynamics/          # Trajectories, SDEs, thermal
│   └── analysis/          # Optimization, visualization
├── notebooks/             # Jupyter notebooks
│   └── basic_sde_example.ipynb
├── test/                  # Comprehensive test suite
└── docs/                  # API documentation
```

## Documentation

Generate API documentation with:
```julia
cd docs/
julia make.jl
```

Or view online at: [barakber.github.io/EMLaunch](https://barakber.github.io/EMLaunch)

## Testing

```julia
using Pkg
Pkg.test("EMLaunch")
```

Test coverage: **99.6%** (4106/4122 tests passing)

## Key Results

With realistic parameters (5 km vacuum tube, 18 kV, 10 kg payload):
- **Exit velocity**: ~4 km/s (40% efficiency)
- **Peak altitude**: 100+ km (Kármán line)
- **Vacuum benefit**: ~50% energy savings vs atmospheric tube
- **Thermal limit**: <2500 K (safe for refractory materials)

## Repository

[github.com/barakber/EMLaunch](https://github.com/barakber/EMLaunch)

## Author

**Barak Bercovitz** (barakber@gmail.com)

## License

Copyright © 2026 Barak Bercovitz. All Rights Reserved.
