# EMLaunch.jl Documentation

**Electromagnetic launcher simulation with stochastic uncertainty quantification**

A Julia package for physics-based modeling of electromagnetic launch systems (coilguns/mass drivers) with Monte Carlo uncertainty analysis.

## Quick Start

See the [README](https://github.com/barakber/EMLaunch) for installation and quick start guide.

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

# Run Monte Carlo analysis
results = monte_carlo_analysis(launcher, payload, mission, 100)
plot_monte_carlo_results(results)
```

## API Reference

```@autodocs
Modules = [EMLaunch]
Order = [:type, :function]
```

## Index

```@index
```
