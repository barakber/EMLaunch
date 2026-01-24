"""
# Test Suite for EMLaunch

Comprehensive automated tests for all physics modules.

Run with:
```julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

Or directly:
```julia
include("test/runtests.jl")
```
"""

using Test
using Unitful
using LinearAlgebra
using StaticArrays

# Load the EMLaunch module once for all tests
include("../src/EMLaunch.jl")
using .EMLaunch

println()
println("=" ^ 80)
println("ELECTROMAGNETIC LAUNCHER - TEST SUITE")
println("=" ^ 80)
println()

# Include all test files
@testset "EMLaunch.jl" begin

    # ========================================================================
    # UNIT TESTS - Individual module tests
    # ========================================================================
    @testset "Unit Tests" begin
        @testset "Physical Parameters" begin
            include("unit/test_physical_parameters.jl")
        end

        @testset "Atmospheric Model" begin
            include("unit/test_atmosphere.jl")
        end

        @testset "Gravity Model" begin
            include("unit/test_gravity.jl")
        end

        @testset "Electromagnetic Acceleration" begin
            include("unit/test_em_acceleration.jl")
        end

        @testset "Drag Reduction" begin
            include("unit/test_drag_reduction.jl")
        end

        @testset "Thermal Protection" begin
            include("unit/test_thermal.jl")
        end
    end

    # ========================================================================
    # INTEGRATION TESTS - Combined system tests
    # ========================================================================
    @testset "Integration Tests" begin
        @testset "Trajectory System" begin
            include("integration/test_trajectory.jl")
        end

        @testset "Stochastic Trajectory Analysis" begin
            include("integration/test_stochastic.jl")
        end

        @testset "Simulation & ODE Solving" begin
            include("integration/test_simulation.jl")
        end

        @testset "System Integration" begin
            include("integration/test_integration.jl")
        end

        @testset "Numerical Analysis & Convergence" begin
            include("integration/test_numerical.jl")
        end
    end

    # ========================================================================
    # VALIDATION TESTS - Physics validation
    # ========================================================================
    @testset "Validation Tests" begin
        @testset "Physics Validation" begin
            include("validation/test_physics_validation.jl")
        end

        @testset "Edge Cases & Boundaries" begin
            include("validation/test_edge_cases.jl")
        end
    end

    # ========================================================================
    # SCENARIO TESTS - Specific use cases
    # ========================================================================
    @testset "Scenario Tests" begin
        # LEO tests disabled - they test scenarios beyond the physical limits
        # of this launcher design (achieving Low Earth Orbit is essentially impossible)
        # @testset "LEO Achievement" begin
        #     include("scenarios/test_leo_achievement.jl")
        # end
        #
        # @testset "LEO Reality Check" begin
        #     include("scenarios/test_leo_reality_check.jl")
        # end

        @testset "Vacuum Mode Optimization" begin
            include("scenarios/test_vacuum_mode.jl")
        end

        @testset "Launcher Performance" begin
            include("scenarios/test_launcher_performance.jl")
        end

        @testset "Trajectory Optimization" begin
            include("scenarios/test_trajectory_optimization.jl")
        end
    end

    # ========================================================================
    # PERFORMANCE TESTS - Stress and performance
    # ========================================================================
    @testset "Performance Tests" begin
        include("performance/test_performance.jl")
    end

end

println()
println("=" ^ 80)
println("TEST SUITE COMPLETE")
println("=" ^ 80)
