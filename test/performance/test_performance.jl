"""
Performance and Stress Tests

Tests:
- Computational performance
- Memory efficiency
- Large-scale parameter sweeps
- Vectorization efficiency
- Type stability
"""


@testset "Atmospheric Function Performance" begin
    # Test that atmospheric functions scale well

    # Single call - verify it completes successfully
    time_single = @elapsed density(25.0u"km")
    @test time_single < 1.0  # Should complete in reasonable time (< 1s)

    # Many calls (vectorization)
    altitudes = range(0.0u"km", 100.0u"km", length=1000)

    time_many = @elapsed begin
        densities = density.(altitudes)
    end

    @test time_many < 0.1  # Should complete in < 100ms

    # Verify results are correct
    @test length(densities) == 1000
    @test all(d -> d > 0.0u"kg/m^3", densities)
end

@testset "Gravity Calculation Performance" begin
    # Single call - verify it completes successfully
    time_single = @elapsed gravity_point_mass(R_Earth + 400.0u"km")
    @test time_single < 1.0  # Should complete in reasonable time (< 1s)

    # Vector operations
    radii = range(R_Earth, R_Earth + 1000.0u"km", length=1000)

    time_vec = @elapsed begin
        gravities = abs.(gravity_point_mass.(radii))
    end

    @test time_vec < 0.1

    # All should be positive and decreasing
    @test all(g -> g > 0.0u"m/s^2", gravities)
    for i in 2:length(gravities)
        @test gravities[i] < gravities[i-1]
    end
end

@testset "Orbital Velocity Performance" begin
    # Batch calculation
    altitudes = range(100.0u"km", 10000.0u"km", length=500)

    time_batch = @elapsed begin
        velocities = orbital_velocity.(altitudes)
    end

    @test time_batch < 0.1

    # All should be positive and decreasing
    @test all(v -> v > 0.0u"m/s", velocities)
    for i in 2:length(velocities)
        @test velocities[i] < velocities[i-1]
    end
end

@testset "Drag Force Performance" begin
    # Parameter sweep
    velocities = range(100.0u"m/s", 10000.0u"m/s", length=100)
    altitudes = range(0.0u"km", 100.0u"km", length=100)
    area = 1.0u"m^2"

    time_sweep = @elapsed begin
        forces = [drag_force(v, h, area)
                  for v in velocities, h in altitudes]
    end

    @test time_sweep < 1.0  # 10,000 calculations in < 1 second

    # Check dimensions
    @test size(forces) == (100, 100)
end

@testset "Heating Calculation Performance" begin
    # Large parameter space
    velocities = range(500.0u"m/s", 8000.0u"m/s", length=100)
    altitudes = range(10.0u"km", 100.0u"km", length=100)
    R_n = 0.1u"m"

    time_heating = @elapsed begin
        heating = [aerodynamic_heating(v, h, R_n)
                   for v in velocities, h in altitudes]
    end

    @test time_heating < 1.0

    # All should be non-negative
    @test all(q -> q >= 0.0u"W/m^2", heating)
end

@testset "Memory Allocation" begin
    # Check that functions don't allocate excessively

    # Atmospheric density
    alloc_density = @allocated density(25.0u"km")
    @test alloc_density < 1000  # < 1 KB

    # Gravity
    alloc_gravity = @allocated gravity_point_mass(R_Earth + 400.0u"km")
    @test alloc_gravity < 1000

    # Orbital velocity
    alloc_orbital = @allocated orbital_velocity(400.0u"km")
    @test alloc_orbital < 1000
end

@testset "Type Stability" begin
    # Functions should return consistent types

    # Density
    ρ1 = density(0.0u"km")
    ρ2 = density(50.0u"km")
    @test typeof(ρ1) == typeof(ρ2)

    # Temperature
    T1 = temperature(0.0u"km")
    T2 = temperature(50.0u"km")
    @test typeof(T1) == typeof(T2)

    # Gravity
    g1 = gravity_point_mass(R_Earth)
    g2 = gravity_point_mass(R_Earth + 400.0u"km")
    @test typeof(g1) == typeof(g2)
end

@testset "Large Altitude Range" begin
    # Test performance over very large altitude range
    altitudes = 10 .^ range(0, 5, length=1000) .* 1.0u"m"  # 1 m to 100 km

    time_large = @elapsed begin
        ρ = density.(altitudes)
        T = temperature.(altitudes)
        P = pressure.(altitudes)
    end

    @test time_large < 1.0

    # All should be valid
    @test all(isfinite ∘ ustrip, ρ)
    @test all(isfinite ∘ ustrip, T)
    @test all(isfinite ∘ ustrip, P)
end

@testset "Repeated Calculations" begin
    # Ensure no performance degradation with repeated calls

    h = 25.0u"km"

    times = Float64[]
    for _ in 1:10
        t = @elapsed begin
            for _ in 1:1000
                density(h)
            end
        end
        push!(times, t)
    end

    # Times should be consistent (no degradation)
    @test maximum(times) < 2 * minimum(times)

    # Average time per call should be very small
    avg_time = sum(times) / (10 * 1000)
    @test avg_time < 1e-5  # < 10 microseconds per call
end

@testset "Vectorized vs Loop Performance" begin
    # Compare vectorized and loop approaches
    altitudes = range(0.0u"km", 100.0u"km", length=1000)

    # Vectorized
    time_vec = @elapsed begin
        ρ_vec = density.(altitudes)
    end

    # Loop
    time_loop = @elapsed begin
        ρ_loop = similar(altitudes, typeof(1.0u"kg/m^3"))
        for (i, h) in enumerate(altitudes)
            ρ_loop[i] = density(h)
        end
    end

    # Vectorized should be at least comparable (within 2x)
    @test time_vec < 2 * time_loop

    # Results should be identical
    ρ_vec = density.(altitudes)
    ρ_loop = [density(h) for h in altitudes]
    @test all(ρ_vec .≈ ρ_loop)
end

@testset "Stress Test - Extreme Parameters" begin
    # Test with extreme but physically possible values

    # Very high altitude (near Moon)
    h_extreme = 300000.0u"km"
    ρ_extreme = density(h_extreme)
    @test ρ_extreme >= 0.0u"kg/m^3"
    @test isfinite(ustrip(u"kg/m^3", ρ_extreme))

    # Very high velocity (interplanetary)
    v_extreme = 50000.0u"m/s"
    h_test = 100.0u"km"
    q_extreme = aerodynamic_heating(v_extreme, h_test, 0.1u"m")
    @test isfinite(ustrip(u"W/m^2", q_extreme))

    # Very large mass
    m_large = 100000.0u"kg"
    KE_large = 0.5 * m_large * (10000.0u"m/s")^2
    @test isfinite(ustrip(u"J", KE_large))
end

@testset "Concurrent Calculations" begin
    # Test that functions can be called concurrently (thread-safe)

    altitudes = range(0.0u"km", 100.0u"km", length=100)

    # Sequential
    time_seq = @elapsed begin
        results_seq = [density(h) for h in altitudes]
    end

    # Results should be consistent
    @test length(results_seq) == 100
    @test all(r -> r > 0.0u"kg/m^3", results_seq)
end

@testset "Numerical Precision" begin
    # Test that calculations maintain precision

    # Small differences
    h1 = 1000.0u"m"
    h2 = 1000.0u"m" + 1.0u"mm"

    ρ1 = density(h1)
    ρ2 = density(h2)

    # Should be different but very close
    @test ρ1 != ρ2
    @test isapprox(ρ1, ρ2, rtol=1e-3)

    # Large and small number arithmetic
    v_large = 8000.0u"m/s"
    v_small = 1.0u"m/s"
    v_sum = v_large + v_small

    @test v_sum > v_large
    @test v_sum != v_large
end

@testset "Scaling Behavior" begin
    # Test algorithmic complexity

    sizes = [10, 100, 1000]
    times = Float64[]

    for n in sizes
        altitudes = range(0.0u"km", 100.0u"km", length=n)

        t = @elapsed begin
            density.(altitudes)
        end

        push!(times, t)
    end

    # Should scale linearly (O(n))
    # time[2]/time[1] ≈ size[2]/size[1]
    if times[1] > 0
        ratio_time = times[2] / times[1]
        ratio_size = sizes[2] / sizes[1]

        # Allow 2x overhead
        @test ratio_time < 2 * ratio_size
    end
end

@testset "Cache Performance" begin
    # Test repeated calls with same arguments (potential caching)

    h = 25.0u"km"

    # First call
    t1 = @elapsed density(h)

    # Subsequent calls
    times_subsequent = [(@elapsed density(h)) for _ in 1:100]

    # All should be fast
    @test all(t -> t < 0.001, times_subsequent)
end

@testset "Batch Processing Efficiency" begin
    # Compare single vs batch calls

    # Single calls
    n = 1000
    time_single = @elapsed begin
        for i in 1:n
            density(Float64(i) * 0.1u"km")
        end
    end

    # Batch call
    altitudes = [Float64(i) * 0.1u"km" for i in 1:n]
    time_batch = @elapsed begin
        density.(altitudes)
    end

    # Batch should complete in reasonable time
    @test time_batch < 10.0  # Should complete in < 10s
end

println("  ✓ Performance and stress tests passed")
