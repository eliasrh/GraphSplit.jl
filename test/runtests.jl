using Test
using LinearAlgebra
using Random
using TOML

include(joinpath(@__DIR__, "..", "src", "GraphSplit.jl"))
using .GraphSplit

module RunGraphSplitCLI
include(joinpath(@__DIR__, "..", "run_graphsplit.jl"))
end

module BuildTravelTimesCLI
include(joinpath(@__DIR__, "..", "build_travel_times.jl"))
end

module BenchmarkComparator
include(joinpath(@__DIR__, "..", "benchmark", "yifan2025", "compare_catalogs.jl"))
end

function config_leaf_paths(node::AbstractDict, prefix::String="")
    result = String[]
    for (key, value) in node
        path = isempty(prefix) ? String(key) : string(prefix, '.', key)
        if value isa AbstractDict
            append!(result, config_leaf_paths(value, path))
        else
            push!(result, path)
        end
    end
    return result
end

@testset "GraphSplit" begin
    @testset "command-line entry points" begin
        @test RunGraphSplitCLI.main(["--help"]) == 0
        @test BuildTravelTimesCLI.main(["--help"]) == 0

        mktempdir() do directory
            catalog_path = joinpath(directory, "catalog_dd.txt")
            truth_path = joinpath(directory, "truelocs.txt")
            prefix = joinpath(directory, "comparison")
            open(catalog_path, "w") do io
                println(io, "0 0 0 0 0 0 35.4000 -117.9560 4.0 1")
                println(io, "0 0 0 0 0 0 35.4050 -117.9510 4.5 2")
            end
            open(truth_path, "w") do io
                println(io, "35.4000 -117.9560 4.0")
                println(io, "35.4050 -117.9510 4.5")
            end
            metrics = BenchmarkComparator.main([catalog_path, truth_path, prefix])
            @test metrics.neighboring_pairs == 1
            @test metrics.mean_horizontal_accuracy_km <= 1.0e-10
            @test metrics.mean_depth_accuracy_km == 0.0
            @test metrics.mean_horizontal_precision_km <= 1.0e-10
            @test metrics.mean_depth_precision_km == 0.0
            @test metrics.chamfer_distance_km <= 1.0e-10
            @test isfile(prefix * "_metrics.csv")
            @test isfile(prefix * "_event_errors.csv")
        end
    end

    @testset "configuration is documented and unambiguous" begin
        root = joinpath(@__DIR__, "..")
        complete = TOML.parsefile(joinpath(root, "config", "graphsplit_complete.toml"))
        reference = read(joinpath(root, "docs", "CONFIGURATION_REFERENCE.md"), String)
        for path in config_leaf_paths(complete)
            @test occursin("`$path`", reference)
        end

        minimal = TOML.parsefile(joinpath(root, "config", "graphsplit_template.toml"))
        @test minimal["prelocation"]["damping_lambda"] == 1.0e-3
        @test minimal["relocation"]["damping_lambda"] == 5.0e-4

        cfg = GraphSplit.default_config()
        cfg["graph"]["augmentation"]["enabled"] = true
        @test_throws ErrorException GraphSplit.validate_config(cfg)
        cfg["graph"]["augmentation"]["add_edges"] = 10
        @test GraphSplit.validate_config(cfg) === cfg
        cfg["graph"]["augmentation"]["add_edges_per_event"] = 0.2
        @test_throws ErrorException GraphSplit.validate_config(cfg)
        cfg["graph"]["augmentation"]["add_edges"] = 0
        @test GraphSplit.validate_config(cfg) === cfg

        cfg["relocation"]["damping_lambda"] = -1.0
        @test_throws ErrorException GraphSplit.validate_config(cfg)
    end

    @testset "coordinates" begin
        lat = [35.9, 36.0, 36.1]
        lon = [-117.8, -117.7, -117.6]
        x, y = GraphSplit.local_xy(lat, lon, 36.0, -117.7, 6_371_000.0)
        lat2, lon2 = GraphSplit.local_xy_to_ll(x, y, 36.0, -117.7, 6_371_000.0)
        @test lat2 ≈ lat atol=1.0e-12
        @test lon2 ≈ lon atol=1.0e-12

        args = (15_000.0, -8_000.0, -4_000.0, 3_000.0, 36.0, -117.7, 6_371_000.0)
        _, gx, gy = GraphSplit.spherical_delta_gradient(args...)
        h = 0.1
        fxp = GraphSplit.spherical_delta_gradient(args[1] + h, args[2:end]...)[1]
        fxm = GraphSplit.spherical_delta_gradient(args[1] - h, args[2:end]...)[1]
        fyp = GraphSplit.spherical_delta_gradient(args[1], args[2] + h, args[3:end]...)[1]
        fym = GraphSplit.spherical_delta_gradient(args[1], args[2] - h, args[3:end]...)[1]
        @test gx ≈ (fxp - fxm) / (2h) rtol=1.0e-5 atol=1.0e-12
        @test gy ≈ (fyp - fym) / (2h) rtol=1.0e-5 atol=1.0e-12
    end

    @testset "exact k-d tree" begin
        rng = MersenneTwister(1234)
        points = randn(rng, 50, 3)
        tree = GraphSplit.KDTree(points)
        for query in 1:50
            indices, distances = GraphSplit.knn(tree, query, 7)
            brute = sort([(sum(abs2, @view(points[query, :]) .- @view(points[j, :])), j)
                for j in 1:50 if j != query])
            expected = Int32[item[2] for item in brute[1:7]]
            @test indices == expected
            @test distances ≈ sqrt.([item[1] for item in brute[1:7]]) rtol=1.0e-13
        end
    end

    @testset "constant travel-time gradients" begin
        radius = 6_371_000.0
        flat = GraphSplit.ConstantTravelTime(:cartesian, 6000.0, 3500.0, radius, 36.0, -117.7, radius)
        value, gx, gy, gz = GraphSplit.travel_time_gradient(flat, UInt8(1), 3000.0, 4000.0, 2000.0, 0.0, 0.0, 0.0)
        distance = sqrt(3000.0^2 + 4000.0^2 + 2000.0^2)
        @test value ≈ distance / 6000.0
        @test [gx, gy, gz] ≈ [3000.0, 4000.0, 2000.0] ./ (6000.0 * distance)

        radial = GraphSplit.ConstantTravelTime(:radial, 6000.0, 3500.0, radius, 36.0, -117.7, radius)
        point = [20_000.0, 15_000.0, 4000.0]
        station = [-5_000.0, 8_000.0, -500.0]
        value, gx, gy, gz = GraphSplit.travel_time_gradient(radial, UInt8(1), point..., station...)
        h = 0.1
        f(offset) = GraphSplit.travel_time_gradient(radial, UInt8(1),
            point[1] + offset[1], point[2] + offset[2], point[3] + offset[3], station...)[1]
        numeric = [(f(ntuple(k -> k == axis ? h : 0.0, 3)) -
            f(ntuple(k -> k == axis ? -h : 0.0, 3))) / (2h) for axis in 1:3]
        @test [gx, gy, gz] ≈ numeric rtol=2.0e-5 atol=1.0e-10
        @test value > 0.0
    end

    @testset "fast sweeping" begin
        z = collect(-1000.0:250.0:4000.0)
        q = collect(0.0:250.0:5000.0)
        velocity = fill(5000.0, length(z))
        source = findfirst(==(1000.0), z)
        field, _, converged = GraphSplit.fast_sweep_field(velocity, z, q, source,
            :cartesian, 6_371_000.0, 32, 1.0e-10)
        @test converged
        @test field[source, :] ≈ q ./ 5000.0 atol=2.0e-7
        @test field[:, 1] ≈ abs.(z .- z[source]) ./ 5000.0 atol=2.0e-7

        delta = q ./ 6_371_000.0
        radial_field, _, radial_converged = GraphSplit.fast_sweep_field(velocity, z, delta, source,
            :radial, 6_371_000.0, 32, 1.0e-10)
        @test radial_converged
        @test radial_field[source, :] ≈ ((6_371_000.0 - z[source]) .* delta) ./ 5000.0 atol=2.0e-7
    end

    @testset "native travel-time table" begin
        mktempdir() do directory
            path = joinpath(directory, "tiny.gstt")
            q = [0.0, 1000.0]; z = [0.0, 1000.0]; zs = [-1000.0, 1000.0]
            table = GraphSplit.create_empty_table(path, :cartesian, 6_371_000.0, 6_371_000.0,
                zeros(UInt8, 32), q, z, zs, 36.0, -117.7, true)
            for is in eachindex(zs), iz in eachindex(z), iq in eachindex(q)
                value = q[iq] / 6000.0 + z[iz] / 10_000.0 + zs[is] / 20_000.0
                table.p[iq, iz, is] = Float32(value)
                table.s[iq, iz, is] = Float32(2value)
            end
            GraphSplit.Mmap.sync!(table.p); GraphSplit.Mmap.sync!(table.s)
            GraphSplit.mark_table_complete!(table)
            close(table.io)

            loaded = GraphSplit.load_travel_time_table(path, 36.0, -117.7, true)
            @test loaded.geometry == :cartesian
            @test loaded.model_hash == zeros(UInt8, 32)
            value, gx, gy, gz = GraphSplit.travel_time_gradient(loaded, UInt8(1),
                200.0, 0.0, 300.0, 0.0, 0.0, 400.0)
            @test value ≈ 200.0 / 6000.0 + 300.0 / 10_000.0 + 400.0 / 20_000.0 rtol=2.0e-6
            @test gx ≈ 1.0 / 6000.0 rtol=2.0e-6
            @test gy == 0.0
            @test gz ≈ 1.0 / 10_000.0 rtol=2.0e-6
            close(loaded.io)
        end
    end

    @testset "serial-ID theta joins" begin
        raw = zeros(3, 10)
        raw[:, 7] .= [36.0, 36.1, 36.2]; raw[:, 8] .= [-117.8, -117.7, -117.6]
        raw[:, 9] .= [4.0, 5.0, 6.0]; raw[:, 10] .= [30, 10, 20]
        catalog = GraphSplit.Catalog("memory", raw, Int64[30, 10, 20], raw[:, 7], raw[:, 8], raw[:, 9],
            7, 8, 9, 10, 36.0, -117.7, 6_371_000.0)
        stations = GraphSplit.Stations(["STA"], [36.0], [-117.7], [0.0], [0.0], [0.0], [0.0],
            36.0, -117.7, 6_371_000.0)
        entries = Dict{Int64,GraphSplit.ThetaEntry}(
            30 => GraphSplit.ThetaEntry(0.3, 20, 0.01, 10),
            10 => GraphSplit.ThetaEntry(0.1, 20, 0.01, 10),
            20 => GraphSplit.ThetaEntry(0.0, 20, 0.01, 10),
        )
        groups = [GraphSplit.ThetaGroup("STA", UInt8(1), "theta_STA_P.txt", Int64[30, 10, 20], entries, true, true)]
        cfg = GraphSplit.default_config()
        cfg["observations"]["minimum_theta_degree"] = 0
        cfg["observations"]["minimum_observations_per_pair"] = 1
        graph = GraphSplit.EventGraph(Int32[1], Int32[2], [1.0], [1.0], [1.0], Int32[1, 1, 0],
            Int32[1, 1, 2], Int32[2, 2, 1], 2, falses(1))
        star = GraphSplit.build_star_observations(groups, stations, catalog, cfg)
        @test Set(zip(star.i, star.j)) == Set([(Int32(1), Int32(3)), (Int32(2), Int32(3))])
        dd = GraphSplit.build_dd_observations(groups, stations, catalog, graph, cfg)
        @test length(dd) == 1
        @test (dd.i[1], dd.j[1]) == (Int32(1), Int32(2))
        @test dd.dt[1] ≈ 0.2
    end

    @testset "matrix-free adjoint and PCG" begin
        system = GraphSplit.LinearizedSystem(Int32[1, 2], Int32[2, 3],
            [0.1, -0.2], [0.2, 0.1], [0.3, 0.4], [-0.1, 0.25], [0.05, -0.3],
            [0.2, 0.15], [1.2, 0.8], 3, :xyzt0_scaled, 2.0, 5000.0)
        rng = MersenneTwister(9)
        x = randn(rng, 12); y = randn(rng, 6)
        @test dot(GraphSplit.apply_A(system, x), y) ≈ dot(x, GraphSplit.apply_At(system, y)) rtol=1.0e-13

        matrix = [4.0 1.0; 1.0 3.0]
        truth = [1.0, -2.0]
        solution, iterations, converged = GraphSplit.pcg(v -> matrix * v, matrix * truth, copy;
            tolerance=1.0e-12, maximum_iterations=10)
        @test converged
        @test iterations <= 2
        @test solution ≈ truth atol=1.0e-12
    end

    @testset "bundled benchmark inputs" begin
        input = joinpath(@__DIR__, "..", "benchmark", "yifan2025", "input")
        cfg = GraphSplit.default_config()
        catalog = GraphSplit.read_catalog(joinpath(input, "catalog.txt"), cfg)
        stations = GraphSplit.read_stations(joinpath(input, "stations.txt"))
        model = GraphSplit.read_velocity_model(joinpath(input, "vm.txt"))
        @test length(catalog) == 1000
        @test catalog.event_id == Int64.(1:1000)
        @test length(stations) == 27
        @test length(model.depth_m) == 10
    end
end
