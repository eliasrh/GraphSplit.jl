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

        cfg = GraphSplit.default_config()
        cfg["uncertainty"]["method"] = "both"
        @test GraphSplit.validate_config(cfg) === cfg
        cfg["uncertainty"]["bootstrap"]["resampling_unit"] = "event"
        @test_throws ErrorException GraphSplit.validate_config(cfg)
        cfg["uncertainty"]["bootstrap"]["resampling_unit"] = "station_phase"
        cfg["uncertainty"]["bootstrap"]["confidence_level"] = 1.0
        @test_throws ErrorException GraphSplit.validate_config(cfg)

        cfg = GraphSplit.default_config()
        cfg["initialization"]["mode"] = "common_centroid"
        @test GraphSplit.validate_config(cfg) === cfg
        cfg["run"]["prelocation"] = false
        @test_throws ErrorException GraphSplit.validate_config(cfg)
        cfg["io"]["restart_shift_file"] = "pass1/catalog_dd_dxdydzt0.txt"
        @test_throws ErrorException GraphSplit.validate_config(cfg)
        cfg["initialization"]["mode"] = "catalog"
        @test GraphSplit.validate_config(cfg) === cfg
        cfg["run"]["prelocation"] = true
        @test_throws ErrorException GraphSplit.validate_config(cfg)
        cfg["io"]["restart_shift_file"] = ""
        cfg["run"]["prelocation"] = true
        cfg["initialization"]["mode"] = "common_manual"
        cfg["initialization"]["latitude"] = 91.0
        @test_throws ErrorException GraphSplit.validate_config(cfg)
        cfg["initialization"]["latitude"] = 35.9
        @test GraphSplit.validate_config(cfg) === cfg
        cfg["initialization"]["mode"] = "one_event"
        @test_throws ErrorException GraphSplit.validate_config(cfg)

        cfg = GraphSplit.default_config()
        cfg["constraints"]["depth_bound"]["enabled"] = true
        cfg["constraints"]["depth_bound"]["minimum_depth_km"] = -0.8
        @test GraphSplit.validate_config(cfg) === cfg
        cfg["constraints"]["depth_bound"]["scope"] = "event_ids"
        @test_throws ErrorException GraphSplit.validate_config(cfg)
        cfg["constraints"]["depth_bound"]["event_ids"] = [10, 20]
        cfg["constraints"]["depth_bound"]["reflected_prelocation_restart"] = true
        @test GraphSplit.validate_config(cfg) === cfg
        cfg["constraints"]["depth_bound"]["pilot_shallow_margin_km"] = -1.0
        @test_throws ErrorException GraphSplit.validate_config(cfg)
        cfg["constraints"]["depth_bound"]["pilot_shallow_margin_km"] = 10.0
        cfg["run"]["prelocation"] = false
        @test_throws ErrorException GraphSplit.validate_config(cfg)

        cfg = GraphSplit.default_config()
        cfg["constraints"]["fixed_depth"]["enabled"] = true
        cfg["constraints"]["fixed_depth"]["scope"] = "all"
        cfg["constraints"]["fixed_depth"]["event_ids"] = [10]
        @test_throws ErrorException GraphSplit.validate_config(cfg)
        cfg["constraints"]["fixed_depth"]["event_ids"] = Int[]
        cfg["constraints"]["fixed_depth"]["depth_km"] = -0.8
        @test GraphSplit.validate_config(cfg) === cfg

        cfg = GraphSplit.default_config()
        cfg["catalog"]["origin_time_columns"] = [1, 2, 3]
        @test_throws ErrorException GraphSplit.validate_config(cfg)
        cfg["catalog"]["origin_time_columns"] = Int[]
        @test GraphSplit.validate_config(cfg) === cfg
        cfg["catalog"]["origin_time_decimals"] = 13
        @test_throws ErrorException GraphSplit.validate_config(cfg)

        cfg = GraphSplit.default_config()
        cfg["run"]["prelocation"] = false
        cfg["run"]["build_travel_times_only"] = true
        @test GraphSplit.validate_config(cfg) === cfg
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

        raw = zeros(3, 10)
        raw[:, 7] .= [35.9, 36.0, 36.1]
        raw[:, 8] .= [-117.8, -117.7, -117.6]
        raw[:, 9] .= [4.0, 5.0, 9.0]
        raw[:, 10] .= [10, 20, 30]
        catalog = GraphSplit.Catalog("memory", raw, Int64[10, 20, 30], raw[:, 7], raw[:, 8], raw[:, 9],
            7, 8, 9, 10, NaN, NaN, 6_371_000.0)
        stations = GraphSplit.Stations(["STA"], [36.0], [-117.7], [0.0], [0.0], [0.0], [0.0],
            NaN, NaN, 6_371_000.0)
        cfg = GraphSplit.default_config()
        cfg["coordinates"]["reference"] = "manual"
        cfg["coordinates"]["reference_latitude"] = 36.0
        cfg["coordinates"]["reference_longitude"] = -117.7

        catalog_state = GraphSplit.attach_coordinates!(stations, catalog, cfg)
        expected_centroid = (sum(catalog_state.x) / 3, sum(catalog_state.y) / 3,
            sum(catalog_state.z) / 3)
        cfg["initialization"]["mode"] = "common_centroid"
        centroid_state = GraphSplit.apply_initialization!(GraphSplit.copy_state(catalog_state), cfg)
        @test all(==(expected_centroid[1]), centroid_state.x)
        @test all(==(expected_centroid[2]), centroid_state.y)
        @test all(==(expected_centroid[3]), centroid_state.z)
        @test all(iszero, centroid_state.t0)

        cfg["initialization"]["mode"] = "common_manual"
        cfg["initialization"]["latitude"] = 36.05
        cfg["initialization"]["longitude"] = -117.65
        cfg["initialization"]["depth_km"] = 6.25
        manual_state = GraphSplit.apply_initialization!(GraphSplit.copy_state(catalog_state), cfg)
        manual_x, manual_y = GraphSplit.local_xy(36.05, -117.65, 36.0, -117.7, 6_371_000.0)
        @test manual_state.x == fill(manual_x, 3)
        @test manual_state.y == fill(manual_y, 3)
        @test manual_state.z == fill(6250.0, 3)

        original_range = GraphSplit.required_range(stations, catalog, catalog_state, :cartesian, 6_371_000.0)
        manual_state.x .= 200_000.0
        manual_range = GraphSplit.required_range(stations, catalog, manual_state, :cartesian, 6_371_000.0)
        @test manual_range > original_range
        @test manual_range >= 200_000.0
    end

    @testset "catalog origin times and cumulative shifts" begin
        ids = Int64[10, 20, 30]
        raw = [
            2024.0 12.0 31.0 23.0 59.0 59.9 35.0 -117.0 1.0 10.0
            2025.0 1.0 1.0 0.0 0.0 0.1 35.0 -117.0 2.0 20.0
            0.0 0.0 0.0 0.0 0.0 0.0 35.0 -117.0 3.0 30.0
        ]
        catalog = GraphSplit.Catalog("memory", raw, ids, raw[:, 7], raw[:, 8], raw[:, 9],
            7, 8, 9, 10, 35.0, -117.0, 6_371_000.0)
        catalog_state = GraphSplit.State(ids, zeros(3), zeros(3), 1000.0 .* raw[:, 9],
            zeros(3), 35.0, -117.0, 6_371_000.0)
        state = GraphSplit.copy_state(catalog_state)
        state.x .= [12.0, 24.0, 36.0]
        state.y .= [-2.0, -4.0, -6.0]
        state.z .+= [3.0, 6.0, 9.0]
        state.t0 .= [0.2, -0.2, 0.5]
        activity = BitVector([true, true, false])
        cfg = GraphSplit.default_config()

        mktempdir() do directory
            catalog_path = joinpath(directory, "catalog_dd.txt")
            GraphSplit.write_catalog(catalog_path, catalog, state, cfg;
                origin_time_activity=activity)
            written = GraphSplit.read_numeric_matrix(catalog_path)
            @test written[1, 1:6] == [2025.0, 1.0, 1.0, 0.0, 0.0, 0.1]
            @test written[2, 1:6] == [2024.0, 12.0, 31.0, 23.0, 59.0, 59.9]
            @test written[3, 1:6] == raw[3, 1:6]

            restarted = GraphSplit.copy_state(state)
            restarted.t0 .= [0.25, -0.15, 0.5]
            restart_catalog_path = joinpath(directory, "catalog_dd_restart.txt")
            GraphSplit.write_catalog(restart_catalog_path,
                GraphSplit.Catalog("restart", written, ids, written[:, 7], written[:, 8],
                    written[:, 9], 7, 8, 9, 10, 35.0, -117.0, 6_371_000.0),
                restarted, cfg; origin_time_basis=state.t0,
                origin_time_activity=BitVector([true, true, false]))
            restarted_written = GraphSplit.read_numeric_matrix(restart_catalog_path)
            @test restarted_written[1, 6] == 0.15
            @test restarted_written[2, 6] == 59.95
            @test restarted_written[3, 1:6] == raw[3, 1:6]

            cfg_without_calendar = deepcopy(cfg)
            cfg_without_calendar["catalog"]["origin_time_columns"] = Int[]
            no_calendar_path = joinpath(directory, "catalog_no_calendar.txt")
            GraphSplit.write_catalog(no_calendar_path, catalog, state, cfg_without_calendar;
                origin_time_activity=activity)
            @test GraphSplit.read_numeric_matrix(no_calendar_path)[:, 1:6] == raw[:, 1:6]

            invalid_activity = BitVector([false, false, true])
            @test_logs (:warn, r"Origin-time corrections could not be applied")
                GraphSplit.write_catalog(joinpath(directory, "catalog_invalid_calendar.txt"),
                    catalog, state, cfg; origin_time_activity=invalid_activity)

            first_shift = GraphSplit.cumulative_catalog_shift(
                GraphSplit.empty_catalog_shift(ids), catalog_state, state, activity)
            @test first_shift.dx_m == [12.0, 24.0, 0.0]
            @test first_shift.dy_m == [-2.0, -4.0, 0.0]
            @test first_shift.dz_m == [3.0, 6.0, 0.0]
            @test first_shift.dt0_s == [0.2, -0.2, 0.0]

            shift_path = joinpath(directory, "catalog_dd_dxdydzt0.txt")
            GraphSplit.write_catalog_shift(shift_path, first_shift)
            reordered = GraphSplit.read_catalog_shift(shift_path, Int64[30, 10, 20])
            @test reordered.event_id == Int64[30, 10, 20]
            @test reordered.dt0_s == [0.0, 0.2, -0.2]

            second_state = GraphSplit.copy_state(catalog_state)
            second_state.x .= state.x .+ [1.0, 2.0, 3.0]
            second_state.y .= state.y
            second_state.z .= state.z
            second_state.t0 .= [0.25, -0.15, 0.0]
            second_shift = GraphSplit.cumulative_catalog_shift(first_shift, state,
                second_state, BitVector([true, false, false]))
            @test second_shift.dx_m == [13.0, 24.0, 0.0]
            @test second_shift.dt0_s == [0.25, -0.2, 0.0]

            stations = GraphSplit.Stations(["STA"], [35.0], [-117.0], [0.0],
                [0.0], [0.0], [0.0], 35.0, -117.0, 6_371_000.0)
            obs = GraphSplit.Observations(Int32[1], Int32[2], zeros(1),
                ones(1), Int32[1], UInt8[1], Int32[1])
            travel_time = GraphSplit.ConstantTravelTime(:cartesian, 6000.0, 3500.0,
                6_371_000.0, 35.0, -117.0, 6_371_000.0)
            prediction = GraphSplit.predict_and_gradients(state, stations, obs, travel_time)[1]
            reconstructed = GraphSplit.copy_state(state)
            reconstructed.t0 .= GraphSplit.read_catalog_shift(shift_path, ids).dt0_s
            @test GraphSplit.predict_and_gradients(reconstructed, stations, obs, travel_time)[1] == prediction
        end
    end

    @testset "DD-only pipeline restart preserves origin-time state" begin
        mktempdir() do directory
            theta_directory = joinpath(directory, "theta")
            mkpath(theta_directory)
            catalog_path = joinpath(directory, "catalog.txt")
            stations_path = joinpath(directory, "stations.txt")
            theta_path = joinpath(theta_directory, "theta_STA_P.txt")
            open(catalog_path, "w") do io
                println(io, "2025 1 2 3 4 10.000000 35.000000 -117.000000 4.0000 1.0 1")
                println(io, "2025 1 2 3 4 10.000000 35.000000 -117.000000 4.0000 1.0 2")
                println(io, "2025 1 2 3 4 10.000000 35.000000 -117.000000 4.0000 1.0 3")
            end
            open(stations_path, "w") do io
                println(io, "STA 35.000000 -117.000000 0.0")
            end
            open(theta_path, "w") do io
                println(io, "1 0.010000 1")
                println(io, "2 -0.010000 1")
                println(io, "3 0.000000 1")
            end

            cfg = GraphSplit.default_config()
            cfg["io"]["catalog_file"] = catalog_path
            cfg["io"]["stations_file"] = stations_path
            cfg["io"]["theta_dir"] = theta_directory
            cfg["io"]["thetastd_dir"] = joinpath(directory, "thetastd")
            cfg["io"]["output_dir"] = joinpath(directory, "pass1")
            cfg["coordinates"]["reference"] = "manual"
            cfg["coordinates"]["reference_latitude"] = 35.0
            cfg["coordinates"]["reference_longitude"] = -117.0
            cfg["travel_time"]["type"] = "constant"
            cfg["travel_time"]["geometry"] = "cartesian"
            cfg["graph"]["neighbors"] = 2
            cfg["graph"]["maximum_degree"] = 2
            cfg["graph"]["mutual"] = false
            cfg["observations"]["minimum_theta_degree"] = 0
            cfg["observations"]["minimum_observations_per_pair"] = 1
            cfg["prelocation"]["max_outer_iterations"] = 6
            cfg["prelocation"]["min_outer_iterations"] = 1
            cfg["prelocation"]["verbose"] = false
            cfg["relocation"]["max_outer_iterations"] = 6
            cfg["relocation"]["min_outer_iterations"] = 1
            cfg["relocation"]["verbose"] = false
            @test GraphSplit.validate_config(cfg) === cfg

            pass1 = GraphSplit.run(cfg)
            pass1_catalog = joinpath(cfg["io"]["output_dir"], "catalog_dd.txt")
            pass1_shift = joinpath(cfg["io"]["output_dir"], "catalog_dd_dxdydzt0.txt")
            @test isfile(pass1_catalog)
            @test isfile(pass1_shift)
            @test maximum(abs.(pass1.dd_shift.dt0_s)) > 0.005

            cfg2 = deepcopy(cfg)
            cfg2["io"]["catalog_file"] = pass1_catalog
            cfg2["io"]["restart_shift_file"] = pass1_shift
            cfg2["io"]["output_dir"] = joinpath(directory, "pass2")
            cfg2["run"]["prelocation"] = false
            @test GraphSplit.validate_config(cfg2) === cfg2
            pass2 = GraphSplit.run(cfg2)

            @test pass2.pre_state.t0 ≈ pass1.dd_shift.dt0_s atol=1.0e-12
            output2 = GraphSplit.read_numeric_matrix(joinpath(cfg2["io"]["output_dir"],
                "catalog_dd.txt"))
            @test output2[:, 6] ≈ 10.0 .+ pass2.dd_shift.dt0_s atol=1.1e-6
            @test GraphSplit.read_numeric_matrix(joinpath(cfg2["io"]["output_dir"],
                "catalog_dd_filt.txt"))[:, 6] ≈ output2[:, 6] atol=1.0e-12
        end
    end

    @testset "depth constraints" begin
        ids = Int64[10, 20, 30]
        state = GraphSplit.State(ids, zeros(3), zeros(3), [1000.0, 2000.0, 3000.0],
            zeros(3), 36.0, -117.0, 6_371_000.0)
        stations = GraphSplit.Stations(["STA", "STB"], [36.0, 36.1], [-117.0, -117.1],
            [800.0, 600.0], zeros(2), zeros(2), [-800.0, -600.0],
            36.0, -117.0, 6_371_000.0)

        cfg = GraphSplit.default_config()
        fixed = cfg["constraints"]["fixed_depth"]
        fixed["enabled"] = true
        fixed["scope"] = "event_ids"
        fixed["event_ids"] = [20]
        fixed["depth_km"] = -0.8
        GraphSplit.apply_fixed_depth_constraints!(state, cfg)
        @test state.z == [1000.0, -800.0, 3000.0]
        active = GraphSplit.parameter_active_mask(state, cfg)
        @test !active[2length(state) + 2]
        @test active[2length(state) + 1]

        fixed["enabled"] = false
        bound = cfg["constraints"]["depth_bound"]
        bound["enabled"] = true
        bound["minimum_depth_km"] = -0.8
        bound["scope"] = "all"
        bound["reflected_prelocation_restart"] = true
        bound["mirror_plane"] = "manual"
        bound["mirror_depth_km"] = -0.8
        state.z .= [-1200.0, 0.0, 3000.0]
        reflected = GraphSplit.reflect_depth_violations!(state, stations, cfg)
        @test reflected == BitVector([true, false, false])
        @test state.z[1] == -400.0
        @test !any(GraphSplit.depth_bound_violation_mask(state, cfg))
        @test -10_800.0 in GraphSplit.constraint_required_depths(state, stations, cfg)

        cfg["coordinates"]["event_vertical"] = "negative_depth"
        bound["minimum_depth_km"] = 1.0
        state.z .= [-500.0, -1500.0, -2000.0]
        @test GraphSplit.depth_bound_violation_mask(state, cfg) == BitVector([true, false, false])

        cfg = GraphSplit.default_config()
        bound = cfg["constraints"]["depth_bound"]
        bound["enabled"] = true
        bound["minimum_depth_km"] = 0.0
        bound["scope"] = "event_ids"
        bound["event_ids"] = [2]
        linear_state = GraphSplit.State(Int64[1, 2], zeros(2), zeros(2), [1.0, 1.0],
            zeros(2), 0.0, 0.0, 6_371_000.0)
        system = GraphSplit.LinearizedSystem(Int32[1], Int32[2], [0.0], [0.0], [0.0],
            [0.0], [0.0], [-1.0], [1.0], 2, :none, 1.0, 5000.0)
        rhs_outward = GraphSplit.apply_At(system, [-2.0])
        base_active = falses(8)
        base_active[6] = true
        solver_cfg = deepcopy(cfg["relocation"])
        solver_cfg["step_damping"] = 1.0
        update, _, converged, engaged, contact = GraphSplit.solve_model_update(system,
            rhs_outward, 1.0e-8, linear_state, solver_cfg, cfg, base_active,
            "diagonal", "direct", 1.0e-12, 20)
        @test converged
        @test update[6] ≈ -1.0 atol=1.0e-10
        @test engaged == BitVector([false, true])
        @test contact == BitVector([false, true])

        linear_state.z[2] = 0.0
        rhs_inward = GraphSplit.apply_At(system, [2.0])
        update, _, converged, engaged, contact = GraphSplit.solve_model_update(system,
            rhs_inward, 1.0e-8, linear_state, solver_cfg, cfg, base_active,
            "diagonal", "direct", 1.0e-12, 20)
        @test converged
        @test update[6] > 0.0
        @test !any(engaged)
        @test !any(contact)
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

    @testset "uncertainty resampling and sidecars" begin
        empty_entries = Dict{Int64,GraphSplit.ThetaEntry}()
        groups = [
            GraphSplit.ThetaGroup("STA", UInt8(1), "theta_STA_P.txt", Int64[], copy(empty_entries), false, false),
            GraphSplit.ThetaGroup("STA", UInt8(2), "theta_STA_S.txt", Int64[], copy(empty_entries), false, false),
            GraphSplit.ThetaGroup("STB", UInt8(1), "theta_STB_P.txt", Int64[], copy(empty_entries), false, false),
        ]
        obs = GraphSplit.Observations(Int32[1, 1, 1, 2], Int32[2, 2, 2, 3],
            [0.1, 0.2, 0.3, 0.4], fill(0.01, 4), Int32[1, 1, 1, 2],
            UInt8[1, 1, 2, 1], Int32[1, 1, 2, 3])
        station_phase_blocks, station_phase_labels = GraphSplit.bootstrap_blocks(obs, groups, "station_phase")
        @test station_phase_blocks == Int32[1, 1, 2, 3]
        @test station_phase_labels == ["theta_STA_P.txt", "theta_STA_S.txt", "theta_STB_P.txt"]
        station_blocks, station_labels = GraphSplit.bootstrap_blocks(obs, groups, "station")
        @test station_blocks == Int32[1, 1, 1, 2]
        @test station_labels == ["STA", "STB"]

        sampled = GraphSplit.resample_observations(obs, station_phase_blocks, [2, 0, 1])
        @test length(sampled) == 5
        @test count(==(Int32(1)), sampled.group) == 4
        @test count(==(Int32(2)), sampled.group) == 0
        @test count(==(Int32(3)), sampled.group) == 1
        counts = GraphSplit.draw_bootstrap_counts(MersenneTwister(12), 7)
        @test sum(counts) == 7
        @test length(counts) == 7

        @test GraphSplit.sample_quantile([0.0, 10.0], 0.25) == 2.5
        covariance = GraphSplit.sample_covariance3([-1.0, 1.0], [-2.0, 2.0], [-3.0, 3.0])
        @test covariance ≈ [2.0 4.0 6.0; 4.0 8.0 12.0; 6.0 12.0 18.0]

        mktempdir() do directory
            ids = Int64[10, 20]
            mask = BitVector([false, true])
            linear_covariance = fill(NaN, 3, 3, 2)
            linear_covariance[:, :, 2] .= Diagonal([1.0, 4.0, 9.0])
            estimate = (covariance=linear_covariance,)
            linear_path = joinpath(directory, "linerrxyz.txt")
            GraphSplit.write_linearized_uncertainty(linear_path, ids, mask, estimate)
            linear_lines = readlines(linear_path)
            @test length(linear_lines) == 2
            @test startswith(linear_lines[2], "20 1 2 3 ")

            nominal = GraphSplit.State(ids, [0.0, 100.0], [0.0, 200.0], [0.0, 300.0],
                zeros(2), 36.0, -117.0, 6_371_000.0)
            samples = (x=[99.0 100.0 102.0], y=[198.0 200.0 204.0],
                z=[297.0 300.0 306.0], t0=zeros(1, 3))
            cfg = GraphSplit.default_config()
            summary_path = joinpath(directory, "booterrxyz.txt")
            GraphSplit.write_bootstrap_summary(summary_path, ids, [2], nominal, samples, 3, cfg)
            summary_lines = readlines(summary_path)
            @test length(summary_lines) == 2
            @test split(summary_lines[2])[1:3] == ["20", "3", "1.00000000"]
        end
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
