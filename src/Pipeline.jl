function prepare_output_directory(cfg::AbstractDict)
    output = String(cfgget(cfg, "io", "output_dir"))
    mkpath(output)
    primary = ("catalog_dd.txt", "catalog_dd_filt.txt", "catalog_preloc.txt", "catalog_preloc_filt.txt",
        "catalog_dd_dxdydzt0.txt", "catalog_preloc_dxdydzt0.txt")
    if !Bool(cfgget(cfg, "run", "overwrite"; default=true)) && any(isfile(joinpath(output, file)) for file in primary)
        error("Output directory already contains GraphSplit results and run.overwrite=false: $output")
    end
    return output
end

"Run the branch-search pilot without hiding departures from its lookup grid."
function solve_reflected_prelocation_pilot(initial::State, stations::Stations,
        observations::Observations, travel_time::AbstractTravelTimeModel,
        solver_cfg::AbstractDict, pilot_cfg::AbstractDict)
    travel_time isa TravelTimeTable ||
        return solve_relocation(initial, stations, observations, travel_time, solver_cfg, pilot_cfg)
    previous_clamp = travel_time.clamp
    travel_time.clamp = false
    try
        return solve_relocation(initial, stations, observations, travel_time, solver_cfg, pilot_cfg)
    catch exception
        if occursin("Travel-time lookup returned a non-finite value", sprint(showerror, exception))
            error("The unconstrained reflected-depth pilot left the travel-time table. Increase constraints.depth_bound.pilot_shallow_margin_km and rebuild the table; clamping is deliberately disabled for this branch search")
        end
        rethrow()
    finally
        travel_time.clamp = previous_clamp
    end
end

function apply_pin_reference!(state::State, stations::Stations, cfg::AbstractDict)
    reference_path = String(cfgget(cfg, "gauge", "pin_reference_catalog"; default=""))
    isempty(reference_path) && return state
    reference = read_catalog(reference_path, cfg)
    reference_x, reference_y = local_xy(reference.lat, reference.lon,
        state.ref_lat, state.ref_lon, state.ref_radius_m)
    reference_z = catalog_depth_to_internal(reference.depth_km, cfg)
    ref_rows = Dict(id => row for (row, id) in enumerate(reference.event_id))
    state_rows = Dict(id => row for (row, id) in enumerate(state.event_id))
    pin_ids = Int64.(cfgget(cfg, "gauge", "pin_event_ids"; default=Int[]))
    fields = lowercase(String(cfgget(cfg, "gauge", "pin_fields"; default="xyz")))
    for id in pin_ids
        haskey(ref_rows, id) || error("Pinned ID $id is absent from gauge.pin_reference_catalog")
        haskey(state_rows, id) || error("Pinned ID $id is absent from the working catalog")
        source, target = ref_rows[id], state_rows[id]
        occursin('x', fields) && (state.x[target] = reference_x[source])
        occursin('y', fields) && (state.y[target] = reference_y[source])
        occursin('z', fields) && (state.z[target] = reference_z[source])
    end
    return state
end

function write_summary(path::AbstractString, cfg::AbstractDict, catalog::Catalog, stations::Stations,
        groups, graph::EventGraph, pre::Observations, dd::Observations, pre_stats::SolveStats, dd_stats::SolveStats,
        travel_time::AbstractTravelTimeModel; reflected_depths::BitVector=falses(length(catalog)),
        prelocation_pilot::Union{Nothing,SolveStats}=nothing)
    geometry = string(travel_time.geometry)
    summary_state = State(copy(catalog.event_id), zeros(length(catalog)), zeros(length(catalog)),
        zeros(length(catalog)), zeros(length(catalog)), catalog.ref_lat, catalog.ref_lon,
        catalog.ref_radius_m)
    summary = Dict{String,Any}(
        "run" => Dict("completed_utc" => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ"),
            "config_file" => String(get(cfg, "_config_file", ""))),
        "inputs" => Dict("events" => length(catalog), "stations" => length(stations), "theta_groups" => length(groups)),
        "travel_time" => Dict{String,Any}("geometry" => geometry, "type" => travel_time isa TravelTime3D ? "3d" : travel_time isa TravelTimeTable ? "lookup" : "constant"),
        "initialization" => Dict("mode" => String(cfgget(cfg, "initialization", "mode"; default="catalog"))),
        "catalog_output" => Dict(
            "origin_time_columns" => Int.(cfgget(cfg, "catalog", "origin_time_columns"; default=Int[])),
            "origin_time_decimals" => Int(cfgget(cfg, "catalog", "origin_time_decimals"; default=6)),
            "shift_spatial_units" => "m",
            "shift_time_units" => "s",
            "shift_semantics" => "cumulative_from_first_pass_input",
            "shift_horizontal_axes" => "local_east_north",
            "shift_vertical_convention" => String(cfgget(cfg, "coordinates", "event_vertical"; default="positive_depth")),
            "shift_reference_latitude" => catalog.ref_lat,
            "shift_reference_longitude" => catalog.ref_lon,
            "restart_shift_file" => String(cfgget(cfg, "io", "restart_shift_file"; default="")),
        ),
        "prelocation" => Dict("observations" => length(pre), "iterations" => pre_stats.iterations,
            "final_robust_rms_s" => isempty(pre_stats.rms_s) ? NaN : pre_stats.rms_s[end]),
        "relocation" => Dict("graph_edges" => length(graph), "graph_components" => graph.ncomp,
            "observations" => length(dd), "iterations" => dd_stats.iterations,
            "final_robust_rms_s" => isempty(dd_stats.rms_s) ? NaN : dd_stats.rms_s[end]),
        "uncertainty" => Dict("method" => String(cfgget(cfg, "uncertainty", "method"; default="none"))),
        "depth_constraints" => Dict(
            "bound_enabled" => Bool(get(depth_bound_options(cfg), "enabled", false)),
            "bound_events" => count(depth_bound_mask(summary_state, cfg)),
            "fixed_depth_enabled" => Bool(get(fixed_depth_options(cfg), "enabled", false)),
            "fixed_depth_events" => count(fixed_depth_mask(summary_state, cfg)),
            "reflected_prelocation_events" => count(reflected_depths),
            "prelocation_bound_active_events" => count(pre_stats.depth_bound_active),
            "relocation_bound_active_events" => count(dd_stats.depth_bound_active),
            "unconstrained_prelocation_pilot_ran" => prelocation_pilot !== nothing,
        ),
    )
    travel_time isa TravelTime3D && merge!(summary["travel_time"], travel_time.metadata)
    open(path, "w") do io
        TOML.print(io, summary; sorted=true)
    end
end

function build_travel_times(cfg::Dict{String,Any}; force::Bool=true)
    catalog = read_catalog(String(cfgget(cfg, "io", "catalog_file")), cfg)
    stations = read_stations(String(cfgget(cfg, "io", "stations_file")))
    state = attach_coordinates!(stations, catalog, cfg)
    apply_initialization!(state, cfg)
    apply_pin_reference!(state, stations, cfg)
    apply_fixed_depth_constraints!(state, cfg)
    validate_depth_constraint_state!(state, cfg; allow_free_outside=true)
    table = prepare_travel_time(cfg, stations, catalog, state; force_build=force)
    table isa TravelTimeTable && @printf("Travel-time table ready: %s (%s)\n", table.file, string(table.geometry))
    return table
end

"Run one complete TOML-driven GraphSplit relocation."
function run(cfg::Dict{String,Any})
    output = prepare_output_directory(cfg)
    catalog = read_catalog(String(cfgget(cfg, "io", "catalog_file")), cfg)
    stations = read_stations(String(cfgget(cfg, "io", "stations_file")))
    catalog_state = attach_coordinates!(stations, catalog, cfg)
    initial = copy_state(catalog_state)
    do_prelocation = Bool(cfgget(cfg, "run", "prelocation"; default=true))
    build_only = Bool(cfgget(cfg, "run", "build_travel_times_only"; default=false))
    previous_shift = empty_catalog_shift(catalog.event_id)
    if !do_prelocation && !build_only
        restart_path = String(cfgget(cfg, "io", "restart_shift_file"; default=""))
        isempty(restart_path) && error("run.prelocation=false requires io.restart_shift_file from the catalog used as the Stage-2 seed")
        previous_shift = read_catalog_shift(restart_path, catalog.event_id)
        initial.t0 .= previous_shift.dt0_s
        @printf("Loaded cumulative shifts for %d events from %s\n", length(previous_shift), restart_path)
    end
    apply_initialization!(initial, cfg)
    apply_pin_reference!(initial, stations, cfg)
    apply_fixed_depth_constraints!(initial, cfg)
    reflected_restart = Bool(get(depth_bound_options(cfg), "reflected_prelocation_restart", false))
    validate_depth_constraint_state!(initial, cfg; allow_free_outside=reflected_restart)
    travel_time = prepare_travel_time(cfg, stations, catalog, initial)
    if build_only
        result = (geometry=travel_time.geometry,
            table_file=travel_time isa TravelTime3D ? travel_time.metadata["cache_dir"] : travel_time isa TravelTimeTable ? travel_time.file : "",
            type=travel_time isa TravelTime3D ? :grid3d : travel_time isa TravelTimeTable ? :lookup : :constant)
        travel_time isa TravelTimeTable && close(travel_time.io)
        return result
    end
    groups = load_theta_folder(String(cfgget(cfg, "io", "theta_dir")), String(cfgget(cfg, "io", "thetastd_dir")))
    bias_model = String(cfgget(cfg, "experimental", "bias", "apply_model_file"; default=""))
    isempty(bias_model) || apply_bias_model!(groups, bias_model, initial, catalog)

    reflected_depths = falses(length(catalog))
    prelocation_pilot_stats = nothing
    if do_prelocation
        println("\n=== Stage 1: theta prelocation ===")
        pre_observations = build_star_observations(groups, stations, catalog, cfg)
        if isempty(pre_observations)
            lowercase(String(cfgget(cfg, "initialization", "mode"; default="catalog"))) == "catalog" ||
                error("No Stage-1 observations survived, so the common initial hypocenter cannot be separated into an event graph")
            @warn "No Stage-1 observations survived; using the input catalog as the Stage-2 seed"
            pre_state = copy_state(initial)
            if reflected_restart
                reflected_depths .= reflect_depth_violations!(pre_state, stations, cfg)
                any(reflected_depths) && @printf("Reflected %d initial depths into the admissible branch before graph construction\n",
                    count(reflected_depths))
            end
            pre_stats = SolveStats(0, Float64[], Float64[], Float64[], Int[], true,
                zeros(Int, length(catalog)), falses(length(catalog)))
        else
            if reflected_restart
                println("--- Unconstrained Stage-1 pilot for mirrored-depth branch search ---")
                pilot_cfg = deepcopy(cfg)
                pilot_cfg["constraints"]["depth_bound"]["enabled"] = false
                pilot_cfg["constraints"]["depth_bound"]["reflected_prelocation_restart"] = false
                pilot_state, pilot_stats = solve_reflected_prelocation_pilot(initial, stations,
                    pre_observations, travel_time, cfg["prelocation"], pilot_cfg)
                restart_state = copy_state(pilot_state)
                reflected_depths .= reflect_depth_violations!(restart_state, stations, cfg)
                if any(reflected_depths)
                    prelocation_pilot_stats = pilot_stats
                    @printf("Reflected %d Stage-1 pilot depths; rerunning Stage 1 with the active bound\n",
                        count(reflected_depths))
                    pre_state, pre_stats = solve_relocation(restart_state, stations, pre_observations,
                        travel_time, cfg["prelocation"], cfg)
                else
                    println("Stage-1 pilot ended inside the depth bound; no reflected restart was needed")
                    pre_state, pre_stats = pilot_state, pilot_stats
                end
            else
                pre_state, pre_stats = solve_relocation(initial, stations, pre_observations, travel_time,
                    cfg["prelocation"], cfg)
            end
        end
    else
        println("\n=== Stage 1 skipped: input catalog is the Stage-2 seed ===")
        pre_observations = Observations()
        pre_state = copy_state(initial)
        pre_stats = SolveStats(0, Float64[], Float64[], Float64[], Int[], true,
            zeros(Int, length(catalog)), falses(length(catalog)))
    end
    pre_activity = do_prelocation && !isempty(pre_observations) ?
        event_activity(pre_observations, length(catalog)) : falses(length(catalog))
    pre_shift = cumulative_catalog_shift(previous_shift, catalog_state, pre_state, pre_activity)
    write_catalog(joinpath(output, "catalog_preloc.txt"), catalog, pre_state, cfg;
        origin_time_basis=previous_shift.dt0_s, origin_time_activity=pre_activity)
    write_catalog_shift(joinpath(output, "catalog_preloc_dxdydzt0.txt"), pre_shift)
    pre_mask = do_prelocation && !isempty(pre_observations) ? event_activity(pre_observations, length(catalog)) : trues(length(catalog))
    Bool(cfgget(cfg, "output", "write_filtered_catalogs"; default=true)) &&
        write_catalog(joinpath(output, "catalog_preloc_filt.txt"), catalog, pre_state, cfg;
            mask=pre_mask, origin_time_basis=previous_shift.dt0_s,
            origin_time_activity=pre_activity)

    println("\n=== Stage 2: sparse graph and DD relocation ===")
    graph = build_event_graph(pre_state, cfg)
    dd_observations = build_dd_observations(groups, stations, catalog, graph, cfg)
    isempty(dd_observations) && error("No Stage-2 observations survived. Check station names, thetaStd degree filtering, pair support, and graph radius")
    dd_state, dd_stats = solve_relocation(pre_state, stations, dd_observations, travel_time,
        cfg["relocation"], cfg)
    dd_mask = event_activity(dd_observations, length(catalog))
    final_activity = pre_activity .| dd_mask
    dd_shift = cumulative_catalog_shift(previous_shift, catalog_state, dd_state, final_activity)
    write_catalog(joinpath(output, "catalog_dd.txt"), catalog, dd_state, cfg;
        origin_time_basis=previous_shift.dt0_s, origin_time_activity=final_activity)
    write_catalog_shift(joinpath(output, "catalog_dd_dxdydzt0.txt"), dd_shift)
    Bool(cfgget(cfg, "output", "write_filtered_catalogs"; default=true)) &&
        write_catalog(joinpath(output, "catalog_dd_filt.txt"), catalog, dd_state, cfg;
            mask=dd_mask, origin_time_basis=previous_shift.dt0_s,
            origin_time_activity=final_activity)

    if Bool(cfgget(cfg, "output", "write_graph_metadata"; default=true))
        write_graph_metadata(joinpath(output, "catalog_dd_graphmeta.csv"), catalog, graph, pre_observations, dd_observations)
    end
    if Bool(cfgget(cfg, "output", "write_solver_history"; default=true))
        write_solver_history(joinpath(output, "solver_history.csv"), pre_stats, dd_stats;
            prelocation_pilot=prelocation_pilot_stats)
    end
    if Bool(get(depth_bound_options(cfg), "enabled", false)) || Bool(get(fixed_depth_options(cfg), "enabled", false))
        write_depth_constraint_status(joinpath(output, "depth_constraint_status.csv"), catalog,
            dd_state, cfg, pre_stats, dd_stats, reflected_depths)
    end
    uncertainty_method = lowercase(String(cfgget(cfg, "uncertainty", "method"; default="none")))
    linearized_uncertainty = nothing
    bootstrap_uncertainty = nothing
    if uncertainty_method in ("linearized", "both")
        println("\n=== Uncertainty: regularized linearized covariance ===")
        linearized_uncertainty = estimate_linearized_uncertainty(dd_state, stations, dd_observations,
            travel_time, dd_mask, cfg; depth_bound_active=dd_stats.depth_bound_active)
        write_linearized_uncertainty(joinpath(output, "linerrxyz.txt"), catalog.event_id,
            dd_mask, linearized_uncertainty)
    end
    if uncertainty_method in ("bootstrap", "both")
        println("\n=== Uncertainty: station-phase block bootstrap ===")
        bootstrap_uncertainty = run_bootstrap_uncertainty(catalog, groups, dd_state, stations,
            dd_observations, travel_time, dd_mask, output, cfg;
            origin_time_basis=previous_shift.dt0_s)
    end
    bias_report = NamedTuple[]
    if Bool(cfgget(cfg, "experimental", "bias", "enabled"; default=false))
        println("\n=== Experimental diagnostic: theta bias scan ===")
        bias_report = scan_theta_bias(groups, dd_state, stations, catalog, travel_time, cfg)
        write_bias_report(joinpath(output, "theta_bias_report.csv"), bias_report)
    end
    if Bool(cfgget(cfg, "output", "write_run_summary"; default=true))
        write_summary(joinpath(output, "run_summary.toml"), cfg, catalog, stations, groups, graph,
            pre_observations, dd_observations, pre_stats, dd_stats, travel_time;
            reflected_depths=reflected_depths, prelocation_pilot=prelocation_pilot_stats)
    end
    travel_time isa TravelTimeTable && close(travel_time.io)
    @printf("\nGraphSplit complete. Outputs: %s\n", output)
    return (catalog=catalog, stations=stations, groups=groups, pre_state=pre_state, dd_state=dd_state,
        graph=graph, pre_observations=pre_observations, dd_observations=dd_observations,
        pre_stats=pre_stats, dd_stats=dd_stats, linearized_uncertainty=linearized_uncertainty,
        bootstrap_uncertainty=bootstrap_uncertainty, bias_report=bias_report,
        pre_shift=pre_shift, dd_shift=dd_shift,
        reflected_depths=reflected_depths, prelocation_pilot_stats=prelocation_pilot_stats,
        output_dir=output)
end
