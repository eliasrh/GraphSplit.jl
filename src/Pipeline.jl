function prepare_output_directory(cfg::AbstractDict)
    output = String(cfgget(cfg, "io", "output_dir"))
    mkpath(output)
    primary = ("catalog_dd.txt", "catalog_dd_filt.txt", "catalog_preloc.txt", "catalog_preloc_filt.txt")
    if !Bool(cfgget(cfg, "run", "overwrite"; default=true)) && any(isfile(joinpath(output, file)) for file in primary)
        error("Output directory already contains GraphSplit results and run.overwrite=false: $output")
    end
    return output
end

function apply_pin_reference!(state::State, stations::Stations, cfg::AbstractDict)
    reference_path = String(cfgget(cfg, "gauge", "pin_reference_catalog"; default=""))
    isempty(reference_path) && return state
    reference = read_catalog(reference_path, cfg)
    reference_x, reference_y = local_xy(reference.lat, reference.lon,
        state.ref_lat, state.ref_lon, state.ref_radius_m)
    event_vertical = lowercase(String(cfgget(cfg, "coordinates", "event_vertical"; default="positive_depth")))
    z0 = Float64(cfgget(cfg, "coordinates", "event_z0_m"; default=0.0))
    reference_z = event_vertical == "positive_depth" ? reference.depth_km .* 1000.0 :
        event_vertical == "negative_depth" ? -reference.depth_km .* 1000.0 :
        event_vertical == "positive_depth_plus_z0" ? reference.depth_km .* 1000.0 .+ z0 :
        error("Unknown coordinates.event_vertical: $event_vertical")
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
        travel_time::AbstractTravelTimeModel)
    geometry = string(travel_time.geometry)
    summary = Dict{String,Any}(
        "run" => Dict("completed_utc" => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ"),
            "config_file" => String(get(cfg, "_config_file", ""))),
        "inputs" => Dict("events" => length(catalog), "stations" => length(stations), "theta_groups" => length(groups)),
        "travel_time" => Dict("geometry" => geometry, "type" => travel_time isa TravelTimeTable ? "lookup" : "constant"),
        "prelocation" => Dict("observations" => length(pre), "iterations" => pre_stats.iterations,
            "final_robust_rms_s" => isempty(pre_stats.rms_s) ? NaN : pre_stats.rms_s[end]),
        "relocation" => Dict("graph_edges" => length(graph), "graph_components" => graph.ncomp,
            "observations" => length(dd), "iterations" => dd_stats.iterations,
            "final_robust_rms_s" => isempty(dd_stats.rms_s) ? NaN : dd_stats.rms_s[end]),
    )
    open(path, "w") do io
        TOML.print(io, summary; sorted=true)
    end
end

function build_travel_times(cfg::Dict{String,Any}; force::Bool=true)
    catalog = read_catalog(String(cfgget(cfg, "io", "catalog_file")), cfg)
    stations = read_stations(String(cfgget(cfg, "io", "stations_file")))
    state = attach_coordinates!(stations, catalog, cfg)
    table = prepare_travel_time(cfg, stations, catalog, state; force_build=force)
    table isa TravelTimeTable && @printf("Travel-time table ready: %s (%s)\n", table.file, string(table.geometry))
    return table
end

"Run one complete TOML-driven GraphSplit relocation."
function run(cfg::Dict{String,Any})
    output = prepare_output_directory(cfg)
    catalog = read_catalog(String(cfgget(cfg, "io", "catalog_file")), cfg)
    stations = read_stations(String(cfgget(cfg, "io", "stations_file")))
    initial = attach_coordinates!(stations, catalog, cfg)
    apply_pin_reference!(initial, stations, cfg)
    travel_time = prepare_travel_time(cfg, stations, catalog, initial)
    if Bool(cfgget(cfg, "run", "build_travel_times_only"; default=false))
        result = (geometry=travel_time.geometry,
            table_file=travel_time isa TravelTimeTable ? travel_time.file : "",
            type=travel_time isa TravelTimeTable ? :lookup : :constant)
        travel_time isa TravelTimeTable && close(travel_time.io)
        return result
    end
    groups = load_theta_folder(String(cfgget(cfg, "io", "theta_dir")), String(cfgget(cfg, "io", "thetastd_dir")))
    bias_model = String(cfgget(cfg, "experimental", "bias", "apply_model_file"; default=""))
    isempty(bias_model) || apply_bias_model!(groups, bias_model, initial, catalog)

    do_prelocation = Bool(cfgget(cfg, "run", "prelocation"; default=true))
    if do_prelocation
        println("\n=== Stage 1: theta prelocation ===")
        pre_observations = build_star_observations(groups, stations, catalog, cfg)
        if isempty(pre_observations)
            @warn "No Stage-1 observations survived; using the input catalog as the Stage-2 seed"
            pre_state = copy_state(initial)
            pre_stats = SolveStats(0, Float64[], Float64[], Float64[], Int[], true)
        else
            pre_state, pre_stats = solve_relocation(initial, stations, pre_observations, travel_time,
                cfg["prelocation"], cfg)
        end
    else
        println("\n=== Stage 1 skipped: input catalog is the Stage-2 seed ===")
        pre_observations = Observations()
        pre_state = copy_state(initial)
        pre_stats = SolveStats(0, Float64[], Float64[], Float64[], Int[], true)
    end
    write_catalog(joinpath(output, "catalog_preloc.txt"), catalog, pre_state, cfg)
    pre_mask = do_prelocation && !isempty(pre_observations) ? event_activity(pre_observations, length(catalog)) : trues(length(catalog))
    Bool(cfgget(cfg, "output", "write_filtered_catalogs"; default=true)) &&
        write_catalog(joinpath(output, "catalog_preloc_filt.txt"), catalog, pre_state, cfg; mask=pre_mask)

    println("\n=== Stage 2: sparse graph and DD relocation ===")
    graph = build_event_graph(pre_state, cfg)
    dd_observations = build_dd_observations(groups, stations, catalog, graph, cfg)
    isempty(dd_observations) && error("No Stage-2 observations survived. Check station names, thetaStd degree filtering, pair support, and graph radius")
    dd_state, dd_stats = solve_relocation(pre_state, stations, dd_observations, travel_time,
        cfg["relocation"], cfg)
    write_catalog(joinpath(output, "catalog_dd.txt"), catalog, dd_state, cfg)
    dd_mask = event_activity(dd_observations, length(catalog))
    Bool(cfgget(cfg, "output", "write_filtered_catalogs"; default=true)) &&
        write_catalog(joinpath(output, "catalog_dd_filt.txt"), catalog, dd_state, cfg; mask=dd_mask)

    if Bool(cfgget(cfg, "output", "write_graph_metadata"; default=true))
        write_graph_metadata(joinpath(output, "catalog_dd_graphmeta.csv"), catalog, graph, pre_observations, dd_observations)
    end
    if Bool(cfgget(cfg, "output", "write_solver_history"; default=true))
        write_solver_history(joinpath(output, "solver_history.csv"), pre_stats, dd_stats)
    end
    bias_report = NamedTuple[]
    if Bool(cfgget(cfg, "experimental", "bias", "enabled"; default=false))
        println("\n=== Experimental diagnostic: theta bias scan ===")
        bias_report = scan_theta_bias(groups, dd_state, stations, catalog, travel_time, cfg)
        write_bias_report(joinpath(output, "theta_bias_report.csv"), bias_report)
    end
    if Bool(cfgget(cfg, "output", "write_run_summary"; default=true))
        write_summary(joinpath(output, "run_summary.toml"), cfg, catalog, stations, groups, graph,
            pre_observations, dd_observations, pre_stats, dd_stats, travel_time)
    end
    travel_time isa TravelTimeTable && close(travel_time.io)
    @printf("\nGraphSplit complete. Outputs: %s\n", output)
    return (catalog=catalog, stations=stations, groups=groups, pre_state=pre_state, dd_state=dd_state,
        graph=graph, pre_observations=pre_observations, dd_observations=dd_observations,
        pre_stats=pre_stats, dd_stats=dd_stats, bias_report=bias_report, output_dir=output)
end
