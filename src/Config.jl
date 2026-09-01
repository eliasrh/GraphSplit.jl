"Return the complete set of GraphSplit defaults. TOML files override only what they contain."
function default_config()
    solver_common = Dict{String,Any}(
        "max_outer_iterations" => 20,
        "min_outer_iterations" => 3,
        "linear_solver" => "pcg",
        "preconditioner" => "block_jacobi",
        "inner_tolerance" => 1.0e-5,
        "inner_max_iterations" => 250,
        "huber_k" => 1.345,
        "min_sigma_s" => 0.002,
        "step_damping" => 0.75,
        "damping_lambda" => 1.0e-3,
        "stop_step_rms_m" => 0.10,
        "stop_step_rms_s" => 1.0e-4,
        "stop_rms_improvement_s" => 1.0e-5,
        "stop_stall_iterations" => 5,
        "max_event_step_m" => 0.0,
        "max_origin_step_s" => 0.0,
        "verbose" => true,
    )
    dd = deepcopy(solver_common)
    dd["max_outer_iterations"] = 60
    dd["inner_tolerance"] = 1.0e-4
    dd["inner_max_iterations"] = 300
    dd["damping_lambda"] = 5.0e-4

    return Dict{String,Any}(
        "io" => Dict{String,Any}(
            "catalog_file" => "catalog.txt",
            "stations_file" => "stations.txt",
            "theta_dir" => "theta",
            "thetastd_dir" => "thetastd",
            "output_dir" => "graphsplit_output",
        ),
        "run" => Dict{String,Any}(
            "prelocation" => true,
            "overwrite" => true,
            "build_travel_times_only" => false,
        ),
        "catalog" => Dict{String,Any}(
            "latitude_column" => 7,
            "longitude_column" => 8,
            "depth_column" => 9,
            "event_id_column" => -1,
        ),
        "initialization" => Dict{String,Any}(
            "mode" => "catalog",
            "latitude" => 0.0,
            "longitude" => 0.0,
            "depth_km" => 0.0,
        ),
        "coordinates" => Dict{String,Any}(
            "reference" => "stations_mean",
            "reference_latitude" => 0.0,
            "reference_longitude" => 0.0,
            "earth_radius_m" => 6_371_000.0,
            "station_vertical" => "depth_from_elevation",
            "event_vertical" => "positive_depth",
            "event_z0_m" => 0.0,
        ),
        "travel_time" => Dict{String,Any}(
            "type" => "lookup",
            "table_file" => "lookuptable/graphsplit_tt.gstt",
            "velocity_model_file" => "vm.txt",
            "geometry" => "auto",
            "build_geometry" => "cartesian",
            "auto_build" => true,
            "auto_rebuild" => true,
            "clamp_to_grid" => true,
            "vp_ms" => 6000.0,
            "vs_ms" => 3464.0,
            "earth_radius_m" => 0.0,
        ),
        "lookup" => Dict{String,Any}(
            "horizontal_step_m" => 125.0,
            "depth_step_m" => 50.0,
            "station_depth_step_m" => 100.0,
            "maximum_distance_km" => 0.0,
            "minimum_depth_km" => 0.0,
            "maximum_depth_km" => 0.0,
            "distance_margin_km" => 1.0,
            "depth_margin_km" => 50.0,
            "maximum_sweeps" => 32,
            "sweep_tolerance_s" => 1.0e-7,
            "threaded" => false,
            "maximum_table_gib" => 8.0,
            "allow_large_table" => false,
        ),
        "prelocation" => solver_common,
        "relocation" => dd,
        "observations" => Dict{String,Any}(
            "minimum_theta_degree" => 6,
            "maximum_sigma_s" => 0.0,
            "minimum_observations_per_pair" => 6,
            "minimum_observations_per_event" => 0,
            "minimum_star_observations_per_event" => 0,
            "minimum_component_size" => 0,
        ),
        "gauge" => Dict{String,Any}(
            "mode" => "zero_mean",
            "zero_mean" => "origin_time",
            "constraint_weight" => 10.0,
            "reference_velocity_ms" => 5000.0,
            "pin_event_ids" => Int[],
            "pin_fields" => "xyz",
            "pin_reference_catalog" => "",
        ),
        "graph" => Dict{String,Any}(
            "neighbors" => 20,
            "maximum_degree" => 30,
            "metric" => "xyz_scaled",
            "depth_scale" => 0.5,
            "mutual" => true,
            "maximum_distance_km" => 0.0,
            "minimum_degree" => 0,
            "ensure_connected" => false,
            "maximum_bridge_distance_km" => 0.0,
            "augmentation" => Dict{String,Any}(
                "enabled" => false,
                "candidate_neighbors" => 40,
                "add_edges" => 0,
                "add_edges_per_event" => 0.0,
                "maximum_added_per_event" => 3,
                "power_iterations" => 30,
            ),
        ),
        "output" => Dict{String,Any}(
            "write_filtered_catalogs" => true,
            "write_graph_metadata" => true,
            "write_solver_history" => true,
            "write_run_summary" => true,
        ),
        "uncertainty" => Dict{String,Any}(
            "method" => "none",
            "linearized" => Dict{String,Any}(
                "probes" => 12,
                "seed" => 24680,
                "inner_tolerance" => 1.0e-3,
                "inner_max_iterations" => 150,
            ),
            "bootstrap" => Dict{String,Any}(
                "replicates" => 100,
                "resampling_unit" => "station_phase",
                "seed" => 12345,
                "summary_method" => "percentile",
                "confidence_level" => 0.95,
                "standard_deviation_multiplier" => 2.0,
                "write_samples" => true,
                "write_catalogs" => false,
            ),
        ),
        "experimental" => Dict{String,Any}(
            "bias" => Dict{String,Any}(
                "enabled" => false,
                "apply_model_file" => "",
                "fit_dimensions" => "xy",
                "minimum_observations" => 30,
                "huber_k" => 1.345,
                "minimum_sigma_s" => 0.002,
                "maximum_irls_iterations" => 6,
                "minimum_improvement_fraction" => 0.05,
            ),
        ),
    )
end

function deep_merge!(dst::Dict{String,Any}, src::AbstractDict)
    for (key0, value) in src
        key = String(key0)
        if value isa AbstractDict && get(dst, key, nothing) isa Dict{String,Any}
            deep_merge!(dst[key], value)
        elseif value isa AbstractDict
            child = Dict{String,Any}()
            deep_merge!(child, value)
            dst[key] = child
        else
            dst[key] = value
        end
    end
    return dst
end

function warn_unknown_keys(src::AbstractDict, defaults::AbstractDict, prefix::String="")
    for (key0, value) in src
        key = String(key0)
        full = isempty(prefix) ? key : string(prefix, ".", key)
        if !haskey(defaults, key)
            @warn "Unknown configuration key" key=full
        elseif value isa AbstractDict && defaults[key] isa AbstractDict
            warn_unknown_keys(value, defaults[key], full)
        end
    end
end

function cfgget(cfg::AbstractDict, keys...; default=nothing)
    node = cfg
    for key0 in keys
        key = String(key0)
        if !(node isa AbstractDict) || !haskey(node, key)
            return default
        end
        node = node[key]
    end
    return node
end

function normalize_config_paths!(cfg::Dict{String,Any}, base::String)
    for (section, keys) in (
        ("io", ("catalog_file", "stations_file", "theta_dir", "thetastd_dir", "output_dir")),
        ("travel_time", ("table_file", "velocity_model_file")),
        ("gauge", ("pin_reference_catalog",)),
        ("experimental", ()),
    )
        node = cfg[section]
        for key in keys
            value = String(get(node, key, ""))
            if !isempty(value) && !isabspath(value)
                node[key] = normpath(joinpath(base, value))
            end
        end
    end
    bias = cfg["experimental"]["bias"]
    value = String(get(bias, "apply_model_file", ""))
    if !isempty(value) && !isabspath(value)
        bias["apply_model_file"] = normpath(joinpath(base, value))
    end
    return cfg
end

"Load a TOML file, recursively merge it into defaults, and resolve relative paths from the TOML directory."
function load_config(path::AbstractString)
    isfile(path) || error("Configuration file not found: $path")
    supplied = TOML.parsefile(path)
    defaults = default_config()
    warn_unknown_keys(supplied, defaults)
    cfg = deepcopy(defaults)
    deep_merge!(cfg, supplied)
    normalize_config_paths!(cfg, dirname(abspath(path)))
    validate_config(cfg)
    cfg["_config_file"] = abspath(path)
    return cfg
end

function validate_config(cfg::Dict{String,Any})
    lowercase(String(cfgget(cfg, "travel_time", "type"))) in ("lookup", "constant", "constant_velocity", "constvel") ||
        error("travel_time.type must be lookup or constant")
    if lowercase(String(cfgget(cfg, "travel_time", "type"))) == "lookup"
        endswith(lowercase(String(cfgget(cfg, "travel_time", "table_file"))), ".gstt") ||
            error("travel_time.table_file must use the native .gstt extension")
    end
    geom = lowercase(String(cfgget(cfg, "travel_time", "geometry")))
    geom in ("auto", "cartesian", "flat", "radial", "spherical") ||
        error("travel_time.geometry must be auto, cartesian, or radial")
    build_geom = lowercase(String(cfgget(cfg, "travel_time", "build_geometry")))
    build_geom in ("cartesian", "flat", "radial", "spherical") ||
        error("travel_time.build_geometry must be cartesian or radial")
    metric = lowercase(String(cfgget(cfg, "graph", "metric")))
    metric in ("xy", "xyz", "xyz_scaled") || error("graph.metric must be xy, xyz, or xyz_scaled")
    Int(cfgget(cfg, "graph", "neighbors")) >= 1 || error("graph.neighbors must be at least 1")
    isfinite(Float64(cfgget(cfg, "graph", "depth_scale"))) && Float64(cfgget(cfg, "graph", "depth_scale")) > 0.0 ||
        error("graph.depth_scale must be finite and positive")
    reference = lowercase(String(cfgget(cfg, "coordinates", "reference")))
    reference in ("stations_mean", "catalog_mean", "manual") ||
        error("coordinates.reference must be stations_mean, catalog_mean, or manual")
    initialization = lowercase(String(cfgget(cfg, "initialization", "mode")))
    initialization in ("catalog", "common_centroid", "common_manual") ||
        error("initialization.mode must be catalog, common_centroid, or common_manual")
    if initialization != "catalog" && !Bool(cfgget(cfg, "run", "prelocation"; default=true))
        error("initialization.mode=$initialization requires run.prelocation=true; Stage 1 must separate the common seed before the Stage-2 event graph is built")
    end
    initial_latitude = Float64(cfgget(cfg, "initialization", "latitude"))
    initial_longitude = Float64(cfgget(cfg, "initialization", "longitude"))
    initial_depth = Float64(cfgget(cfg, "initialization", "depth_km"))
    isfinite(initial_latitude) && -90.0 <= initial_latitude <= 90.0 ||
        error("initialization.latitude must be finite and between -90 and 90 degrees")
    isfinite(initial_longitude) && -180.0 <= initial_longitude <= 180.0 ||
        error("initialization.longitude must be finite and between -180 and 180 degrees")
    isfinite(initial_depth) || error("initialization.depth_km must be finite")
    Float64(cfgget(cfg, "coordinates", "earth_radius_m")) > 0.0 || error("coordinates.earth_radius_m must be positive")
    lowercase(String(cfgget(cfg, "coordinates", "station_vertical"))) in ("depth_from_elevation", "elevation") ||
        error("coordinates.station_vertical must be depth_from_elevation or elevation")
    lowercase(String(cfgget(cfg, "coordinates", "event_vertical"))) in
        ("positive_depth", "negative_depth", "positive_depth_plus_z0") ||
        error("coordinates.event_vertical is invalid")
    travel_radius = Float64(cfgget(cfg, "travel_time", "earth_radius_m"))
    isfinite(travel_radius) && travel_radius >= 0.0 || error("travel_time.earth_radius_m must be zero (auto) or positive")
    for key in ("horizontal_step_m", "depth_step_m", "station_depth_step_m")
        value = Float64(cfgget(cfg, "lookup", key))
        isfinite(value) && value > 0.0 || error("lookup.$key must be finite and positive")
    end
    Float64(cfgget(cfg, "lookup", "station_depth_step_m")) >= Float64(cfgget(cfg, "lookup", "depth_step_m")) ||
        error("lookup.station_depth_step_m must be at least lookup.depth_step_m")
    gauge_mode = lowercase(String(cfgget(cfg, "gauge", "mode")))
    gauge_mode in ("zero_mean", "pin") || error("gauge.mode must be zero_mean or pin")
    Float64(cfgget(cfg, "gauge", "constraint_weight")) > 0.0 || error("gauge.constraint_weight must be positive")
    augmentation = cfgget(cfg, "graph", "augmentation")
    candidate_neighbors = Int(get(augmentation, "candidate_neighbors", 40))
    add_edges = Int(get(augmentation, "add_edges", 0))
    add_edges_per_event = Float64(get(augmentation, "add_edges_per_event", 0.0))
    maximum_added_per_event = Int(get(augmentation, "maximum_added_per_event", 3))
    power_iterations = Int(get(augmentation, "power_iterations", 30))
    candidate_neighbors >= 1 || error("graph.augmentation.candidate_neighbors must be at least 1")
    add_edges >= 0 || error("graph.augmentation.add_edges cannot be negative")
    isfinite(add_edges_per_event) && add_edges_per_event >= 0.0 ||
        error("graph.augmentation.add_edges_per_event must be finite and nonnegative")
    maximum_added_per_event >= 1 || error("graph.augmentation.maximum_added_per_event must be at least 1")
    power_iterations >= 1 || error("graph.augmentation.power_iterations must be at least 1")
    if Bool(get(augmentation, "enabled", false))
        (add_edges > 0) != (add_edges_per_event > 0.0) ||
            error("Enabled graph augmentation requires exactly one positive target: graph.augmentation.add_edges or add_edges_per_event")
    end
    for section in ("prelocation", "relocation")
        lowercase(String(cfgget(cfg, section, "linear_solver"))) in ("pcg", "direct") ||
            error("$section.linear_solver must be pcg or direct")
        lowercase(String(cfgget(cfg, section, "preconditioner"))) in ("block_jacobi", "diagonal", "none") ||
            error("$section.preconditioner must be block_jacobi, diagonal, or none")
        damping = Float64(cfgget(cfg, section, "damping_lambda"))
        isfinite(damping) && damping >= 0.0 || error("$section.damping_lambda must be finite and nonnegative")
        Float64(cfgget(cfg, section, "min_sigma_s")) > 0.0 || error("$section.min_sigma_s must be positive")
        0.0 < Float64(cfgget(cfg, section, "step_damping")) <= 1.0 ||
            error("$section.step_damping must be in (0, 1]")
    end
    uncertainty_method = lowercase(String(cfgget(cfg, "uncertainty", "method")))
    uncertainty_method in ("none", "linearized", "bootstrap", "both") ||
        error("uncertainty.method must be none, linearized, bootstrap, or both")
    linearized = cfgget(cfg, "uncertainty", "linearized")
    Int(get(linearized, "probes", 12)) >= 1 || error("uncertainty.linearized.probes must be at least 1")
    Float64(get(linearized, "inner_tolerance", 1.0e-3)) > 0.0 ||
        error("uncertainty.linearized.inner_tolerance must be positive")
    Int(get(linearized, "inner_max_iterations", 150)) >= 1 ||
        error("uncertainty.linearized.inner_max_iterations must be at least 1")
    bootstrap = cfgget(cfg, "uncertainty", "bootstrap")
    Int(get(bootstrap, "replicates", 100)) >= 1 || error("uncertainty.bootstrap.replicates must be at least 1")
    lowercase(String(get(bootstrap, "resampling_unit", "station_phase"))) in ("station_phase", "station") ||
        error("uncertainty.bootstrap.resampling_unit must be station_phase or station")
    summary_method = lowercase(String(get(bootstrap, "summary_method", "percentile")))
    summary_method in ("percentile", "standard_deviation", "std", "2std") ||
        error("uncertainty.bootstrap.summary_method must be percentile or standard_deviation")
    confidence_level = Float64(get(bootstrap, "confidence_level", 0.95))
    0.0 < confidence_level < 1.0 || error("uncertainty.bootstrap.confidence_level must be between 0 and 1")
    Float64(get(bootstrap, "standard_deviation_multiplier", 2.0)) > 0.0 ||
        error("uncertainty.bootstrap.standard_deviation_multiplier must be positive")
    lowercase(String(cfgget(cfg, "experimental", "bias", "fit_dimensions"))) in ("xy", "xyz") ||
        error("experimental.bias.fit_dimensions must be xy or xyz")
    return cfg
end
