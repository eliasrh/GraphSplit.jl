"Build the final robust, damped normal-equation operator used by uncertainty diagnostics."
function final_normal_equations(state::State, stations::Stations, obs::Observations,
        travel_time::AbstractTravelTimeModel, solver_cfg::AbstractDict, cfg::AbstractDict;
        depth_bound_active::BitVector=falses(length(state)))
    prediction, gix, giy, giz, gjx, gjy, gjz = predict_and_gradients(state, stations, obs, travel_time)
    residual = obs.dt .- prediction
    minimum_sigma = Float64(get(solver_cfg, "min_sigma_s", 0.002))
    huber_k = Float64(get(solver_cfg, "huber_k", 1.345))
    sigma = max.(obs.sigma, minimum_sigma)
    robust = [huber_weight(residual[row] / sigma[row], huber_k) for row in eachindex(residual)]
    sqrt_weight = sqrt.(robust ./ (sigma .^ 2))
    n = length(state)
    system = LinearizedSystem(obs.i, obs.j, gix, giy, giz, gjx, gjy, gjz, sqrt_weight, n,
        gauge_definition(cfg), Float64(cfgget(cfg, "gauge", "constraint_weight"; default=10.0)),
        Float64(cfgget(cfg, "gauge", "reference_velocity_ms"; default=5000.0)))
    length(depth_bound_active) == n || error("Depth-bound active-mask length mismatch")
    active = parameter_active_mask(state, cfg)
    for event in findall(depth_bound_active)
        active[2n + event] = false
    end
    damping = max(Float64(get(solver_cfg, "damping_lambda", 0.0)), 1.0e-8)
    operator = function (vector::Vector{Float64})
        work = copy(vector)
        work[.!active] .= 0.0
        result = apply_At(system, apply_A(system, work)) .+ damping .* work
        result[.!active] .= vector[.!active]
        return result
    end
    preconditioner_name = lowercase(String(get(solver_cfg, "preconditioner", "block_jacobi")))
    apply_prec = preconditioner(system, damping, active, preconditioner_name)
    return system, operator, apply_prec, active
end

"Estimate spatial blocks of the regularized inverse Hessian with block Hutchinson probes."
function estimate_linearized_uncertainty(state::State, stations::Stations, obs::Observations,
        travel_time::AbstractTravelTimeModel, event_mask::AbstractVector{Bool}, cfg::AbstractDict;
        depth_bound_active::BitVector=falses(length(state)))
    length(event_mask) == length(state) || error("Linearized uncertainty event-mask length mismatch")
    options = cfgget(cfg, "uncertainty", "linearized")
    probes = Int(get(options, "probes", 12))
    seed = Int(get(options, "seed", 24680))
    tolerance = Float64(get(options, "inner_tolerance", 1.0e-3))
    maximum_iterations = Int(get(options, "inner_max_iterations", 150))
    system, operator, apply_prec, active = final_normal_equations(state, stations, obs,
        travel_time, cfg["relocation"], cfg; depth_bound_active=depth_bound_active)
    n = length(state)
    accumulated = zeros(3, 3, n)
    successful_by_axis = zeros(Int, 3)
    rng = MersenneTwister(seed)
    rhs = zeros(4n)
    signs = zeros(Float64, n)
    for probe in 1:probes, column in 1:3
        fill!(rhs, 0.0)
        for event in 1:n
            signs[event] = rand(rng, Bool) ? 1.0 : -1.0
            index = (column - 1) * n + event
            active[index] && (rhs[index] = signs[event])
        end
        solution, _, converged = pcg(operator, rhs, apply_prec;
            tolerance=tolerance, maximum_iterations=maximum_iterations)
        converged || continue
        successful_by_axis[column] += 1
        for event in 1:n, row in 1:3
            accumulated[row, column, event] += solution[(row - 1) * n + event] * signs[event]
        end
    end
    covariance = fill(NaN, 3, 3, n)
    if any(==(0), successful_by_axis)
        @warn "Some linearized uncertainty probe directions had no converged solves" successful_by_axis
    end
    for event in 1:n
        event_mask[event] || continue
        raw = Matrix{Float64}(undef, 3, 3)
        for row in 1:3, column in 1:3
            count = successful_by_axis[column]
            raw[row, column] = count > 0 ? accumulated[row, column, event] / count : NaN
        end
        all(isfinite, raw) || continue
        for dimension in 1:3
            if !active[(dimension - 1) * n + event]
                raw[dimension, :] .= 0.0
                raw[:, dimension] .= 0.0
            end
        end
        block = 0.5 .* (raw .+ raw')
        decomposition = eigen(Symmetric(block))
        values = max.(decomposition.values, 0.0)
        covariance[:, :, event] .= decomposition.vectors * Diagonal(values) * decomposition.vectors'
        if depth_bound_active[event]
            covariance[3, :, event] .= NaN
            covariance[:, 3, event] .= NaN
        end
    end
    @printf("Linearized uncertainty: %d/%d randomized solves converged\n",
        sum(successful_by_axis), 3probes)
    return (covariance=covariance, probes=probes,
        successful_solves=sum(successful_by_axis), attempted_solves=3probes,
        successful_by_axis=successful_by_axis)
end

function write_linearized_uncertainty(path::AbstractString, event_id::Vector{Int64},
        event_mask::AbstractVector{Bool}, estimate)
    covariance = estimate.covariance
    open(path, "w") do io
        println(io, "EventID std_x_m std_y_m std_z_m cov_xx_m2 cov_xy_m2 cov_xz_m2 cov_yy_m2 cov_yz_m2 cov_zz_m2")
        for event in eachindex(event_id)
            event_mask[event] || continue
            block = @view covariance[:, :, event]
            stdx = isfinite(block[1, 1]) ? sqrt(max(block[1, 1], 0.0)) : NaN
            stdy = isfinite(block[2, 2]) ? sqrt(max(block[2, 2], 0.0)) : NaN
            stdz = isfinite(block[3, 3]) ? sqrt(max(block[3, 3], 0.0)) : NaN
            @printf(io, "%d %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n",
                event_id[event], stdx, stdy, stdz, block[1, 1], block[1, 2], block[1, 3],
                block[2, 2], block[2, 3], block[3, 3])
        end
    end
    return path
end

"Map each retained observation to its station-phase or station resampling block."
function bootstrap_blocks(obs::Observations, groups::Vector{ThetaGroup}, unit::String)
    used_groups = sort!(unique(Int.(obs.group)))
    isempty(used_groups) && error("No retained station-phase groups are available for bootstrap resampling")
    group_to_block = Dict{Int,Int}()
    labels = String[]
    if unit == "station_phase"
        for group_index in used_groups
            1 <= group_index <= length(groups) || error("Observation contains an invalid theta-group index")
            push!(labels, groups[group_index].name)
            group_to_block[group_index] = length(labels)
        end
    elseif unit == "station"
        station_to_block = Dict{String,Int}()
        for group_index in used_groups
            1 <= group_index <= length(groups) || error("Observation contains an invalid theta-group index")
            station = groups[group_index].station
            block = get(station_to_block, station, 0)
            if block == 0
                push!(labels, station)
                block = length(labels)
                station_to_block[station] = block
            end
            group_to_block[group_index] = block
        end
    else
        error("Unknown bootstrap resampling unit: $unit")
    end
    block_index = Vector{Int32}(undef, length(obs))
    for row in eachindex(obs.group)
        block_index[row] = Int32(group_to_block[Int(obs.group[row])])
    end
    return block_index, labels
end

function draw_bootstrap_counts(rng::AbstractRNG, number_of_blocks::Int)
    counts = zeros(Int, number_of_blocks)
    for _ in 1:number_of_blocks
        counts[rand(rng, 1:number_of_blocks)] += 1
    end
    return counts
end

"Apply integer bootstrap multiplicities without altering robust sigma normalization."
function resample_observations(obs::Observations, block_index::AbstractVector{<:Integer}, counts::Vector{Int})
    length(block_index) == length(obs) || error("Bootstrap block-index length mismatch")
    total = 0
    for row in eachindex(block_index)
        total += counts[block_index[row]]
    end
    result = Observations()
    for values in (result.i, result.j, result.dt, result.sigma, result.station, result.phase, result.group)
        sizehint!(values, total)
    end
    for row in eachindex(obs.i)
        for _ in 1:counts[block_index[row]]
            push_observation!(result, obs.i[row], obs.j[row], obs.dt[row], obs.sigma[row],
                obs.station[row], obs.phase[row], obs.group[row])
        end
    end
    return result
end

function fully_spatially_pinned(state::State, cfg::AbstractDict)
    pinned = falses(length(state))
    lowercase(String(cfgget(cfg, "gauge", "mode"; default="zero_mean"))) == "pin" || return pinned
    fields = lowercase(String(cfgget(cfg, "gauge", "pin_fields"; default="xyz")))
    all(field -> occursin(field, fields), ('x', 'y', 'z')) || return pinned
    rows = Dict(id => row for (row, id) in enumerate(state.event_id))
    for id in Int64.(cfgget(cfg, "gauge", "pin_event_ids"; default=Int[]))
        haskey(rows, id) && (pinned[rows[id]] = true)
    end
    return pinned
end

function create_bootstrap_store(output::AbstractString, rows::Int, replicates::Int)
    path = joinpath(output, ".bootstrap_samples_work.bin")
    io = open(path, "w+")
    values_per_field = Base.Checked.checked_mul(rows, replicates)
    bytes_per_field = Base.Checked.checked_mul(sizeof(Float64), values_per_field)
    total_bytes = Base.Checked.checked_mul(4, bytes_per_field)
    seek(io, total_bytes - 1)
    write(io, UInt8(0))
    flush(io)
    dims = (rows, replicates)
    x = Mmap.mmap(io, Array{Float64,2}, dims, 0)
    y = Mmap.mmap(io, Array{Float64,2}, dims, bytes_per_field)
    z = Mmap.mmap(io, Array{Float64,2}, dims, 2bytes_per_field)
    t0 = Mmap.mmap(io, Array{Float64,2}, dims, 3bytes_per_field)
    return (path=path, io=io, x=x, y=y, z=z, t0=t0)
end

function write_bootstrap_sample_tables(output::AbstractString, catalog::Catalog, nominal::State,
        selected_rows::Vector{Int}, samples, replicates::Int, cfg::AbstractDict)
    nominal_lat, nominal_lon, nominal_depth = catalog_with_state(catalog, nominal, cfg)
    paths = [joinpath(output, "bootstrap_samples_lon.txt"),
        joinpath(output, "bootstrap_samples_lat.txt"),
        joinpath(output, "bootstrap_samples_depth_km.txt"),
        joinpath(output, "bootstrap_samples_t0_s.txt")]
    ios = [open(path, "w") for path in paths]
    try
        for io in ios
            print(io, "EventID all_data")
            for replicate in 1:replicates
                @printf(io, " sample_%04d", replicate)
            end
            println(io)
        end
        for (sample_row, event) in enumerate(selected_rows)
            id = catalog.event_id[event]
            @printf(ios[1], "%d %.10f", id, nominal_lon[event])
            @printf(ios[2], "%d %.10f", id, nominal_lat[event])
            @printf(ios[3], "%d %.8f", id, nominal_depth[event])
            @printf(ios[4], "%d %.12g", id, nominal.t0[event])
            for replicate in 1:replicates
                x = samples.x[sample_row, replicate]
                y = samples.y[sample_row, replicate]
                z = samples.z[sample_row, replicate]
                if isfinite(x) && isfinite(y)
                    lat, lon = local_xy_to_ll(x, y, nominal.ref_lat, nominal.ref_lon, nominal.ref_radius_m)
                    @printf(ios[1], " %.10f", lon)
                    @printf(ios[2], " %.10f", lat)
                else
                    print(ios[1], " NaN")
                    print(ios[2], " NaN")
                end
                isfinite(z) ? @printf(ios[3], " %.8f", internal_depth_to_km(z, cfg)) : print(ios[3], " NaN")
                t0 = samples.t0[sample_row, replicate]
                isfinite(t0) ? @printf(ios[4], " %.12g", t0) : print(ios[4], " NaN")
            end
            for io in ios
                println(io)
            end
        end
    finally
        foreach(close, ios)
    end
    return paths
end

function sample_quantile(values::Vector{Float64}, probability::Float64)
    isempty(values) && return NaN
    sorted = sort(values)
    length(sorted) == 1 && return sorted[1]
    position = 1.0 + (length(sorted) - 1) * clamp(probability, 0.0, 1.0)
    lower = floor(Int, position)
    upper = ceil(Int, position)
    lower == upper && return sorted[lower]
    fraction = position - lower
    return (1.0 - fraction) * sorted[lower] + fraction * sorted[upper]
end

function sample_covariance3(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64})
    n = length(x)
    n >= 2 || return fill(NaN, 3, 3)
    means = (sum(x) / n, sum(y) / n, sum(z) / n)
    covariance = zeros(3, 3)
    for k in 1:n
        delta = (x[k] - means[1], y[k] - means[2], z[k] - means[3])
        for row in 1:3, column in row:3
            covariance[row, column] += delta[row] * delta[column]
        end
    end
    covariance ./= n - 1
    covariance[2, 1] = covariance[1, 2]
    covariance[3, 1] = covariance[1, 3]
    covariance[3, 2] = covariance[2, 3]
    return covariance
end

function bootstrap_axis_summary(values::Vector{Float64}, method::String,
        confidence_level::Float64, standard_deviation_multiplier::Float64)
    isempty(values) && return (NaN, NaN, NaN)
    if method == "percentile"
        tail = (1.0 - confidence_level) / 2.0
        return sample_quantile(values, 0.5), sample_quantile(values, tail),
            sample_quantile(values, 1.0 - tail)
    end
    center = sum(values) / length(values)
    standard_deviation = length(values) >= 2 ? sqrt(sum((value - center)^2 for value in values) /
        (length(values) - 1)) : NaN
    half_width = standard_deviation_multiplier * standard_deviation
    return center, center - half_width, center + half_width
end

function write_bootstrap_summary(path::AbstractString, event_id::Vector{Int64}, selected_rows::Vector{Int},
        nominal::State, samples, replicates::Int, cfg::AbstractDict)
    options = cfgget(cfg, "uncertainty", "bootstrap")
    method0 = lowercase(String(get(options, "summary_method", "percentile")))
    method = method0 == "percentile" ? "percentile" : "standard_deviation"
    confidence = Float64(get(options, "confidence_level", 0.95))
    multiplier = Float64(get(options, "standard_deviation_multiplier", 2.0))
    open(path, "w") do io
        println(io, "EventID n_valid valid_fraction center_dx_m lower_dx_m upper_dx_m center_dy_m lower_dy_m upper_dy_m center_dz_m lower_dz_m upper_dz_m cov_xx_m2 cov_xy_m2 cov_xz_m2 cov_yy_m2 cov_yz_m2 cov_zz_m2")
        for (sample_row, event) in enumerate(selected_rows)
            dx, dy, dz = Float64[], Float64[], Float64[]
            for replicate in 1:replicates
                x = samples.x[sample_row, replicate]
                y = samples.y[sample_row, replicate]
                z = samples.z[sample_row, replicate]
                isfinite(x) && isfinite(y) && isfinite(z) || continue
                push!(dx, x - nominal.x[event])
                push!(dy, y - nominal.y[event])
                push!(dz, z - nominal.z[event])
            end
            sx = bootstrap_axis_summary(dx, method, confidence, multiplier)
            sy = bootstrap_axis_summary(dy, method, confidence, multiplier)
            sz = bootstrap_axis_summary(dz, method, confidence, multiplier)
            covariance = sample_covariance3(dx, dy, dz)
            @printf(io, "%d %d %.8f %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n",
                event_id[event], length(dx), length(dx) / replicates,
                sx[1], sx[2], sx[3], sy[1], sy[2], sy[3], sz[1], sz[2], sz[3],
                covariance[1, 1], covariance[1, 2], covariance[1, 3],
                covariance[2, 2], covariance[2, 3], covariance[3, 3])
        end
    end
    return path
end

"Run a conditional Stage-2 block bootstrap and write samples separately from catalogs."
function run_bootstrap_uncertainty(catalog::Catalog, groups::Vector{ThetaGroup}, nominal::State,
        stations::Stations, obs::Observations, travel_time::AbstractTravelTimeModel,
        event_mask::AbstractVector{Bool}, output::AbstractString, cfg::AbstractDict)
    options = cfgget(cfg, "uncertainty", "bootstrap")
    replicates = Int(get(options, "replicates", 100))
    unit = lowercase(String(get(options, "resampling_unit", "station_phase")))
    seed = Int(get(options, "seed", 12345))
    write_samples = Bool(get(options, "write_samples", true))
    write_catalogs = Bool(get(options, "write_catalogs", false))
    block_index, labels = bootstrap_blocks(obs, groups, unit)
    selected_rows = findall(event_mask)
    samples = create_bootstrap_store(output, length(selected_rows), replicates)
    rng = MersenneTwister(seed)
    solver_cfg = deepcopy(cfg["relocation"])
    solver_cfg["verbose"] = false
    pinned = fully_spatially_pinned(nominal, cfg)
    replicate_converged = falses(replicates)
    replicate_iterations = zeros(Int, replicates)
    replicate_active = zeros(Int, replicates)
    replicate_unique_blocks = zeros(Int, replicates)
    valid_depth_samples = zeros(Int, length(selected_rows))
    active_depth_bound_samples = zeros(Int, length(selected_rows))
    block_counts = zeros(Int32, length(labels), replicates)
    catalog_directory = joinpath(output, "bootstrap_catalogs")
    write_catalogs && mkpath(catalog_directory)
    try
        for replicate in 1:replicates
            @views samples.x[:, replicate] .= NaN
            @views samples.y[:, replicate] .= NaN
            @views samples.z[:, replicate] .= NaN
            @views samples.t0[:, replicate] .= NaN
            counts = draw_bootstrap_counts(rng, length(labels))
            block_counts[:, replicate] .= counts
            replicate_unique_blocks[replicate] = count(>(0), counts)
            sampled_obs = resample_observations(obs, block_index, counts)
            active_events = event_activity(sampled_obs, length(catalog))
            active_events .|= pinned
            replicate_active[replicate] = count(active_events .& event_mask)
            state, stats = solve_relocation(nominal, stations, sampled_obs, travel_time, solver_cfg, cfg)
            finite_state = all(isfinite, state.x) && all(isfinite, state.y) && all(isfinite, state.z) && all(isfinite, state.t0)
            replicate_converged[replicate] = stats.converged && finite_state
            replicate_iterations[replicate] = stats.iterations
            if replicate_converged[replicate]
                for (sample_row, event) in enumerate(selected_rows)
                    active_events[event] || continue
                    samples.x[sample_row, replicate] = state.x[event]
                    samples.y[sample_row, replicate] = state.y[event]
                    samples.z[sample_row, replicate] = state.z[event]
                    samples.t0[sample_row, replicate] = state.t0[event]
                    valid_depth_samples[sample_row] += 1
                    stats.depth_bound_active[event] && (active_depth_bound_samples[sample_row] += 1)
                end
                if write_catalogs
                    mask = event_mask .& active_events
                    write_catalog(joinpath(catalog_directory, @sprintf("catalog_dd_filt_boot_%04d.txt", replicate)),
                        catalog, state, cfg; mask=mask)
                end
            end
            @printf("  bootstrap %4d/%d: blocks %d/%d; active events %d; %s\n", replicate, replicates,
                replicate_unique_blocks[replicate], length(labels), replicate_active[replicate],
                replicate_converged[replicate] ? "converged" : "not converged")
        end
        Mmap.sync!(samples.x); Mmap.sync!(samples.y); Mmap.sync!(samples.z); Mmap.sync!(samples.t0)
        write_bootstrap_summary(joinpath(output, "booterrxyz.txt"), catalog.event_id, selected_rows,
            nominal, samples, replicates, cfg)
        write_samples && write_bootstrap_sample_tables(output, catalog, nominal, selected_rows,
            samples, replicates, cfg)
        open(joinpath(output, "bootstrap_replicates.txt"), "w") do io
            println(io, "replicate converged iterations unique_blocks total_blocks active_filtered_events")
            for replicate in 1:replicates
                @printf(io, "%d %d %d %d %d %d\n", replicate, Int(replicate_converged[replicate]),
                    replicate_iterations[replicate], replicate_unique_blocks[replicate], length(labels),
                    replicate_active[replicate])
            end
        end
        open(joinpath(output, "bootstrap_block_counts.txt"), "w") do io
            print(io, "block label")
            for replicate in 1:replicates
                @printf(io, " sample_%04d", replicate)
            end
            println(io)
            for block in eachindex(labels)
                @printf(io, "%d %s", block, labels[block])
                for replicate in 1:replicates
                    @printf(io, " %d", block_counts[block, replicate])
                end
                println(io)
            end
        end
        if Bool(get(depth_bound_options(cfg), "enabled", false))
            open(joinpath(output, "bootstrap_depth_bound_status.txt"), "w") do io
                println(io, "EventID n_valid n_bound_active bound_active_fraction")
                for (sample_row, event) in enumerate(selected_rows)
                    valid = valid_depth_samples[sample_row]
                    fraction = valid > 0 ? active_depth_bound_samples[sample_row] / valid : NaN
                    @printf(io, "%d %d %d %.12g\n", catalog.event_id[event], valid,
                        active_depth_bound_samples[sample_row], fraction)
                end
            end
        end
        metadata = Dict{String,Any}(
            "method" => "conditional_stage2_block_bootstrap",
            "resampling_unit" => unit,
            "blocks" => length(labels),
            "replicates" => replicates,
            "converged_replicates" => count(identity, replicate_converged),
            "seed" => seed,
            "fixed_nominal_graph" => true,
            "prelocation_resampled" => false,
            "write_samples" => write_samples,
            "write_catalogs" => write_catalogs,
            "summary_method" => String(get(options, "summary_method", "percentile")),
            "confidence_level" => Float64(get(options, "confidence_level", 0.95)),
            "standard_deviation_multiplier" => Float64(get(options, "standard_deviation_multiplier", 2.0)),
            "depth_bound_enabled" => Bool(get(depth_bound_options(cfg), "enabled", false)),
        )
        open(joinpath(output, "bootstrap_metadata.toml"), "w") do io
            TOML.print(io, metadata; sorted=true)
        end
    finally
        close(samples.io)
        try
            rm(samples.path; force=true)
        catch exception
            @warn "Could not remove temporary bootstrap sample store" path=samples.path exception=exception
        end
    end
    return (replicates=replicates, converged=count(identity, replicate_converged), blocks=length(labels),
        resampling_unit=unit, samples_written=write_samples, catalogs_written=write_catalogs)
end
