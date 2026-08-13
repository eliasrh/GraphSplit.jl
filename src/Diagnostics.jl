function weighted_robust_rms(residual::Vector{Float64}, sigma::Vector{Float64}, huber_k::Float64)
    weights = [huber_weight(residual[i] / sigma[i], huber_k) / sigma[i]^2 for i in eachindex(residual)]
    return sqrt(sum(weights .* residual .^ 2) / max(sum(weights), eps(Float64)))
end

function robust_spatial_fit(matrix::Matrix{Float64}, values::Vector{Float64}, sigma::Vector{Float64},
        huber_k::Float64, maximum_iterations::Int)
    parameters = zeros(size(matrix, 2))
    for _ in 1:maximum_iterations
        normalized = (values .- matrix * parameters) ./ sigma
        weights = [huber_weight(value, huber_k) / sigma[i]^2 for (i, value) in enumerate(normalized)]
        weighted = matrix .* sqrt.(weights)
        rhs = values .* sqrt.(weights)
        hessian = weighted' * weighted + 1.0e-10I
        candidate = hessian \ (weighted' * rhs)
        norm(candidate - parameters) <= 1.0e-12 * (1.0 + norm(parameters)) && return candidate
        parameters = candidate
    end
    return parameters
end

"Fit the experimental per-station-phase coherent theta bias model."
function scan_theta_bias(groups::Vector{ThetaGroup}, state::State, stations::Stations,
        catalog::Catalog, travel_time::AbstractTravelTimeModel, cfg::AbstractDict)
    options = cfgget(cfg, "experimental", "bias")
    station_rows = station_row_map(stations)
    event_rows = event_row_map(catalog)
    use_z = lowercase(String(get(options, "fit_dimensions", "xy"))) == "xyz"
    minimum_observations = Int(get(options, "minimum_observations", 30))
    minimum_sigma = Float64(get(options, "minimum_sigma_s", 0.002))
    minimum_degree = Int(cfgget(cfg, "observations", "minimum_theta_degree"; default=6))
    huber_k = Float64(get(options, "huber_k", 1.345))
    maximum_iterations = Int(get(options, "maximum_irls_iterations", 6))
    minimum_improvement = Float64(get(options, "minimum_improvement_fraction", 0.05))
    result = NamedTuple[]
    for (group_index, group) in enumerate(groups)
        haskey(station_rows, group.station) || continue
        station = station_rows[group.station]
        rows_i, rows_j = Int32[], Int32[]
        theta_difference, sigma = Float64[], Float64[]
        for event_id in group.event_id
            first = group.entries[event_id]
            event_id == first.ref_id && continue
            haskey(event_rows, event_id) && haskey(event_rows, first.ref_id) && haskey(group.entries, first.ref_id) || continue
            group.has_degree && minimum_degree > 0 && first.degree < minimum_degree && continue
            second = group.entries[first.ref_id]
            push!(rows_i, event_rows[event_id]); push!(rows_j, event_rows[first.ref_id])
            push!(theta_difference, first.theta - second.theta)
            push!(sigma, max(combined_sigma(first, second, minimum_sigma, group.has_sigma), minimum_sigma))
        end
        length(rows_i) >= minimum_observations || continue
        residual = Vector{Float64}(undef, length(rows_i))
        design = Matrix{Float64}(undef, length(rows_i), use_z ? 3 : 2)
        for row in eachindex(rows_i)
            i, j = rows_i[row], rows_j[row]
            ti = travel_time_gradient(travel_time, group.phase, state.x[i], state.y[i], state.z[i],
                stations.x_m[station], stations.y_m[station], stations.z_m[station])[1]
            tj = travel_time_gradient(travel_time, group.phase, state.x[j], state.y[j], state.z[j],
                stations.x_m[station], stations.y_m[station], stations.z_m[station])[1]
            predicted = state.t0[i] - state.t0[j] + ti - tj
            residual[row] = theta_difference[row] - predicted
            design[row, 1] = state.x[i] - state.x[j]
            design[row, 2] = state.y[i] - state.y[j]
            use_z && (design[row, 3] = state.z[i] - state.z[j])
        end
        coefficients_fit = robust_spatial_fit(design, residual, sigma, huber_k, maximum_iterations)
        coefficients = zeros(3); coefficients[1:length(coefficients_fit)] .= coefficients_fit
        before = weighted_robust_rms(residual, sigma, huber_k)
        after = weighted_robust_rms(residual .- design * coefficients_fit, sigma, huber_k)
        improvement = before > 0.0 ? (before - after) / before : 0.0
        apply = improvement >= minimum_improvement
        push!(result, (group_index=group_index, station=group.station, phase=group.phase,
            observations=length(rows_i), apply=apply, ax=coefficients[1], ay=coefficients[2], az=coefficients[3],
            rms_before=before, rms_after=after, improvement=improvement))
    end
    return result
end

function write_bias_report(path::AbstractString, report)
    open(path, "w") do io
        println(io, "group_index,station,phase,nobs,apply_bias,ax_s_per_m,ay_s_per_m,az_s_per_m,rms_before_s,rms_after_s,improvement_fraction")
        for row in report
            @printf(io, "%d,%s,%s,%d,%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
                row.group_index, row.station, row.phase == 1 ? "P" : "S", row.observations, Int(row.apply),
                row.ax, row.ay, row.az, row.rms_before, row.rms_after, row.improvement)
        end
    end
end

function read_bias_model(path::AbstractString)
    isfile(path) || error("Bias model file not found: $path")
    model = Dict{Tuple{String,UInt8},NTuple{3,Float64}}()
    for (line_number, line) in enumerate(eachline(path))
        line_number == 1 && continue
        fields = split(strip(line), ',')
        length(fields) >= 8 || continue
        apply = lowercase(fields[5]) in ("1", "true")
        apply || continue
        phase = uppercase(fields[3]) == "P" ? UInt8(1) : UInt8(2)
        coefficients = (parse(Float64, fields[6]), parse(Float64, fields[7]), parse(Float64, fields[8]))
        model[(fields[2], phase)] = coefficients
    end
    return model
end

function apply_bias_model!(groups::Vector{ThetaGroup}, path::AbstractString, state::State, catalog::Catalog)
    model = read_bias_model(path)
    rows = event_row_map(catalog)
    corrected = 0
    for group in groups
        coefficients = get(model, (group.station, group.phase), nothing)
        coefficients === nothing && continue
        for event_id in group.event_id
            haskey(rows, event_id) || continue
            row = rows[event_id]
            old = group.entries[event_id]
            correction = coefficients[1] * state.x[row] + coefficients[2] * state.y[row] + coefficients[3] * state.z[row]
            group.entries[event_id] = ThetaEntry(old.theta - correction, old.ref_id, old.sigma, old.degree)
        end
        corrected += 1
    end
    @printf("Applied bias model to %d theta groups\n", corrected)
    return groups
end

function write_solver_history(path::AbstractString, pre::SolveStats, dd::SolveStats)
    open(path, "w") do io
        println(io, "stage,outer_iteration,robust_rms_s,step_rms_m,step_rms_s,inner_iterations")
        for (name, stats) in (("prelocation", pre), ("relocation", dd))
            for iteration in 1:stats.iterations
                @printf(io, "%s,%d,%.12g,%.12g,%.12g,%d\n", name, iteration,
                    stats.rms_s[iteration], stats.step_rms_m[iteration], stats.step_rms_s[iteration], stats.inner_iterations[iteration])
            end
        end
    end
end

function write_graph_metadata(path::AbstractString, catalog::Catalog, graph::EventGraph,
        pre::Observations, dd::Observations)
    n = length(catalog)
    pre_active, dd_active = event_activity(pre, n), event_activity(dd, n)
    pre_count, dd_count = event_observation_counts(pre, n), event_observation_counts(dd, n)
    pair_degree = zeros(Int, n)
    pair_seen = Set{UInt64}()
    for row in eachindex(dd.i)
        key = edge_key(dd.i[row], dd.j[row])
        key in pair_seen && continue
        push!(pair_seen, key); pair_degree[dd.i[row]] += 1; pair_degree[dd.j[row]] += 1
    end
    open(path, "w") do io
        println(io, "EventID,component_id,component_size,geometry_degree,dd_pair_degree,prelocation_observations,dd_observations,prelocation_used,relocation_used")
        for event in 1:n
            @printf(io, "%d,%d,%d,%d,%d,%d,%d,%d,%d\n", catalog.event_id[event], graph.comp_id[event],
                graph.comp_size[event], graph.degree[event], pair_degree[event], pre_count[event], dd_count[event],
                Int(pre_active[event]), Int(dd_active[event]))
        end
    end
end
