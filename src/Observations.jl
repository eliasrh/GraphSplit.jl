function event_row_map(catalog::Catalog)
    return Dict{Int64,Int32}(id => Int32(row) for (row, id) in enumerate(catalog.event_id))
end

function station_row_map(stations::Stations)
    return Dict{String,Int32}(id => Int32(row) for (row, id) in enumerate(stations.id))
end

function combined_sigma(first::ThetaEntry, second::ThetaEntry, minimum_sigma::Float64, has_sigma::Bool)
    if has_sigma
        s1 = isfinite(first.sigma) ? abs(first.sigma) : minimum_sigma
        s2 = isfinite(second.sigma) ? abs(second.sigma) : minimum_sigma
        return hypot(s1, s2)
    end
    return sqrt(2.0) * minimum_sigma
end

function push_observation!(obs::Observations, i, j, dt, sigma, station, phase, group)
    push!(obs.i, Int32(i)); push!(obs.j, Int32(j)); push!(obs.dt, Float64(dt)); push!(obs.sigma, Float64(sigma))
    push!(obs.station, Int32(station)); push!(obs.phase, UInt8(phase)); push!(obs.group, Int32(group))
    return obs
end

function subset_observations(obs::Observations, keep::AbstractVector{Bool})
    length(keep) == length(obs) || error("Observation mask length mismatch")
    return Observations(obs.i[keep], obs.j[keep], obs.dt[keep], obs.sigma[keep],
        obs.station[keep], obs.phase[keep], obs.group[keep])
end

"Build Stage-1 reference-centered theta constraints using serial-ID lookup."
function build_star_observations(groups::Vector{ThetaGroup}, stations::Stations, catalog::Catalog, cfg::AbstractDict)
    event_rows = event_row_map(catalog)
    station_rows = station_row_map(stations)
    minimum_degree = Int(cfgget(cfg, "observations", "minimum_theta_degree"; default=6))
    minimum_sigma = Float64(cfgget(cfg, "prelocation", "min_sigma_s"; default=0.002))
    maximum_sigma = Float64(cfgget(cfg, "observations", "maximum_sigma_s"; default=0.0))
    obs = Observations()
    missing_stations = Set{String}()
    for (group_index, group) in enumerate(groups)
        if !haskey(station_rows, group.station)
            push!(missing_stations, group.station)
            continue
        end
        station = station_rows[group.station]
        for event_id in group.event_id
            entry = group.entries[event_id]
            event_id == entry.ref_id && continue
            haskey(event_rows, event_id) && haskey(event_rows, entry.ref_id) || continue
            haskey(group.entries, entry.ref_id) || continue
            group.has_degree && minimum_degree > 0 && entry.degree < minimum_degree && continue
            reference = group.entries[entry.ref_id]
            sigma = combined_sigma(entry, reference, minimum_sigma, group.has_sigma)
            maximum_sigma > 0.0 && sigma > maximum_sigma && continue
            push_observation!(obs, event_rows[event_id], event_rows[entry.ref_id],
                entry.theta - reference.theta, sigma, station, group.phase, group_index)
        end
    end
    !isempty(missing_stations) && @warn "Theta groups ignored because stations were absent from stations.txt" stations=join(sort!(collect(missing_stations)), ",")
    minimum_per_event = Int(cfgget(cfg, "observations", "minimum_star_observations_per_event"; default=0))
    if minimum_per_event > 0 && !isempty(obs)
        counts = zeros(Int, length(catalog))
        for row in obs.i
            counts[row] += 1
        end
        keep = BitVector([counts[row] >= minimum_per_event for row in obs.i])
        obs = subset_observations(obs, keep)
    end
    @printf("Built %d Stage-1 star observations from %d theta groups\n", length(obs), length(groups))
    return obs
end

function graph_adjacency(graph::EventGraph, n::Int)
    adjacency = [Int32[] for _ in 1:n]
    for edge in eachindex(graph.i)
        i, j = graph.i[edge], graph.j[edge]
        push!(adjacency[i], j); push!(adjacency[j], i)
    end
    for list in adjacency
        sort!(list)
    end
    return adjacency
end

function filter_minimum_pair_support(obs::Observations, minimum::Int)
    minimum <= 1 && return obs
    counts = Dict{UInt64,Int}()
    for row in eachindex(obs.i)
        key = edge_key(obs.i[row], obs.j[row])
        counts[key] = get(counts, key, 0) + 1
    end
    keep = BitVector([get(counts, edge_key(obs.i[row], obs.j[row]), 0) >= minimum for row in eachindex(obs.i)])
    return subset_observations(obs, keep)
end

function peel_low_support_events(obs::Observations, n::Int, minimum::Int; maximum_rounds::Int=10)
    minimum <= 0 && return obs
    current = obs
    for round in 1:maximum_rounds
        counts = zeros(Int, n)
        for row in eachindex(current.i)
            counts[current.i[row]] += 1; counts[current.j[row]] += 1
        end
        keep = BitVector([counts[current.i[row]] >= minimum && counts[current.j[row]] >= minimum for row in eachindex(current.i)])
        all(keep) && break
        before = length(current)
        current = subset_observations(current, keep)
        @printf("  event support peel %d: %d -> %d observations\n", round, before, length(current))
        isempty(current) && break
    end
    return current
end

"Build Stage-2 theta-derived DD observations in O(sum(group nodes × graph degree))."
function build_dd_observations(groups::Vector{ThetaGroup}, stations::Stations, catalog::Catalog,
        graph::EventGraph, cfg::AbstractDict)
    n = length(catalog)
    event_rows = event_row_map(catalog)
    station_rows = station_row_map(stations)
    adjacency = graph_adjacency(graph, n)
    minimum_degree = Int(cfgget(cfg, "observations", "minimum_theta_degree"; default=6))
    minimum_sigma = Float64(cfgget(cfg, "relocation", "min_sigma_s"; default=0.002))
    maximum_sigma = Float64(cfgget(cfg, "observations", "maximum_sigma_s"; default=0.0))
    minimum_component = Int(cfgget(cfg, "observations", "minimum_component_size"; default=0))
    obs = Observations()
    for (group_index, group) in enumerate(groups)
        haskey(station_rows, group.station) || continue
        station = station_rows[group.station]
        for event_id in group.event_id
            haskey(event_rows, event_id) || continue
            i = event_rows[event_id]
            first = group.entries[event_id]
            group.has_degree && minimum_degree > 0 && first.degree < minimum_degree && continue
            for j in adjacency[i]
                j > i || continue
                minimum_component > 0 && (graph.comp_size[i] < minimum_component || graph.comp_size[j] < minimum_component) && continue
                other_id = catalog.event_id[j]
                haskey(group.entries, other_id) || continue
                second = group.entries[other_id]
                first.ref_id == second.ref_id || continue
                group.has_degree && minimum_degree > 0 && second.degree < minimum_degree && continue
                sigma = combined_sigma(first, second, minimum_sigma, group.has_sigma)
                maximum_sigma > 0.0 && sigma > maximum_sigma && continue
                push_observation!(obs, i, j, first.theta - second.theta, sigma, station, group.phase, group_index)
            end
        end
    end
    raw_count = length(obs)
    minimum_pair = Int(cfgget(cfg, "observations", "minimum_observations_per_pair"; default=6))
    obs = filter_minimum_pair_support(obs, minimum_pair)
    minimum_event = Int(cfgget(cfg, "observations", "minimum_observations_per_event"; default=0))
    obs = peel_low_support_events(obs, n, minimum_event)
    @printf("Built %d Stage-2 DD observations (raw=%d) on %d event edges\n", length(obs), raw_count, length(graph))
    return obs
end

function event_activity(obs::Observations, n::Int)
    active = falses(n)
    for row in eachindex(obs.i)
        active[obs.i[row]] = true; active[obs.j[row]] = true
    end
    return active
end

function event_observation_counts(obs::Observations, n::Int)
    counts = zeros(Int, n)
    for row in eachindex(obs.i)
        counts[obs.i[row]] += 1; counts[obs.j[row]] += 1
    end
    return counts
end
