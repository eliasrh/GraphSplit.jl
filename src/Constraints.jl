"Resolve an explicit all-events or serial-ID constraint scope."
function constraint_event_mask(state::State, options::AbstractDict, label::String)
    Bool(get(options, "enabled", false)) || return falses(length(state))
    scope = lowercase(String(get(options, "scope", "all")))
    scope == "all" && return trues(length(state))
    scope == "event_ids" || error("$label.scope must be all or event_ids")
    ids = Int64.(get(options, "event_ids", Int[]))
    rows = Dict(id => row for (row, id) in enumerate(state.event_id))
    missing = [id for id in ids if !haskey(rows, id)]
    isempty(missing) || error("$label contains EventIDs absent from the working catalog: $(join(missing, ','))")
    mask = falses(length(state))
    for id in ids
        mask[rows[id]] = true
    end
    return mask
end

depth_bound_options(cfg::AbstractDict) = cfgget(cfg, "constraints", "depth_bound"; default=Dict{String,Any}())
fixed_depth_options(cfg::AbstractDict) = cfgget(cfg, "constraints", "fixed_depth"; default=Dict{String,Any}())

depth_bound_mask(state::State, cfg::AbstractDict) =
    constraint_event_mask(state, depth_bound_options(cfg), "constraints.depth_bound")
fixed_depth_mask(state::State, cfg::AbstractDict) =
    constraint_event_mask(state, fixed_depth_options(cfg), "constraints.fixed_depth")

"Sign mapping internal z increments to increasing catalog depth."
function internal_depth_direction(cfg::AbstractDict)
    vertical = lowercase(String(cfgget(cfg, "coordinates", "event_vertical"; default="positive_depth")))
    return vertical == "negative_depth" ? -1.0 : 1.0
end

function depth_bound_internal(cfg::AbstractDict)
    value = Float64(get(depth_bound_options(cfg), "minimum_depth_km", 0.0))
    return Float64(catalog_depth_to_internal(value, cfg))
end

function depth_is_feasible(value::Real, boundary::Real, direction::Real; tolerance_m::Float64=1.0e-6)
    return direction * (Float64(value) - Float64(boundary)) >= -tolerance_m
end

function depth_bound_violation_mask(state::State, cfg::AbstractDict; tolerance_m::Float64=1.0e-6)
    scope = depth_bound_mask(state, cfg)
    any(scope) || return scope
    boundary = depth_bound_internal(cfg)
    direction = internal_depth_direction(cfg)
    feasible = BitVector([depth_is_feasible(value, boundary, direction; tolerance_m=tolerance_m)
        for value in state.z])
    return scope .& .!feasible
end

function depth_bound_contact_mask(state::State, cfg::AbstractDict; tolerance_m::Float64=1.0e-5)
    scope = depth_bound_mask(state, cfg)
    any(scope) || return scope
    boundary = depth_bound_internal(cfg)
    direction = internal_depth_direction(cfg)
    return BitVector(scope .& [abs(direction * (value - boundary)) <= tolerance_m for value in state.z])
end

"Events whose z coordinate is frozen by the gauge pin rather than the depth constraint."
function gauge_depth_pin_mask(state::State, cfg::AbstractDict)
    mask = falses(length(state))
    lowercase(String(cfgget(cfg, "gauge", "mode"; default="zero_mean"))) == "pin" || return mask
    fields = lowercase(String(cfgget(cfg, "gauge", "pin_fields"; default="xyz")))
    occursin('z', fields) || return mask
    rows = Dict(id => row for (row, id) in enumerate(state.event_id))
    ids = Int64.(cfgget(cfg, "gauge", "pin_event_ids"; default=Int[]))
    missing = [id for id in ids if !haskey(rows, id)]
    isempty(missing) || error("Pinned event IDs were not found in the catalog: $(join(missing, ','))")
    for id in ids
        mask[rows[id]] = true
    end
    return mask
end

"Return the exact internal-z target for every enabled fixed-depth event."
function fixed_depth_targets(state::State, cfg::AbstractDict)
    options = fixed_depth_options(cfg)
    mask = fixed_depth_mask(state, cfg)
    targets = fill(NaN, length(state))
    any(mask) || return mask, targets
    reference_path = String(get(options, "reference_catalog", ""))
    if isempty(reference_path)
        targets[mask] .= Float64(catalog_depth_to_internal(Float64(get(options, "depth_km", 0.0)), cfg))
        return mask, targets
    end
    reference = read_catalog(reference_path, cfg)
    reference_rows = Dict(id => row for (row, id) in enumerate(reference.event_id))
    for event in findall(mask)
        id = state.event_id[event]
        haskey(reference_rows, id) || error("Fixed-depth ID $id is absent from constraints.fixed_depth.reference_catalog")
        targets[event] = Float64(catalog_depth_to_internal(reference.depth_km[reference_rows[id]], cfg))
    end
    return mask, targets
end

"Apply exact fixed depths once to the seed; the solver subsequently removes their z parameters."
function apply_fixed_depth_constraints!(state::State, cfg::AbstractDict)
    mask, targets = fixed_depth_targets(state, cfg)
    state.z[mask] .= targets[mask]
    return state
end

"Check overlap and feasibility before constructing travel times or solving."
function validate_depth_constraint_state!(state::State, cfg::AbstractDict; allow_free_outside::Bool=false)
    fixed = fixed_depth_mask(state, cfg)
    gauge_pinned = gauge_depth_pin_mask(state, cfg)
    overlap = fixed .& gauge_pinned
    if any(overlap)
        ids = state.event_id[findall(overlap)]
        error("The same z parameter is constrained by both gauge pinning and fixed depth for EventIDs $(join(ids, ',')). Remove z from gauge.pin_fields or remove the overlapping fixed-depth IDs")
    end
    violations = depth_bound_violation_mask(state, cfg)
    frozen_violations = violations .& (fixed .| gauge_pinned)
    if any(frozen_violations)
        ids = state.event_id[findall(frozen_violations)]
        error("Frozen depths violate constraints.depth_bound for EventIDs $(join(ids, ','))")
    end
    if !allow_free_outside && any(violations)
        ids = state.event_id[findall(violations)]
        error("Initial depths violate constraints.depth_bound for EventIDs $(join(ids, ',')). Correct the seeds or enable reflected_prelocation_restart")
    end
    return state
end

function median_value(values::AbstractVector{<:Real})
    isempty(values) && error("Cannot define a station-median mirror plane without stations")
    ordered = sort(Float64.(values))
    middle = length(ordered) ÷ 2
    return isodd(length(ordered)) ? ordered[middle + 1] : 0.5 * (ordered[middle] + ordered[middle + 1])
end

function mirror_plane_internal(stations::Stations, cfg::AbstractDict)
    options = depth_bound_options(cfg)
    plane = lowercase(String(get(options, "mirror_plane", "stations_median")))
    plane == "stations_median" && return median_value(stations.z_m)
    plane == "manual" && return Float64(catalog_depth_to_internal(Float64(get(options, "mirror_depth_km", 0.0)), cfg))
    error("constraints.depth_bound.mirror_plane must be stations_median or manual")
end

"Additional event-depth coordinates that a compatible lookup table must contain."
function constraint_required_depths(state::State, stations::Stations, cfg::AbstractDict)
    values = Float64[]
    options = depth_bound_options(cfg)
    if Bool(get(options, "enabled", false))
        boundary = depth_bound_internal(cfg)
        push!(values, boundary)
        if Bool(get(options, "reflected_prelocation_restart", false))
            margin_m = 1000.0 * Float64(get(options, "pilot_shallow_margin_km", 10.0))
            push!(values, boundary - internal_depth_direction(cfg) * margin_m)
            plane = mirror_plane_internal(stations, cfg)
            push!(values, 2.0 * plane - minimum(state.z))
            push!(values, 2.0 * plane - maximum(state.z))
        end
    end
    return values
end

"Reflect only forbidden free depths into the admissible branch before the bounded restart."
function reflect_depth_violations!(state::State, stations::Stations, cfg::AbstractDict)
    reflected = depth_bound_violation_mask(state, cfg)
    any(reflected) || return reflected
    plane = mirror_plane_internal(stations, cfg)
    state.z[reflected] .= 2.0 .* plane .- state.z[reflected]
    remaining = depth_bound_violation_mask(state, cfg)
    if any(remaining)
        ids = state.event_id[findall(remaining)]
        error("The configured mirror plane did not place EventIDs $(join(ids, ',')) inside the depth bound. Set constraints.depth_bound.mirror_plane=manual and choose a deeper mirror_depth_km")
    end
    return reflected
end
