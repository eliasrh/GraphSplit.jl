const TT_MAGIC = UInt8[codeunits("GraphSplitTT1")...; 0x00; 0x00; 0x00]
const TT_VERSION = UInt16(1)
const TT_STATUS_OFFSET = 19 # zero-based byte offset after magic, version, and geometry

struct VelocityModel
    file::String
    depth_m::Vector{Float64}
    vp_ms::Vector{Float64}
    vs_ms::Vector{Float64}
    earth_radius_m::Float64
    hash::Vector{UInt8}
end

function travel_time_earth_radius(cfg::AbstractDict, model::Union{Nothing,VelocityModel}, coordinate_radius_m::Float64)
    configured = Float64(cfgget(cfg, "travel_time", "earth_radius_m"; default=0.0))
    configured > 0.0 && isfinite(configured) && return configured
    model !== nothing && isfinite(model.earth_radius_m) && model.earth_radius_m > 0.0 && return model.earth_radius_m
    return coordinate_radius_m
end

function normalize_geometry(value::AbstractString)
    text = lowercase(value)
    text in ("cartesian", "flat") && return :cartesian
    text in ("radial", "spherical") && return :radial
    text == "auto" && return :auto
    error("Unknown travel-time geometry: $value")
end

"Read a 1-D `depth Vp Vs` or `depth radius Vp Vs` model, retaining repeated-depth jumps."
function read_velocity_model(path::AbstractString)
    matrix = read_numeric_matrix(path)
    size(matrix, 2) >= 3 || error("Velocity model requires depth, Vp, and Vs columns")
    radius_layout = false
    if size(matrix, 2) >= 4
        finite2 = filter(isfinite, matrix[:, 2])
        radius_layout = !isempty(finite2) && sort(finite2)[cld(length(finite2), 2)] > 1000.0
    end
    depth = matrix[:, 1]
    radius = radius_layout ? matrix[:, 2] : fill(NaN, size(matrix, 1))
    vp = radius_layout ? matrix[:, 3] : matrix[:, 2]
    vs = radius_layout ? matrix[:, 4] : matrix[:, 3]
    keep = isfinite.(depth) .& isfinite.(vp) .& isfinite.(vs) .& (vp .> 0.0) .& (vs .> 0.0)
    radius_layout && (keep .&= isfinite.(radius) .& (radius .> 0.0))
    depth, radius, vp, vs = depth[keep], radius[keep], vp[keep], vs[keep]
    length(depth) >= 2 || error("Velocity model must contain at least two valid rows with positive P and S velocity")
    order = sortperm(eachindex(depth); by=i -> (depth[i], i))
    depth, radius, vp, vs = depth[order], radius[order], vp[order], vs[order]
    any(diff(depth) .< 0.0) && error("Velocity model depths could not be sorted")
    earth_radius = radius_layout ? sum(depth .+ radius) / length(depth) * 1000.0 : NaN
    return VelocityModel(String(path), depth .* 1000.0, vp .* 1000.0, vs .* 1000.0,
        earth_radius, collect(SHA.sha256(read(path))))
end

"Piecewise-linear layered velocity with the lower-side value at repeated-depth interfaces."
function sample_layered_velocity(depth::Vector{Float64}, velocity::Vector{Float64}, query::Vector{Float64})
    length(depth) == length(velocity) || error("Velocity-model vector length mismatch")
    unique_depth = Float64[]
    above = Float64[]
    below = Float64[]
    i = 1
    while i <= length(depth)
        j = i
        while j < length(depth) && depth[j + 1] == depth[i]
            j += 1
        end
        push!(unique_depth, depth[i])
        push!(above, velocity[i])
        push!(below, velocity[j])
        i = j + 1
    end
    result = Vector{Float64}(undef, length(query))
    for iq in eachindex(query)
        x = query[iq]
        if x < unique_depth[1]
            result[iq] = above[1]
        elseif x >= unique_depth[end]
            result[iq] = below[end]
        else
            k = searchsortedlast(unique_depth, x)
            if x == unique_depth[k]
                result[iq] = below[k]
            else
                fraction = (x - unique_depth[k]) / (unique_depth[k + 1] - unique_depth[k])
                result[iq] = (1.0 - fraction) * below[k] + fraction * above[k + 1]
            end
        end
    end
    return result
end

function great_circle_delta(lat1, lon1, lat2, lon2)
    phi1, phi2 = Float64(lat1) * DEG2RAD, Float64(lat2) * DEG2RAD
    dlambda = atan(sin((Float64(lon1) - Float64(lon2)) * DEG2RAD), cos((Float64(lon1) - Float64(lon2)) * DEG2RAD))
    u = clamp(sin(phi1) * sin(phi2) + cos(phi1) * cos(phi2) * cos(dlambda), -1.0, 1.0)
    return atan(sqrt(max(0.0, 1.0 - u * u)), u)
end

function required_range(stations::Stations, catalog::Catalog, state::State,
        geometry::Symbol, earth_radius_m::Float64)
    maximum = 0.0
    if geometry == :cartesian
        cx, cy = local_xy(catalog.lat, catalog.lon, stations.ref_lat, stations.ref_lon, stations.ref_radius_m)
        for s in eachindex(stations.id), e in eachindex(cx)
            maximum = max(maximum, hypot(cx[e] - stations.x_m[s], cy[e] - stations.y_m[s]))
        end
        for s in eachindex(stations.id), e in eachindex(state.x)
            maximum = max(maximum, hypot(state.x[e] - stations.x_m[s], state.y[e] - stations.y_m[s]))
        end
    else
        for s in eachindex(stations.id), e in eachindex(catalog.lat)
            maximum = max(maximum, earth_radius_m * great_circle_delta(catalog.lat[e], catalog.lon[e], stations.lat[s], stations.lon[s]))
        end
        state_lat, state_lon = local_xy_to_ll(state.x, state.y,
            state.ref_lat, state.ref_lon, state.ref_radius_m)
        for s in eachindex(stations.id), e in eachindex(state_lat)
            maximum = max(maximum, earth_radius_m * great_circle_delta(state_lat[e], state_lon[e], stations.lat[s], stations.lon[s]))
        end
    end
    return maximum
end

function auto_number(value, fallback::Float64; positive::Bool=false)
    if value isa Real
        number = Float64(value)
        if isfinite(number) && (!positive || number > 0.0)
            return number
        end
    end
    return fallback
end

function lookup_grids(cfg::AbstractDict, stations::Stations, catalog::Catalog, state::State,
        model::VelocityModel, geometry::Symbol, earth_radius_m::Float64)
    dq_surface = Float64(cfgget(cfg, "lookup", "horizontal_step_m"))
    dz = Float64(cfgget(cfg, "lookup", "depth_step_m"))
    dzs_requested = Float64(cfgget(cfg, "lookup", "station_depth_step_m"))
    dq_surface > 0.0 && dz > 0.0 && dzs_requested > 0.0 || error("Lookup grid steps must be positive")
    required = required_range(stations, catalog, state, geometry, earth_radius_m)
    margin = 1000.0 * Float64(cfgget(cfg, "lookup", "distance_margin_km"; default=1.0))
    configured_max = cfgget(cfg, "lookup", "maximum_distance_km"; default=0.0)
    rmax = auto_number(configured_max, required + max(margin, 2dq_surface); positive=true) *
        (configured_max isa Real && Float64(configured_max) > 0.0 ? 1000.0 : 1.0)
    rmax + 1.0e-6 >= required || error(@sprintf("Configured lookup maximum distance %.3f km is smaller than the required %.3f km", rmax / 1000, required / 1000))

    required_constraint_depths = constraint_required_depths(state, stations, cfg)
    auto_zmin = min(minimum(model.depth_m), minimum(stations.z_m), minimum(state.z),
        isempty(required_constraint_depths) ? Inf : minimum(required_constraint_depths))
    configured_zmin = cfgget(cfg, "lookup", "minimum_depth_km"; default=0.0)
    zmin = configured_zmin isa Real && isfinite(Float64(configured_zmin)) && Float64(configured_zmin) != 0.0 ?
        1000.0 * Float64(configured_zmin) : auto_zmin
    depth_margin = 1000.0 * Float64(cfgget(cfg, "lookup", "depth_margin_km"; default=50.0))
    auto_zmax = max(maximum(state.z) + depth_margin, 0.35rmax,
        min(maximum(model.depth_m), 100_000.0),
        isempty(required_constraint_depths) ? -Inf : maximum(required_constraint_depths))
    configured_zmax = cfgget(cfg, "lookup", "maximum_depth_km"; default=0.0)
    zmax = configured_zmax isa Real && isfinite(Float64(configured_zmax)) && Float64(configured_zmax) > 0.0 ?
        1000.0 * Float64(configured_zmax) : auto_zmax
    zmin <= minimum(stations.z_m) && zmax >= maximum(stations.z_m) || error("Lookup depth range does not contain every station")
    zmin <= minimum(state.z) && zmax >= maximum(state.z) || error("Lookup depth range does not contain every input event")
    all(zmin <= required_depth <= zmax for required_depth in required_constraint_depths) ||
        error("Lookup depth range does not contain the configured physical/reflected depth-constraint coverage")

    izmin, izmax = floor(Int, zmin / dz), ceil(Int, zmax / dz)
    z = collect((izmin:izmax) .* dz)
    length(z) >= 2 || error("Lookup depth grid requires at least two nodes")
    q = geometry == :cartesian ? collect(0.0:dq_surface:(ceil(rmax / dq_surface) * dq_surface)) :
        collect(0.0:(dq_surface / earth_radius_m):(ceil(rmax / dq_surface) * dq_surface / earth_radius_m))
    length(q) >= 2 || error("Lookup range grid requires at least two nodes")

    source_stride = max(1, floor(Int, dzs_requested / dz))
    first_source = max(1, searchsortedlast(z, minimum(stations.z_m)))
    last_source = min(length(z), searchsortedfirst(z, maximum(stations.z_m)))
    source_indices = collect(first_source:source_stride:last_source)
    isempty(source_indices) && push!(source_indices, first_source)
    source_indices[end] != last_source && push!(source_indices, last_source)
    if length(source_indices) == 1
        push!(source_indices, source_indices[1] < length(z) ? source_indices[1] + 1 : source_indices[1] - 1)
        sort!(source_indices)
    end
    zs = z[source_indices]
    return q, z, zs, required
end

function local_eikonal_update(a::Float64, b::Float64, hz::Float64, hq::Float64, slowness::Float64)
    !isfinite(a) && !isfinite(b) && return Inf
    !isfinite(a) && return b + slowness * hq
    !isfinite(b) && return a + slowness * hz
    one_sided = min(a + slowness * hz, b + slowness * hq)
    az, aq = inv(hz * hz), inv(hq * hq)
    aa = az + aq
    bb = -2.0 * (a * az + b * aq)
    cc = a * a * az + b * b * aq - slowness * slowness
    discriminant = bb * bb - 4.0 * aa * cc
    discriminant < 0.0 && return one_sided
    root = (-bb + sqrt(max(0.0, discriminant))) / (2.0 * aa)
    return root >= max(a, b) ? root : one_sided
end

"First-order Godunov fast sweeping in flat cylindrical or axisymmetric spherical coordinates."
function fast_sweep_field(velocity::Vector{Float64}, z::Vector{Float64}, q::Vector{Float64},
        source_z_index::Int, geometry::Symbol, earth_radius_m::Float64, max_sweeps::Int, tolerance::Float64)
    nz, nq = length(z), length(q)
    dz = z[2] - z[1]
    dq = q[2] - q[1]
    travel = fill(Inf, nz, nq)
    travel[source_z_index, 1] = 0.0
    slowness = 1.0 ./ velocity
    radius = earth_radius_m .- z
    offsets = ((-1, 0), (1, 0), (0, 1), (-1, 1), (1, 1))
    for (oz, oq) in offsets
        iz, iq = source_z_index + oz, 1 + oq
        1 <= iz <= nz && 1 <= iq <= nq || continue
        vertical = abs(z[iz] - z[source_z_index])
        horizontal = geometry == :cartesian ? abs(q[iq] - q[1]) :
            0.5 * (radius[iz] + radius[source_z_index]) * abs(q[iq] - q[1])
        travel[iz, iq] = hypot(vertical, horizontal) * 0.5 * (slowness[iz] + slowness[source_z_index])
    end

    for sweep in 1:max_sweeps
        maximum_change = 0.0
        for zforward in (true, false), qforward in (true, false)
            zrange = zforward ? (1:nz) : (nz:-1:1)
            qrange = qforward ? (1:nq) : (nq:-1:1)
            for iz in zrange
                hq = geometry == :cartesian ? dq : radius[iz] * dq
                radius[iz] > 0.0 || error("Radial lookup reached or crossed the Earth center")
                for iq in qrange
                    iz == source_z_index && iq == 1 && continue
                    a = iz == 1 ? travel[2, iq] : iz == nz ? travel[nz - 1, iq] : min(travel[iz - 1, iq], travel[iz + 1, iq])
                    b = iq == 1 ? travel[iz, 2] : iq == nq ? travel[iz, nq - 1] : min(travel[iz, iq - 1], travel[iz, iq + 1])
                    candidate = local_eikonal_update(a, b, dz, hq, slowness[iz])
                    previous = travel[iz, iq]
                    if candidate < previous
                        travel[iz, iq] = candidate
                        maximum_change = isfinite(previous) ? max(maximum_change, previous - candidate) : Inf
                    end
                end
            end
        end
        all(isfinite, travel) && isfinite(maximum_change) && maximum_change < tolerance && return travel, sweep, true
    end
    return travel, max_sweeps, false
end

function write_table_header(io::IO, geometry::Symbol, earth_radius_m::Float64, coordinate_radius_m::Float64,
        model_hash::Vector{UInt8}, q::Vector{Float64}, z::Vector{Float64}, zs::Vector{Float64}; complete::Bool=false)
    length(TT_MAGIC) == 16 || error("Internal table magic length error")
    length(model_hash) == 32 || error("Model hash must contain 32 bytes")
    write(io, TT_MAGIC)
    write(io, TT_VERSION)
    write(io, geometry == :cartesian ? UInt8(1) : UInt8(2))
    write(io, complete ? UInt8(1) : UInt8(0)) # 0=incomplete, 1=complete Float32 table
    write(io, earth_radius_m)
    write(io, coordinate_radius_m)
    write(io, Int64(length(q)), Int64(length(z)), Int64(length(zs)))
    write(io, model_hash)
    write(io, q, z, zs)
    return position(io)
end

function read_table_header(io::IO)
    seekstart(io)
    magic = Vector{UInt8}(undef, 16)
    read!(io, magic)
    magic == TT_MAGIC || error("Not a GraphSplit native travel-time table")
    version = read(io, UInt16)
    version == TT_VERSION || error("Unsupported travel-time table version: $version")
    geometry_code = read(io, UInt8)
    scalar_code = read(io, UInt8)
    scalar_code == 0 && error("Travel-time table is incomplete (a previous build did not finish)")
    scalar_code == 1 || error("Unsupported travel-time table scalar type")
    geometry = geometry_code == 1 ? :cartesian : geometry_code == 2 ? :radial : error("Invalid geometry code in travel-time table")
    earth_radius = read(io, Float64)
    coordinate_radius = read(io, Float64)
    dims = (Int(read(io, Int64)), Int(read(io, Int64)), Int(read(io, Int64)))
    all(dimension -> dimension > 1, dims) || error("Travel-time table axes must each have at least two nodes")
    model_hash = Vector{UInt8}(undef, 32)
    read!(io, model_hash)
    q, z, zs = Vector{Float64}(undef, dims[1]), Vector{Float64}(undef, dims[2]), Vector{Float64}(undef, dims[3])
    read!(io, q); read!(io, z); read!(io, zs)
    return geometry, earth_radius, coordinate_radius, dims, model_hash, q, z, zs, position(io)
end

function mark_table_complete!(table::TravelTimeTable)
    seek(table.io, TT_STATUS_OFFSET)
    write(table.io, UInt8(1))
    flush(table.io)
    return table
end

function create_empty_table(path::AbstractString, geometry::Symbol, earth_radius_m::Float64,
        coordinate_radius_m::Float64, model_hash::Vector{UInt8}, q, z, zs, ref_lat, ref_lon, clamp)
    mkpath(dirname(path))
    io = open(path, "w+")
    p_offset = write_table_header(io, geometry, earth_radius_m, coordinate_radius_m, model_hash, q, z, zs)
    dims = (length(q), length(z), length(zs))
    values_per_phase = prod(dims)
    s_offset = p_offset + sizeof(Float32) * values_per_phase
    total_bytes = s_offset + sizeof(Float32) * values_per_phase
    seek(io, total_bytes - 1)
    write(io, UInt8(0)); flush(io)
    p = Mmap.mmap(io, Array{Float32,3}, dims, p_offset)
    s = Mmap.mmap(io, Array{Float32,3}, dims, s_offset)
    return TravelTimeTable(String(path), geometry, earth_radius_m, coordinate_radius_m, copy(model_hash),
        collect(q), collect(z), collect(zs), p, s, ref_lat, ref_lon, clamp, io)
end

function load_travel_time_table(path::AbstractString, ref_lat::Float64, ref_lon::Float64, clamp::Bool)
    io = open(path, "r")
    geometry, earth_radius, coordinate_radius, dims, model_hash, q, z, zs, p_offset = read_table_header(io)
    values_per_phase = prod(dims)
    s_offset = p_offset + sizeof(Float32) * values_per_phase
    expected = s_offset + sizeof(Float32) * values_per_phase
    filesize(path) >= expected || (close(io); error("Travel-time table is truncated: $path"))
    p = Mmap.mmap(io, Array{Float32,3}, dims, p_offset)
    s = Mmap.mmap(io, Array{Float32,3}, dims, s_offset)
    return TravelTimeTable(String(path), geometry, earth_radius, coordinate_radius, model_hash, q, z, zs,
        p, s, ref_lat, ref_lon, clamp, io)
end

function compare_table_coverage(table::TravelTimeTable, cfg::AbstractDict, stations::Stations,
        catalog::Catalog, state::State, model::VelocityModel)
    requested = normalize_geometry(String(cfgget(cfg, "travel_time", "geometry")))
    requested != :auto && requested != table.geometry && return false, "geometry mismatch"
    model.hash == table.model_hash || return false, "velocity model changed"
    intended_radius = travel_time_earth_radius(cfg, model, stations.ref_radius_m)
    abs(intended_radius - table.earth_radius_m) <= 1.0e-10 * intended_radius || return false, "travel-time Earth radius changed"
    abs(stations.ref_radius_m - table.coordinate_radius_m) <= 1.0e-10 * stations.ref_radius_m || return false, "coordinate Earth radius changed"
    required = required_range(stations, catalog, state, table.geometry, table.earth_radius_m)
    qmax_m = table.geometry == :cartesian ? table.q[end] : table.q[end] * table.earth_radius_m
    required <= qmax_m + 1.0e-6 || return false, "distance coverage is too small"
    minimum(stations.z_m) >= table.z[1] - 1.0e-6 && maximum(stations.z_m) <= table.z[end] + 1.0e-6 ||
        return false, "station depth coverage is too small"
    minimum(state.z) >= table.z[1] - 1.0e-6 && maximum(state.z) <= table.z[end] + 1.0e-6 ||
        return false, "event depth coverage is too small"
    for required_depth in constraint_required_depths(state, stations, cfg)
        table.z[1] - 1.0e-6 <= required_depth <= table.z[end] + 1.0e-6 ||
            return false, "physical/reflected depth-constraint coverage is too small"
    end
    minimum(stations.z_m) >= table.zs[1] - 1.0e-6 && maximum(stations.z_m) <= table.zs[end] + 1.0e-6 ||
        return false, "station-depth interpolation coverage is too small"
    requested_step = Float64(cfgget(cfg, "lookup", "horizontal_step_m"))
    actual_step = table.geometry == :cartesian ? table.q[2] - table.q[1] : (table.q[2] - table.q[1]) * table.earth_radius_m
    abs(requested_step - actual_step) <= 1.0e-8 * max(1.0, requested_step) || return false, "horizontal grid spacing changed"
    requested_dz = Float64(cfgget(cfg, "lookup", "depth_step_m"))
    abs(requested_dz - (table.z[2] - table.z[1])) <= 1.0e-8 * max(1.0, requested_dz) || return false, "depth grid spacing changed"
    requested_dzs = Float64(cfgget(cfg, "lookup", "station_depth_step_m"))
    maximum(diff(table.zs)) <= requested_dzs + 1.0e-8 * max(1.0, requested_dzs) ||
        return false, "station-depth grid is coarser than requested"
    return true, "compatible"
end

function fill_phase!(cube::Array{Float32,3}, velocity::Vector{Float64}, table::TravelTimeTable,
        max_sweeps::Int, tolerance::Float64, threaded::Bool, phase_name::String)
    function build_slice(ks)
        source_index = searchsortedfirst(table.z, table.zs[ks])
        field, iterations, converged = fast_sweep_field(velocity, table.z, table.q, source_index,
            table.geometry, table.earth_radius_m, max_sweeps, tolerance)
        if !converged
            @warn @sprintf("%s travel-time slice at station depth %.1f m did not converge in %d sweeps",
                phase_name, table.zs[ks], iterations)
        end
        @inbounds for iz in eachindex(table.z), iq in eachindex(table.q)
            cube[iq, iz, ks] = Float32(field[iz, iq])
        end
    end
    if threaded && Threads.nthreads() > 1
        Threads.@threads for ks in eachindex(table.zs)
            build_slice(ks)
        end
    else
        for ks in eachindex(table.zs)
            build_slice(ks)
            @printf("  %s station-depth slice %d/%d\n", phase_name, ks, length(table.zs))
        end
    end
    return cube
end

function build_lookup_table(cfg::AbstractDict, stations::Stations, catalog::Catalog, state::State,
        model::VelocityModel, geometry::Symbol)
    earth_radius = travel_time_earth_radius(cfg, model, stations.ref_radius_m)
    q, z, zs, required = lookup_grids(cfg, stations, catalog, state, model, geometry, earth_radius)
    table_path = String(cfgget(cfg, "travel_time", "table_file"))
    clamp = Bool(cfgget(cfg, "travel_time", "clamp_to_grid"; default=true))
    dims = (length(q), length(z), length(zs))
    gib = 2.0 * sizeof(Float32) * prod(dims) / 2.0^30
    @printf("Building %s travel-time table %s\n", geometry == :cartesian ? "Cartesian" : "radial", table_path)
    @printf("  grid: %d x %d x %d; required range %.3f km; saved values %.3f GiB\n",
        dims..., required / 1000.0, gib)
    maximum_gib = Float64(cfgget(cfg, "lookup", "maximum_table_gib"; default=8.0))
    allow_large = Bool(cfgget(cfg, "lookup", "allow_large_table"; default=false))
    gib <= maximum_gib || allow_large || error(@sprintf("Requested table is %.2f GiB, above lookup.maximum_table_gib=%.2f", gib, maximum_gib))
    table = create_empty_table(table_path, geometry, earth_radius, stations.ref_radius_m, model.hash,
        q, z, zs, stations.ref_lat, stations.ref_lon, clamp)
    vp = sample_layered_velocity(model.depth_m, model.vp_ms, z)
    vs = sample_layered_velocity(model.depth_m, model.vs_ms, z)
    max_sweeps = Int(cfgget(cfg, "lookup", "maximum_sweeps"; default=32))
    tolerance = Float64(cfgget(cfg, "lookup", "sweep_tolerance_s"; default=1.0e-7))
    threaded = Bool(cfgget(cfg, "lookup", "threaded"; default=false))
    fill_phase!(table.p, vp, table, max_sweeps, tolerance, threaded, "P")
    fill_phase!(table.s, vs, table, max_sweeps, tolerance, threaded, "S")
    Mmap.sync!(table.p); Mmap.sync!(table.s)
    mark_table_complete!(table)
    return table
end

function prepare_travel_time(cfg::AbstractDict, stations::Stations, catalog::Catalog, state::State;
        force_build::Bool=false)
    kind = lowercase(String(cfgget(cfg, "travel_time", "type"; default="lookup")))
    kind == "3d" && return prepare_travel_time_3d(cfg, stations, catalog, state; force_build=force_build)
    requested = normalize_geometry(String(cfgget(cfg, "travel_time", "geometry"; default="auto")))
    if kind in ("constant", "constant_velocity", "constvel")
        geometry = requested == :auto ? normalize_geometry(String(cfgget(cfg, "travel_time", "build_geometry"))) : requested
        earth_radius = travel_time_earth_radius(cfg, nothing, stations.ref_radius_m)
        return ConstantTravelTime(geometry, Float64(cfgget(cfg, "travel_time", "vp_ms")),
            Float64(cfgget(cfg, "travel_time", "vs_ms")), earth_radius,
            stations.ref_lat, stations.ref_lon, stations.ref_radius_m)
    end
    kind == "lookup" || error("travel_time.type must be lookup or constant")
    model_path = String(cfgget(cfg, "travel_time", "velocity_model_file"))
    model = read_velocity_model(model_path)
    table_path = String(cfgget(cfg, "travel_time", "table_file"))
    if isfile(table_path) && !force_build
        table = nothing
        reason = "unknown incompatibility"
        try
            table = load_travel_time_table(table_path, stations.ref_lat, stations.ref_lon,
                Bool(cfgget(cfg, "travel_time", "clamp_to_grid"; default=true)))
            compatible, reason = compare_table_coverage(table, cfg, stations, catalog, state, model)
            compatible && return table
            close(table.io)
            table = nothing
        catch exception
            table !== nothing && isopen(table.io) && close(table.io)
            reason = sprint(showerror, exception)
        end
        Bool(cfgget(cfg, "travel_time", "auto_rebuild"; default=true)) || error("Existing travel-time table is incompatible ($reason), and auto_rebuild=false")
        @info "Rebuilding incompatible travel-time table" reason=reason
    elseif !isfile(table_path) && !force_build && !Bool(cfgget(cfg, "travel_time", "auto_build"; default=true))
        error("Travel-time table not found and travel_time.auto_build=false: $table_path")
    end
    geometry = requested == :auto ? normalize_geometry(String(cfgget(cfg, "travel_time", "build_geometry"))) : requested
    return build_lookup_table(cfg, stations, catalog, state, model, geometry)
end

function grid_cell(grid::Vector{Float64}, value::Float64, clamp_to_grid::Bool)
    outside = value < grid[1] || value > grid[end]
    outside && !clamp_to_grid && return 0, 0.0, true
    v = clamp(value, grid[1], grid[end])
    index = searchsortedlast(grid, v)
    index >= length(grid) && (index = length(grid) - 1)
    fraction = (v - grid[index]) / (grid[index + 1] - grid[index])
    return index, fraction, outside
end

function trilinear_value_gradient(values::Array{Float32,3}, qgrid, zgrid, zsgrid,
        q::Float64, z::Float64, zs::Float64, clamp_to_grid::Bool)
    iq, tq, oq = grid_cell(qgrid, q, clamp_to_grid)
    iz, tz, oz = grid_cell(zgrid, z, clamp_to_grid)
    is, ts, os = grid_cell(zsgrid, zs, clamp_to_grid)
    (iq == 0 || iz == 0 || is == 0) && return NaN, 0.0, 0.0
    v000 = Float64(values[iq, iz, is]);       v100 = Float64(values[iq + 1, iz, is])
    v010 = Float64(values[iq, iz + 1, is]);   v110 = Float64(values[iq + 1, iz + 1, is])
    v001 = Float64(values[iq, iz, is + 1]);   v101 = Float64(values[iq + 1, iz, is + 1])
    v011 = Float64(values[iq, iz + 1, is + 1]); v111 = Float64(values[iq + 1, iz + 1, is + 1])
    q00 = muladd(tq, v100 - v000, v000); q10 = muladd(tq, v110 - v010, v010)
    q01 = muladd(tq, v101 - v001, v001); q11 = muladd(tq, v111 - v011, v011)
    z0 = muladd(tz, q10 - q00, q00); z1 = muladd(tz, q11 - q01, q01)
    value = muladd(ts, z1 - z0, z0)
    dq0 = ((1 - tz) * (v100 - v000) + tz * (v110 - v010)) / (qgrid[iq + 1] - qgrid[iq])
    dq1 = ((1 - tz) * (v101 - v001) + tz * (v111 - v011)) / (qgrid[iq + 1] - qgrid[iq])
    dz0 = ((1 - tq) * (v010 - v000) + tq * (v110 - v100)) / (zgrid[iz + 1] - zgrid[iz])
    dz1 = ((1 - tq) * (v011 - v001) + tq * (v111 - v101)) / (zgrid[iz + 1] - zgrid[iz])
    dtdq = oq ? 0.0 : muladd(ts, dq1 - dq0, dq0)
    dtdz = oz ? 0.0 : muladd(ts, dz1 - dz0, dz0)
    return value, dtdq, dtdz
end

function travel_time_gradient(model::ConstantTravelTime, phase::UInt8, x, y, z, sx, sy, sz)
    velocity = phase == 1 ? model.vp_ms : model.vs_ms
    if model.geometry == :cartesian
        dx, dy, dz = x - sx, y - sy, z - sz
        distance = max(eps(Float64), sqrt(dx * dx + dy * dy + dz * dz))
        return distance / velocity, dx / (velocity * distance), dy / (velocity * distance), dz / (velocity * distance)
    end
    delta, ddelta_dx, ddelta_dy = spherical_delta_gradient(x, y, sx, sy,
        model.ref_lat, model.ref_lon, model.coordinate_radius_m)
    re, rs = model.earth_radius_m - z, model.earth_radius_m - sz
    distance = max(eps(Float64), sqrt(max(0.0, re * re + rs * rs - 2re * rs * cos(delta))))
    dtdq = re * rs * sin(delta) / (velocity * distance)
    dtdz = -(re - rs * cos(delta)) / (velocity * distance)
    return distance / velocity, dtdq * ddelta_dx, dtdq * ddelta_dy, dtdz
end

function travel_time_gradient(model::TravelTimeTable, phase::UInt8, x, y, z, sx, sy, sz)
    if model.geometry == :cartesian
        dx, dy = x - sx, y - sy
        q = hypot(dx, dy)
        if q > 1.0e-12
            dqdx, dqdy = dx / q, dy / q
        else
            dqdx, dqdy = 0.0, 0.0
        end
    else
        q, dqdx, dqdy = spherical_delta_gradient(x, y, sx, sy,
            model.ref_lat, model.ref_lon, model.coordinate_radius_m)
    end
    values = phase == 1 ? model.p : phase == 2 ? model.s : error("Phase code must be 1 (P) or 2 (S)")
    value, dtdq, dtdz = trilinear_value_gradient(values, model.q, model.z, model.zs,
        Float64(q), Float64(z), Float64(sz), model.clamp)
    return value, dtdq * dqdx, dtdq * dqdy, dtdz
end
