function strip_comment(line::AbstractString)
    stop = lastindex(line)
    for marker in ('#', '%')
        idx = findfirst(==(marker), line)
        if idx !== nothing
            stop = min(stop, prevind(line, idx))
        end
    end
    return stop < firstindex(line) ? "" : strip(line[firstindex(line):stop])
end

"Read a rectangular whitespace/comma-separated numeric text file without external packages."
function read_numeric_matrix(path::AbstractString)
    isfile(path) || error("Input file not found: $path")
    rows = Vector{Vector{Float64}}()
    ncol = 0
    for (line_number, line) in enumerate(eachline(path))
        clean = replace(strip_comment(line), ',' => ' ')
        isempty(clean) && continue
        fields = split(clean)
        values = Float64[]
        valid = true
        for field in fields
            value = tryparse(Float64, field)
            if value === nothing
                valid = false
                break
            end
            push!(values, value)
        end
        if !valid
            isempty(rows) && continue  # permit a single textual header block
            error("Non-numeric field in $path at line $line_number")
        end
        if ncol == 0
            ncol = length(values)
        elseif length(values) != ncol
            error("Inconsistent column count in $path at line $line_number: expected $ncol, found $(length(values))")
        end
        push!(rows, values)
    end
    isempty(rows) && error("No numeric rows found in $path")
    matrix = Matrix{Float64}(undef, length(rows), ncol)
    for i in eachindex(rows)
        @inbounds matrix[i, :] .= rows[i]
    end
    return matrix
end

function resolved_column(requested::Integer, ncol::Integer, name::String)
    column = requested < 0 ? ncol + requested + 1 : requested
    1 <= column <= ncol || error("$name column $requested resolves outside a $ncol-column file")
    return column
end

"Read the DDSync/GraphSplit catalog; serial event IDs are taken from the configured final column."
function read_catalog(path::AbstractString, cfg::AbstractDict=default_config())
    raw = read_numeric_matrix(path)
    ncol = size(raw, 2)
    ncol >= 4 || error("Catalog must contain at least four columns: $path")
    lat_col = resolved_column(Int(cfgget(cfg, "catalog", "latitude_column"; default=7)), ncol, "latitude")
    lon_col = resolved_column(Int(cfgget(cfg, "catalog", "longitude_column"; default=8)), ncol, "longitude")
    dep_col = resolved_column(Int(cfgget(cfg, "catalog", "depth_column"; default=9)), ncol, "depth")
    id_col = resolved_column(Int(cfgget(cfg, "catalog", "event_id_column"; default=-1)), ncol, "event ID")
    length(unique((lat_col, lon_col, dep_col, id_col))) == 4 || error("Catalog latitude, longitude, depth, and ID columns must be distinct")

    lat = copy(raw[:, lat_col])
    lon = copy(raw[:, lon_col])
    depth = copy(raw[:, dep_col])
    all(isfinite, lat) && all(value -> -90.0 <= value <= 90.0, lat) || error("Invalid catalog latitudes in column $lat_col")
    all(isfinite, lon) && all(value -> -180.0 <= value <= 180.0, lon) || error("Invalid catalog longitudes in column $lon_col")
    all(isfinite, depth) || error("Catalog depth column contains non-finite values")

    ids_float = raw[:, id_col]
    all(isfinite, ids_float) || error("Catalog event IDs must be finite")
    all(abs.(ids_float .- round.(ids_float)) .< 1.0e-8) || error("Catalog event IDs must be integer-valued")
    ids = round.(Int64, ids_float)
    all(value -> value > 0, ids) || error("Catalog serial event IDs must be positive")
    length(unique(ids)) == length(ids) || error("Catalog serial event IDs are not unique")

    return Catalog(String(path), raw, ids, lat, lon, depth, lat_col, lon_col, dep_col, id_col,
        NaN, NaN, Float64(cfgget(cfg, "coordinates", "earth_radius_m"; default=6_371_000.0)))
end

"Read stations as `STA lat lon elevation_m`; `NET STA lat lon elevation_m` is also accepted."
function read_stations(path::AbstractString)
    isfile(path) || error("Stations file not found: $path")
    ids = String[]
    lat = Float64[]
    lon = Float64[]
    elev = Float64[]
    for (line_number, line) in enumerate(eachline(path))
        clean = replace(strip_comment(line), ',' => ' ')
        isempty(clean) && continue
        fields = split(clean)
        if length(fields) == 4
            station_field = 1
            value_fields = 2:4
        elseif length(fields) >= 5
            station_field = 2
            value_fields = 3:5
        else
            error("Stations file line $line_number must contain 4 or 5 fields")
        end
        values = map(x -> tryparse(Float64, x), fields[value_fields])
        any(isnothing, values) && error("Invalid station coordinates at $path:$line_number")
        push!(ids, String(fields[station_field]))
        push!(lat, values[1]::Float64)
        push!(lon, values[2]::Float64)
        push!(elev, values[3]::Float64)
    end
    isempty(ids) && error("No stations found in $path")
    length(unique(ids)) == length(ids) || error("Station identifiers must be unique")
    all(value -> -90.0 <= value <= 90.0, lat) || error("Invalid station latitude")
    all(value -> -180.0 <= value <= 180.0, lon) || error("Invalid station longitude")
    n = length(ids)
    return Stations(ids, lat, lon, elev, zeros(n), zeros(n), zeros(n), NaN, NaN, 6_371_000.0)
end

function theta_name(path::AbstractString)
    filename = basename(path)
    startswith(filename, "theta_") && endswith(lowercase(filename), ".txt") ||
        error("Not a DDSync theta filename: $filename")
    body = filename[7:end-4]
    parts = split(body, '_')
    length(parts) >= 2 || error("Could not parse station and phase from $filename")
    phase_text = uppercase(parts[end])
    phase = phase_text == "P" ? UInt8(1) : phase_text == "S" ? UInt8(2) :
        error("Theta filename phase must be P or S: $filename")
    return join(parts[1:end-1], "_"), phase
end

function read_std_entries(path::AbstractString)
    result = Dict{Int64,Tuple{Float64,Int32}}()
    isfile(path) || return result, false, false
    matrix = read_numeric_matrix(path)
    size(matrix, 2) >= 2 || return result, false, false
    has_degree = size(matrix, 2) >= 4
    for row in axes(matrix, 1)
        idf = matrix[row, 1]
        isfinite(idf) || continue
        id = round(Int64, idf)
        sigma = matrix[row, 2]
        degree = has_degree && isfinite(matrix[row, 4]) ? round(Int32, matrix[row, 4]) : Int32(0)
        result[id] = (sigma, degree)
    end
    return result, true, has_degree
end

"Load sparse, ID-keyed DDSync theta/thetaStd groups. Missing thetaStd is supported."
function load_theta_folder(theta_dir::AbstractString, thetastd_dir::AbstractString="")
    isdir(theta_dir) || error("Theta directory not found: $theta_dir")
    paths = sort(filter(path -> startswith(basename(path), "theta_") && endswith(lowercase(path), ".txt"), readdir(theta_dir; join=true)))
    isempty(paths) && error("No theta_*.txt files found in $theta_dir")
    groups = ThetaGroup[]
    for path in paths
        station, phase = theta_name(path)
        matrix = read_numeric_matrix(path)
        size(matrix, 2) >= 3 || error("Theta file must have at least 3 columns: $path")
        phase_text = phase == 1 ? "P" : "S"
        std_path = isempty(thetastd_dir) ? "" : joinpath(thetastd_dir, "std_theta_$(station)_$(phase_text).txt")
        std, has_sigma, has_degree = read_std_entries(std_path)
        event_ids = Int64[]
        entries = Dict{Int64,ThetaEntry}()
        for row in axes(matrix, 1)
            idf, theta, reff = matrix[row, 1], matrix[row, 2], matrix[row, 3]
            isfinite(idf) && isfinite(theta) && isfinite(reff) || continue
            id = round(Int64, idf)
            ref = round(Int64, reff)
            sigma, degree = get(std, id, (NaN, Int32(0)))
            haskey(entries, id) || push!(event_ids, id)
            entries[id] = ThetaEntry(theta, ref, sigma, degree)
        end
        isempty(entries) && continue
        push!(groups, ThetaGroup(station, phase, basename(path), event_ids, entries, has_sigma, has_degree))
    end
    isempty(groups) && error("Theta files contained no finite event entries")
    return groups
end

function catalog_with_state(catalog::Catalog, state::State, cfg::AbstractDict)
    lat, lon = local_xy_to_ll(state.x, state.y, catalog.ref_lat, catalog.ref_lon, catalog.ref_radius_m)
    vertical = lowercase(String(cfgget(cfg, "coordinates", "event_vertical"; default="positive_depth")))
    z0 = Float64(cfgget(cfg, "coordinates", "event_z0_m"; default=0.0))
    depth = vertical == "positive_depth" ? state.z ./ 1000.0 :
        vertical == "negative_depth" ? -state.z ./ 1000.0 :
        vertical == "positive_depth_plus_z0" ? (state.z .- z0) ./ 1000.0 :
        error("Unknown coordinates.event_vertical: $vertical")
    return lat, lon, depth
end

function write_catalog(path::AbstractString, catalog::Catalog, state::State, cfg::AbstractDict; mask::Union{Nothing,AbstractVector{Bool}}=nothing)
    mkpath(dirname(path))
    lat, lon, depth = catalog_with_state(catalog, state, cfg)
    chosen = mask === nothing ? trues(length(catalog)) : mask
    length(chosen) == length(catalog) || error("Catalog output mask length mismatch")
    open(path, "w") do io
        for row in axes(catalog.raw, 1)
            chosen[row] || continue
            values = copy(@view catalog.raw[row, :])
            values[catalog.lat_col] = lat[row]
            values[catalog.lon_col] = lon[row]
            values[catalog.depth_col] = depth[row]
            for col in eachindex(values)
                if col == catalog.id_col
                    @printf(io, "%d", round(Int64, values[col]))
                elseif col <= 5 && abs(values[col] - round(values[col])) < 1.0e-8
                    @printf(io, "%d", round(Int64, values[col]))
                elseif col == 6
                    @printf(io, "%.3f", values[col])
                elseif col == catalog.lat_col
                    @printf(io, "%.6f", values[col])
                elseif col == catalog.lon_col
                    @printf(io, "%.6f", values[col])
                elseif col == catalog.depth_col
                    @printf(io, "%.4f", values[col])
                elseif col == 10 && col != catalog.id_col
                    @printf(io, "%.3f", values[col])
                else
                    @printf(io, "%.8g", values[col])
                end
                col < length(values) && print(io, ' ')
            end
            println(io)
        end
    end
    return path
end
