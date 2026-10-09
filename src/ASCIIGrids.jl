"Metadata for a complete regular velocity grid supplied as four-column text."
struct ASCIIVelocityHeader3D
    file::String
    dimensions::NTuple{3,Int}
    origin_m::NTuple{3,Float64}
    spacing_m::NTuple{3,Float64}
    projection::GridProjection3D
    coordinates::String
    reference::NTuple{3,Float64}
end

grid3d_input_bytes(h::NLLHeader3D) = nll_bytes(h)
# Native Float64 velocities, occupancy checks, axes and parsing allowance.
grid3d_input_bytes(h::ASCIIVelocityHeader3D) = 16prod(Int128.(h.dimensions))

function ascii_velocity_row(line, file, line_number, coordinates, reference)
    fields = split(strip(first(split(line, '#'))))
    isempty(fields) && return nothing
    length(fields) == 4 || error("$file:$line_number: expected exactly four columns: lat lon depth_km velocity_km_s, or x_km y_km depth_km velocity_km_s")
    values = tryparse.(Float64, fields)
    all(v -> v !== nothing && isfinite(v), values) || error("$file:$line_number: all four columns must be finite numbers; prefix comments with #")
    a, b, depth, velocity = Float64.(values)
    velocity > 0 || error("$file:$line_number: velocity must be positive, in km/s")
    if coordinates == "geographic"
        abs(a) < 90 && abs(b) <= 180 || error("$file:$line_number: invalid latitude/longitude in degrees")
        x, y = local_xy(a, b, reference...)
    else
        x, y = 1000a, 1000b
    end
    (x, y, 1000depth, 1000velocity)
end

"Scan axes without retaining all rows; reject scattered/incomplete grids before dense allocation."
function read_ascii_velocity_header(path, cfg)
    coordinates = cfg["grid3d"]["coordinate_system"]
    reference = (Float64(cfg["coordinates"]["reference_latitude"]),
        Float64(cfg["coordinates"]["reference_longitude"]), Float64(cfg["coordinates"]["earth_radius_m"]))
    axes = (Set{Float64}(), Set{Float64}(), Set{Float64}())
    count = 0
    # Protect the metadata scan itself if an unstructured cloud has a distinct
    # coordinate in every row. Dense arrays have a separate combined memory guard.
    axis_limit = max(6, floor(Int, cfg["grid3d"]["maximum_memory_gib"] * 2.0^30 / 128))
    for (line_number, line) in enumerate(eachline(path))
        row = ascii_velocity_row(line, path, line_number, coordinates, reference)
        row === nothing && continue
        for a in 1:3; push!(axes[a], row[a]); end
        sum(length, axes) <= axis_limit || error("ASCII coordinate scan exceeds its memory allowance; supply a regular grid or pre-grid scattered samples")
        count += 1
    end
    dims = length.(axes)
    all(>=(2), dims) || error("ASCII velocity input needs at least two distinct nodes along every axis")
    prod(Int128.(dims)) == count || error("ASCII rows must form a complete rectangular grid with one row at every coordinate combination. Scattered or incomplete point clouds must be gridded explicitly before import.")
    prod(Int128.(dims)) <= typemax(Int32) || error("ASCII model exceeds the supported grid size")
    sorted = map(s -> sort!(collect(s)), axes)
    origin = map(first, sorted)
    step = ntuple(a -> (last(sorted[a]) - first(sorted[a])) / (dims[a]-1), 3)
    for a in 1:3, i in eachindex(sorted[a])
        expected = origin[a] + (i-1)*step[a]
        abs(sorted[a][i] - expected) <= max(1e-5, 1e-7*step[a]) ||
            error("ASCII model axis $a is not regularly spaced. Use a regular grid; no scattered-point interpolation or extrapolation is assumed.")
    end
    ASCIIVelocityHeader3D(abspath(path), dims, origin, step,
        GridProjection3D(:local, 0., 0., 0.), coordinates, reference)
end

function same_nll_grid(a::ASCIIVelocityHeader3D, b::ASCIIVelocityHeader3D)
    a.dimensions == b.dimensions && a.origin_m == b.origin_m && a.spacing_m == b.spacing_m &&
        a.coordinates == b.coordinates && a.reference == b.reference
end

function sample_velocity3d(h::ASCIIVelocityHeader3D, g, method, byte_order; sampling="nodes")
    native = FMM3D.Grid3D((collect(range(h.origin_m[a]; step=h.spacing_m[a], length=h.dimensions[a])) for a in 1:3)...)
    values = fill(NaN, h.dimensions)
    for (line_number, line) in enumerate(eachline(h.file))
        row = ascii_velocity_row(line, h.file, line_number, h.coordinates, h.reference)
        row === nothing && continue
        index = ntuple(a -> round(Int, (row[a]-h.origin_m[a])/h.spacing_m[a])+1, 3)
        all(1 <= index[a] <= h.dimensions[a] for a in 1:3) || error("ASCII model changed while being read")
        for a in 1:3
            abs(row[a] - (h.origin_m[a]+(index[a]-1)*h.spacing_m[a])) <= max(1e-5, 1e-7*h.spacing_m[a]) || error("ASCII model changed while being read")
        end
        isnan(values[index...]) || error("Duplicate ASCII model node at line $line_number of $(h.file)")
        values[index...] = row[4]
    end
    all(isfinite, values) || error("ASCII model has missing nodes")
    result = Array{Float64}(undef, size(g))
    for k in eachindex(g.z), j in eachindex(g.y), i in eachindex(g.x)
        xyz = (g.x[i], g.y[j], g.z[k])
        if method == "nearest"
            index = ntuple(a -> begin
                FMM3D._cell_fraction((native.x,native.y,native.z)[a], xyz[a])
                clamp(floor(Int, (xyz[a]-h.origin_m[a])/h.spacing_m[a]+.5)+1, 1, h.dimensions[a])
            end, 3)
            result[i,j,k] = values[index...]
        else
            ids, weights, _ = FMM3D._interpolation(native, xyz)
            result[i,j,k] = 1 / sum(weights[a]/values[ids[a]] for a in 1:8)
        end
    end
    result
end
sample_velocity3d(h::NLLHeader3D, g, method, byte_order; sampling="cell_centers") = sample_nll_velocity(h,g,method,byte_order;sampling=sampling)
