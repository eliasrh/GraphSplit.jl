#!/usr/bin/env julia

# Standalone, standard-library-only reproduction of the location metrics in
# Yu, Ellsworth & Beroza (2025), following their public error_analysis.py.

using Printf

const WGS84_A_M = 6_378_137.0
const WGS84_F = 1.0 / 298.257_223_563
const WGS84_B_M = WGS84_A_M * (1.0 - WGS84_F)
const WGS84_E2 = WGS84_F * (2.0 - WGS84_F)

function numeric_matrix(path)
    rows = Vector{Vector{Float64}}()
    for line in eachline(path)
        clean = strip(first(split(line, '#'; limit=2)))
        isempty(clean) && continue
        fields = split(replace(clean, ',' => ' '))
        values = try
            parse.(Float64, fields)
        catch
            continue
        end
        push!(rows, values)
    end
    isempty(rows) && error("No numeric rows in $path")
    width = length(rows[1])
    all(row -> length(row) == width, rows) || error("Ragged numeric file: $path")
    matrix = Matrix{Float64}(undef, length(rows), width)
    for row in eachindex(rows)
        matrix[row, :] .= rows[row]
    end
    return matrix
end

average(values) = isempty(values) ? NaN : sum(values) / length(values)
function middle(values)
    isempty(values) && return NaN
    ordered = sort(values); n = length(ordered)
    return isodd(n) ? ordered[(n + 1) >>> 1] : 0.5 * (ordered[n >>> 1] + ordered[(n >>> 1) + 1])
end

"WGS84 inverse-geodesic distance in kilometres (Vincenty, local fallback for a rare non-convergence)."
function horizontal_distance_km(lat1, lon1, lat2, lon2)
    phi1, phi2 = deg2rad(lat1), deg2rad(lat2)
    longitude_difference = atan(sin(deg2rad(lon2 - lon1)), cos(deg2rad(lon2 - lon1)))
    u1, u2 = atan((1.0 - WGS84_F) * tan(phi1)), atan((1.0 - WGS84_F) * tan(phi2))
    sinu1, cosu1 = sincos(u1); sinu2, cosu2 = sincos(u2)
    lambda = longitude_difference
    sigma = 0.0; sinsigma = 0.0; cossigma = 1.0
    sinalpha = 0.0; cos2alpha = 1.0; cos2sigmam = 0.0
    converged = false
    for _ in 1:100
        sinlambda, coslambda = sincos(lambda)
        sinsigma = hypot(cosu2 * sinlambda, cosu1 * sinu2 - sinu1 * cosu2 * coslambda)
        sinsigma == 0.0 && return 0.0
        cossigma = sinu1 * sinu2 + cosu1 * cosu2 * coslambda
        sigma = atan(sinsigma, cossigma)
        sinalpha = cosu1 * cosu2 * sinlambda / sinsigma
        cos2alpha = 1.0 - sinalpha * sinalpha
        cos2sigmam = cos2alpha > eps(Float64) ? cossigma - 2.0 * sinu1 * sinu2 / cos2alpha : 0.0
        coefficient = WGS84_F / 16.0 * cos2alpha * (4.0 + WGS84_F * (4.0 - 3.0 * cos2alpha))
        previous = lambda
        lambda = longitude_difference + (1.0 - coefficient) * WGS84_F * sinalpha *
            (sigma + coefficient * sinsigma * (cos2sigmam + coefficient * cossigma *
            (-1.0 + 2.0 * cos2sigmam^2)))
        if abs(lambda - previous) <= 1.0e-12
            converged = true
            break
        end
    end
    if !converged
        central = 2.0 * asin(sqrt(sin(0.5 * (phi2 - phi1))^2 + cos(phi1) * cos(phi2) *
            sin(0.5 * longitude_difference)^2))
        return 6371.0088 * central
    end
    u2coefficient = cos2alpha * (WGS84_A_M^2 - WGS84_B_M^2) / WGS84_B_M^2
    aa = 1.0 + u2coefficient / 16384.0 * (4096.0 + u2coefficient *
        (-768.0 + u2coefficient * (320.0 - 175.0 * u2coefficient)))
    bb = u2coefficient / 1024.0 * (256.0 + u2coefficient *
        (-128.0 + u2coefficient * (74.0 - 47.0 * u2coefficient)))
    delta_sigma = bb * sinsigma * (cos2sigmam + 0.25 * bb *
        (cossigma * (-1.0 + 2.0 * cos2sigmam^2) - bb / 6.0 * cos2sigmam *
        (-3.0 + 4.0 * sinsigma^2) * (-3.0 + 4.0 * cos2sigmam^2)))
    return WGS84_B_M * aa * (sigma - delta_sigma) / 1000.0
end

function geodetic_ecef(lat_deg, lon_deg, height_m)
    phi, lambda = deg2rad(lat_deg), deg2rad(lon_deg)
    sinphi, cosphi = sincos(phi); sinlambda, coslambda = sincos(lambda)
    prime = WGS84_A_M / sqrt(1.0 - WGS84_E2 * sinphi^2)
    return ((prime + height_m) * cosphi * coslambda,
        (prime + height_m) * cosphi * sinlambda,
        (prime * (1.0 - WGS84_E2) + height_m) * sinphi)
end

"Match the paper code's pymap3d ENU conversion about (35.4, -117.956, 0)."
function benchmark_enu_km(lat, lon, depth_km)
    lat0, lon0 = 35.4, -117.956
    x0, y0, z0 = geodetic_ecef(lat0, lon0, 0.0)
    phi0, lambda0 = deg2rad(lat0), deg2rad(lon0)
    sinphi0, cosphi0 = sincos(phi0); sinlambda0, coslambda0 = sincos(lambda0)
    east, north, down = zeros(length(lat)), zeros(length(lat)), zeros(length(lat))
    for i in eachindex(lat)
        x, y, z = geodetic_ecef(lat[i], lon[i], -1000.0 * depth_km[i])
        dx, dy, dz = x - x0, y - y0, z - z0
        east[i] = (-sinlambda0 * dx + coslambda0 * dy) / 1000.0
        north[i] = (-sinphi0 * coslambda0 * dx - sinphi0 * sinlambda0 * dy + cosphi0 * dz) / 1000.0
        down[i] = (-(cosphi0 * coslambda0 * dx + cosphi0 * sinlambda0 * dy + sinphi0 * dz)) / 1000.0
    end
    return east, north, down
end

"Point Cloud Utils convention: half the sum of the two mean nearest-neighbor Euclidean distances."
function chamfer_distance(x1, y1, z1, x2, y2, z2)
    function directed(ax, ay, az, bx, by, bz)
        total = 0.0
        for i in eachindex(ax)
            best2 = Inf
            for j in eachindex(bx)
                distance2 = (ax[i] - bx[j])^2 + (ay[i] - by[j])^2 + (az[i] - bz[j])^2
                best2 = min(best2, distance2)
            end
            total += sqrt(best2)
        end
        return total / length(ax)
    end
    return 0.5 * (directed(x1, y1, z1, x2, y2, z2) + directed(x2, y2, z2, x1, y1, z1))
end

function main(args=ARGS)
    length(args) >= 2 || error("Usage: julia compare_catalogs.jl relocated_catalog.txt truelocs.txt [output_prefix]")
    length(args) <= 3 || error("Usage: julia compare_catalogs.jl relocated_catalog.txt truelocs.txt [output_prefix]")
    catalog_path, truth_path = args[1], args[2]
    prefix = length(args) >= 3 ? args[3] : joinpath(dirname(abspath(catalog_path)), "yifan_benchmark")
    catalog = numeric_matrix(catalog_path)
    truth_all = numeric_matrix(truth_path)
    size(catalog, 2) >= 9 || error("Relocated catalog needs standard latitude/longitude/depth columns 7/8/9")
    size(truth_all, 2) >= 3 || error("Truth file needs latitude, longitude, and depth")
    ids = round.(Int, catalog[:, end])
    all(id -> 1 <= id <= size(truth_all, 1), ids) || error("Catalog serial IDs must index rows of truelocs.txt")
    length(unique(ids)) == length(ids) || error("Relocated catalog contains duplicate serial IDs")
    truth = truth_all[ids, 1:3]
    n = length(ids)

    horizontal_error = [horizontal_distance_km(truth[i, 1], truth[i, 2], catalog[i, 7], catalog[i, 8]) for i in 1:n]
    depth_error = abs.(catalog[:, 9] .- truth[:, 3])

# The authors' public code selects a pair when both horizontal and vertical
# truth separations are individually below 2 km, then reports the mean of each
# event's RMS neighbor-distance misfit. Every unique pair contributes to both
# endpoint event summaries.
    neighbors = [Tuple{Int,Float64,Float64}[] for _ in 1:n]
    neighboring_pairs = 0
    for i in 1:n-1, j in i+1:n
        truth_horizontal = horizontal_distance_km(truth[i, 1], truth[i, 2], truth[j, 1], truth[j, 2])
        truth_vertical = abs(truth[i, 3] - truth[j, 3])
        truth_horizontal < 2.0 && truth_vertical < 2.0 || continue
        push!(neighbors[i], (j, truth_horizontal, truth_vertical))
        push!(neighbors[j], (i, truth_horizontal, truth_vertical))
        neighboring_pairs += 1
    end
    all(list -> !isempty(list), neighbors) || error("At least one event has no truth neighbor within the benchmark's 2 km horizontal/vertical limits")

    precision_horizontal, precision_depth = zeros(n), zeros(n)
    for i in 1:n
        horizontal_sum2, depth_sum2 = 0.0, 0.0
        for (j, truth_horizontal, truth_vertical) in neighbors[i]
            output_horizontal = horizontal_distance_km(catalog[i, 7], catalog[i, 8], catalog[j, 7], catalog[j, 8])
            output_vertical = abs(catalog[i, 9] - catalog[j, 9])
            horizontal_sum2 += (output_horizontal - truth_horizontal)^2
            depth_sum2 += (output_vertical - truth_vertical)^2
        end
        precision_horizontal[i] = sqrt(horizontal_sum2 / length(neighbors[i]))
        precision_depth[i] = sqrt(depth_sum2 / length(neighbors[i]))
    end

    tx, ty, tz = benchmark_enu_km(truth[:, 1], truth[:, 2], truth[:, 3])
    ox, oy, oz = benchmark_enu_km(catalog[:, 7], catalog[:, 8], catalog[:, 9])
    metrics = (
        mean_horizontal_accuracy_km=average(horizontal_error),
        mean_depth_accuracy_km=average(depth_error),
        median_horizontal_accuracy_km=middle(horizontal_error),
        median_depth_accuracy_km=middle(depth_error),
        chamfer_distance_km=chamfer_distance(ox, oy, oz, tx, ty, tz),
        mean_horizontal_precision_km=average(precision_horizontal),
        mean_depth_precision_km=average(precision_depth),
        neighboring_pairs=neighboring_pairs,
    )

    @printf("Yifan benchmark (%d events; %d unique truth-neighbor pairs)\n", n, metrics.neighboring_pairs)
    @printf("  mean accuracy:   horizontal %.6f km; depth %.6f km\n", metrics.mean_horizontal_accuracy_km, metrics.mean_depth_accuracy_km)
    @printf("  median accuracy: horizontal %.6f km; depth %.6f km\n", metrics.median_horizontal_accuracy_km, metrics.median_depth_accuracy_km)
    @printf("  mean precision:  horizontal %.6f km; depth %.6f km\n", metrics.mean_horizontal_precision_km, metrics.mean_depth_precision_km)
    @printf("  Chamfer distance: %.6f km\n", metrics.chamfer_distance_km)

    mkpath(dirname(prefix))
    open(prefix * "_metrics.csv", "w") do io
        println(io, "events,neighboring_pairs,mean_horizontal_accuracy_km,mean_depth_accuracy_km,median_horizontal_accuracy_km,median_depth_accuracy_km,chamfer_distance_km,mean_horizontal_precision_km,mean_depth_precision_km")
        @printf(io, "%d,%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n", n, metrics.neighboring_pairs,
            metrics.mean_horizontal_accuracy_km, metrics.mean_depth_accuracy_km, metrics.median_horizontal_accuracy_km,
            metrics.median_depth_accuracy_km, metrics.chamfer_distance_km, metrics.mean_horizontal_precision_km, metrics.mean_depth_precision_km)
    end
    open(prefix * "_event_errors.csv", "w") do io
        println(io, "EventID,horizontal_accuracy_error_km,depth_accuracy_error_km,horizontal_precision_error_km,depth_precision_error_km,truth_neighbor_count")
        for i in eachindex(ids)
            @printf(io, "%d,%.12g,%.12g,%.12g,%.12g,%d\n", ids[i], horizontal_error[i], depth_error[i],
                precision_horizontal[i], precision_depth[i], length(neighbors[i]))
        end
    end
    return metrics
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
