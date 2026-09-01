const DEG2RAD = pi / 180.0
const RAD2DEG = 180.0 / pi

function circular_mean_longitude(lon_deg::AbstractVector{<:Real})
    angles = Float64.(lon_deg) .* DEG2RAD
    return atan(sum(sin, angles), sum(cos, angles)) * RAD2DEG
end

function local_xy(lat, lon, lat0::Real, lon0::Real, radius_m::Real)
    phi = Float64.(lat) .* DEG2RAD
    lambda = Float64.(lon) .* DEG2RAD
    phi0 = Float64(lat0) * DEG2RAD
    lambda0 = Float64(lon0) * DEG2RAD
    dlambda = atan.(sin.(lambda .- lambda0), cos.(lambda .- lambda0))
    return Float64(radius_m) .* cos(phi0) .* dlambda, Float64(radius_m) .* (phi .- phi0)
end

function local_xy_to_ll(x, y, lat0::Real, lon0::Real, radius_m::Real)
    phi0 = Float64(lat0) * DEG2RAD
    lambda0 = Float64(lon0) * DEG2RAD
    phi = Float64.(y) ./ Float64(radius_m) .+ phi0
    lambda = Float64.(x) ./ (Float64(radius_m) * cos(phi0)) .+ lambda0
    lambda = atan.(sin.(lambda), cos.(lambda))
    return phi .* RAD2DEG, lambda .* RAD2DEG
end

function choose_reference(stations::Stations, catalog::Catalog, cfg::AbstractDict)
    mode = lowercase(String(cfgget(cfg, "coordinates", "reference"; default="stations_mean")))
    if mode == "stations_mean"
        return sum(stations.lat) / length(stations.lat), circular_mean_longitude(stations.lon)
    elseif mode == "catalog_mean"
        return sum(catalog.lat) / length(catalog.lat), circular_mean_longitude(catalog.lon)
    elseif mode == "manual"
        return Float64(cfgget(cfg, "coordinates", "reference_latitude")), Float64(cfgget(cfg, "coordinates", "reference_longitude"))
    end
    error("Unknown coordinates.reference: $mode")
end

function catalog_depth_to_internal(depth_km, cfg::AbstractDict)
    vertical = lowercase(String(cfgget(cfg, "coordinates", "event_vertical"; default="positive_depth")))
    z0 = Float64(cfgget(cfg, "coordinates", "event_z0_m"; default=0.0))
    depth_m = Float64.(depth_km) .* 1000.0
    return vertical == "positive_depth" ? depth_m :
        vertical == "negative_depth" ? -depth_m :
        vertical == "positive_depth_plus_z0" ? depth_m .+ z0 :
        error("Unknown coordinates.event_vertical: $vertical")
end

function internal_depth_to_km(value::Real, cfg::AbstractDict)
    vertical = lowercase(String(cfgget(cfg, "coordinates", "event_vertical"; default="positive_depth")))
    z0 = Float64(cfgget(cfg, "coordinates", "event_z0_m"; default=0.0))
    return vertical == "positive_depth" ? Float64(value) / 1000.0 :
        vertical == "negative_depth" ? -Float64(value) / 1000.0 :
        vertical == "positive_depth_plus_z0" ? (Float64(value) - z0) / 1000.0 :
        error("Unknown coordinates.event_vertical: $vertical")
end

"Attach a reproducible local coordinate frame and return the initial relocation state."
function attach_coordinates!(stations::Stations, catalog::Catalog, cfg::AbstractDict)
    lat0, lon0 = choose_reference(stations, catalog, cfg)
    radius = Float64(cfgget(cfg, "coordinates", "earth_radius_m"; default=6_371_000.0))
    isfinite(lat0) && isfinite(lon0) || error("Coordinate reference is not finite")
    abs(cos(lat0 * DEG2RAD)) > 1.0e-10 || error("Coordinate reference is too close to a pole")

    stations.x_m, stations.y_m = local_xy(stations.lat, stations.lon, lat0, lon0, radius)
    catalog_x, catalog_y = local_xy(catalog.lat, catalog.lon, lat0, lon0, radius)
    station_vertical = lowercase(String(cfgget(cfg, "coordinates", "station_vertical"; default="depth_from_elevation")))
    stations.z_m = station_vertical == "depth_from_elevation" ? -stations.elev_m :
        station_vertical == "elevation" ? copy(stations.elev_m) :
        error("Unknown coordinates.station_vertical: $station_vertical")

    event_z = catalog_depth_to_internal(catalog.depth_km, cfg)

    stations.ref_lat = catalog.ref_lat = lat0
    stations.ref_lon = catalog.ref_lon = lon0
    stations.ref_radius_m = catalog.ref_radius_m = radius
    return State(copy(catalog.event_id), catalog_x, catalog_y, event_z, zeros(length(catalog)), lat0, lon0, radius)
end

"Override catalog hypocenters with one common Stage-1 seed when requested."
function apply_initialization!(state::State, cfg::AbstractDict)
    mode = lowercase(String(cfgget(cfg, "initialization", "mode"; default="catalog")))
    mode == "catalog" && return state
    if mode == "common_centroid"
        n = length(state)
        common_x = sum(state.x) / n
        common_y = sum(state.y) / n
        common_z = sum(state.z) / n
    elseif mode == "common_manual"
        latitude = Float64(cfgget(cfg, "initialization", "latitude"))
        longitude = Float64(cfgget(cfg, "initialization", "longitude"))
        common_x, common_y = local_xy(latitude, longitude,
            state.ref_lat, state.ref_lon, state.ref_radius_m)
        common_z = catalog_depth_to_internal(Float64(cfgget(cfg, "initialization", "depth_km")), cfg)
    else
        error("Unknown initialization.mode: $mode")
    end
    fill!(state.x, common_x)
    fill!(state.y, common_y)
    fill!(state.z, common_z)
    fill!(state.t0, 0.0)
    return state
end

"Exact spherical central angle and derivatives with respect to event local x/y."
function spherical_delta_gradient(x::Real, y::Real, sx::Real, sy::Real,
        ref_lat_deg::Real, ref_lon_deg::Real, coordinate_radius_m::Real)
    radius = Float64(coordinate_radius_m)
    phi0 = Float64(ref_lat_deg) * DEG2RAD
    lambda0 = Float64(ref_lon_deg) * DEG2RAD
    cosphi0 = cos(phi0)
    phi = phi0 + Float64(y) / radius
    phis = phi0 + Float64(sy) / radius
    lambda = lambda0 + Float64(x) / (radius * cosphi0)
    lambdas = lambda0 + Float64(sx) / (radius * cosphi0)
    dlambda = atan(sin(lambda - lambdas), cos(lambda - lambdas))
    sp, cp = sincos(phi)
    sps, cps = sincos(phis)
    u = clamp(sp * sps + cp * cps * cos(dlambda), -1.0, 1.0)
    sindelta = sqrt(max(0.0, 1.0 - u * u))
    delta = atan(sindelta, u)
    if sindelta > 1.0e-10
        du_dphi = cp * sps - sp * cps * cos(dlambda)
        du_dlambda = -cp * cps * sin(dlambda)
        return delta, (-du_dlambda / sindelta) / (radius * cosphi0), (-du_dphi / sindelta) / radius
    end
    # Coincident points have a direction-dependent derivative. A centered metre-scale
    # difference is stable for the rare exact-coincidence case.
    h = 1.0
    dx = (spherical_delta_only(x + h, y, sx, sy, phi0, lambda0, radius) -
          spherical_delta_only(x - h, y, sx, sy, phi0, lambda0, radius)) / (2h)
    dy = (spherical_delta_only(x, y + h, sx, sy, phi0, lambda0, radius) -
          spherical_delta_only(x, y - h, sx, sy, phi0, lambda0, radius)) / (2h)
    return delta, dx, dy
end

function spherical_delta_only(x, y, sx, sy, phi0, lambda0, radius)
    cosphi0 = cos(phi0)
    phi = phi0 + y / radius
    phis = phi0 + sy / radius
    lambda = lambda0 + x / (radius * cosphi0)
    lambdas = lambda0 + sx / (radius * cosphi0)
    dlambda = atan(sin(lambda - lambdas), cos(lambda - lambdas))
    u = clamp(sin(phi) * sin(phis) + cos(phi) * cos(phis) * cos(dlambda), -1.0, 1.0)
    return atan(sqrt(max(0.0, 1.0 - u * u)), u)
end
