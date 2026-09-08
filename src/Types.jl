"A numeric earthquake catalog whose final column is the DDSync serial ID."
mutable struct Catalog
    file::String
    raw::Matrix{Float64}
    event_id::Vector{Int64}
    lat::Vector{Float64}
    lon::Vector{Float64}
    depth_km::Vector{Float64}
    lat_col::Int
    lon_col::Int
    depth_col::Int
    id_col::Int
    ref_lat::Float64
    ref_lon::Float64
    ref_radius_m::Float64
end

"Station metadata and its local east/north/depth-positive coordinates."
mutable struct Stations
    id::Vector{String}
    lat::Vector{Float64}
    lon::Vector{Float64}
    elev_m::Vector{Float64}
    x_m::Vector{Float64}
    y_m::Vector{Float64}
    z_m::Vector{Float64}
    ref_lat::Float64
    ref_lon::Float64
    ref_radius_m::Float64
end

"The current four-parameter model for every event."
mutable struct State
    event_id::Vector{Int64}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    t0::Vector{Float64}
    ref_lat::Float64
    ref_lon::Float64
    ref_radius_m::Float64
end

"Cumulative catalog displacement in the GraphSplit local frame, keyed by event ID."
struct CatalogShift
    event_id::Vector{Int64}
    dx_m::Vector{Float64}
    dy_m::Vector{Float64}
    dz_m::Vector{Float64}
    dt0_s::Vector{Float64}
end

Base.length(c::Catalog) = length(c.event_id)
Base.length(s::Stations) = length(s.id)
Base.length(s::State) = length(s.event_id)
Base.length(s::CatalogShift) = length(s.event_id)

"One finite DDSync theta entry, keyed by serial event ID."
struct ThetaEntry
    theta::Float64
    ref_id::Int64
    sigma::Float64
    degree::Int32
end

"One DDSync station-phase synchronization group."
mutable struct ThetaGroup
    station::String
    phase::UInt8             # 1=P, 2=S
    name::String
    event_id::Vector{Int64}  # input order, for deterministic allocation-free scans
    entries::Dict{Int64,ThetaEntry}
    has_sigma::Bool
    has_degree::Bool
end

"Compact row-oriented observations for Stage 1 or Stage 2."
mutable struct Observations
    i::Vector{Int32}
    j::Vector{Int32}
    dt::Vector{Float64}
    sigma::Vector{Float64}
    station::Vector{Int32}
    phase::Vector{UInt8}
    group::Vector{Int32}
end

Observations() = Observations(Int32[], Int32[], Float64[], Float64[], Int32[], UInt8[], Int32[])
Base.length(o::Observations) = length(o.dt)
Base.isempty(o::Observations) = isempty(o.dt)

"Sparse undirected event graph plus support diagnostics."
struct EventGraph
    i::Vector{Int32}
    j::Vector{Int32}
    dist_metric_m::Vector{Float64}
    dist_xy_m::Vector{Float64}
    dist_xyz_m::Vector{Float64}
    degree::Vector{Int32}
    comp_id::Vector{Int32}
    comp_size::Vector{Int32}
    ncomp::Int
    is_augmented::BitVector
end

Base.length(g::EventGraph) = length(g.i)

"Per-outer-iteration solver diagnostics."
struct SolveStats
    iterations::Int
    rms_s::Vector{Float64}
    step_rms_m::Vector{Float64}
    step_rms_s::Vector{Float64}
    inner_iterations::Vector{Int}
    converged::Bool
    depth_bound_hits::Vector{Int}
    depth_bound_active::BitVector
end

abstract type AbstractTravelTimeModel end

"Constant velocity model used for tests and diagnostics."
struct ConstantTravelTime <: AbstractTravelTimeModel
    geometry::Symbol
    vp_ms::Float64
    vs_ms::Float64
    earth_radius_m::Float64
    ref_lat::Float64
    ref_lon::Float64
    coordinate_radius_m::Float64
end

"Memory-mapped native GraphSplit travel-time table."
mutable struct TravelTimeTable <: AbstractTravelTimeModel
    file::String
    geometry::Symbol
    earth_radius_m::Float64
    coordinate_radius_m::Float64
    model_hash::Vector{UInt8}
    q::Vector{Float64}
    z::Vector{Float64}
    zs::Vector{Float64}
    p::Array{Float32,3}
    s::Array{Float32,3}
    ref_lat::Float64
    ref_lon::Float64
    clamp::Bool
    io::IOStream
end
