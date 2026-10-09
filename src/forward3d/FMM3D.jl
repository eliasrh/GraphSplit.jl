# Grid geometry, interpolation and upwind update adapted from TomoSplit.jl.
# See TOMOSPLIT_LICENSE and docs/THREE_DIMENSIONAL.md for provenance.
module FMM3D
_time_less(t1, n1, t2, n2) = t1 < t2 || (t1 == t2 && n1 < n2)
"""
    Grid3D(x, y, z)

Cartesian node coordinates in metres. Each increasing axis must have at least
two nodes and uniform spacing; the three axis spacings may differ. The depth
axis convention is set by the caller (TomoSplit uses positive-down depth).
"""
struct Grid3D
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    function Grid3D(x, y, z)
        axes = map(a -> Float64.(collect(a)), (x, y, z))
        for a in axes
            length(a) >= 2 || throw(ArgumentError("each grid axis needs at least two nodes"))
            all(isfinite, a) || throw(ArgumentError("grid coordinates must be finite"))
            h = a[2] - a[1]
            h > 0 || throw(ArgumentError("grid axes must increase strictly"))
            all(d -> d > 0 && isapprox(d, h; rtol=1e-10, atol=1e-10 * h), diff(a)) ||
                throw(ArgumentError("forward grid axes must be uniformly spaced"))
        end
        new(axes...)
    end
end

Base.size(g::Grid3D) = (length(g.x), length(g.y), length(g.z))
Base.length(g::Grid3D) = prod(size(g))
grid_spacing(g::Grid3D) = (g.x[2] - g.x[1], g.y[2] - g.y[1], g.z[2] - g.z[1])

"""Nodal P and S velocities in metres/second on a `Grid3D`."""
struct Model3D
    grid::Grid3D
    vp::Array{Float64,3}
    vs::Array{Float64,3}
    function Model3D(grid::Grid3D, vp::AbstractArray, vs::AbstractArray)
        size(vp) == size(vs) == size(grid) || throw(DimensionMismatch("velocity arrays must match the grid"))
        all(v -> isfinite(v) && v > 0, vp) && all(v -> isfinite(v) && v > 0, vs) ||
            throw(ArgumentError("velocities must be positive and finite"))
        new(grid, Array{Float64,3}(vp), Array{Float64,3}(vs))
    end
end

function _cell_fraction(a::Vector{Float64}, q::Real)
    # Unit conversion (metres -> kilometres -> metres) can move an exact
    # physical boundary by one ulp. Canonicalize only this floating-point
    # ambiguity; this is not a distance-based extrapolation tolerance.
    query=Float64(q); lower=first(a); upper=last(a)
    if isfinite(query)
        query < lower && lower-query <= 8eps(lower) && (query=lower)
        query > upper && query-upper <= 8eps(upper) && (query=upper)
    end
    isfinite(query) && lower <= query <= upper ||
        throw(DomainError(q, "point lies outside the forward grid [$(first(a)), $(last(a))]"))
    # Julia's total-order search distinguishes -0.0 from +0.0 although the
    # physical bounds comparison does not. Surface depths from -elevation can
    # be negative zero, so keep a valid boundary cell after the strict check.
    i = clamp(searchsortedlast(a, query), 1, length(a)-1)
    (i, (query-a[i])/(a[i+1]-a[i]))
end

function _interpolation(g::Grid3D, xyz)
    length(xyz) == 3 || throw(ArgumentError("a Cartesian point needs three coordinates"))
    i, u = _cell_fraction(g.x, xyz[1])
    j, v = _cell_fraction(g.y, xyz[2])
    k, w = _cell_fraction(g.z, xyz[3])
    nx, ny, _ = size(g)
    ids = ntuple(8) do n
        a = (n-1) & 1; b = ((n-1) >> 1) & 1; c = ((n-1) >> 2) & 1
        (i+a) + nx*((j+b)-1) + nx*ny*((k+c)-1)
    end
    values = ntuple(8) do n
        a = (n-1) & 1; b = ((n-1) >> 1) & 1; c = ((n-1) >> 2) & 1
        (a == 0 ? 1-u : u)*(b == 0 ? 1-v : v)*(c == 0 ? 1-w : w)
    end
    ids, values, (u, v, w)
end

"""Nonzero trilinear interpolation weights, keyed by linear grid node index."""
function interpolation_weights(g::Grid3D, xyz)
    ids, ws, _ = _interpolation(g, xyz)
    Pair{Int,Float64}[ids[a] => ws[a] for a in 1:8 if ws[a] != 0.0]
end

function _source_axis_nodes(a, q)
    i, u = _cell_fraction(a, q)
    if u == 0
        return max(1, i-1):min(length(a), i+1)
    elseif u == 1
        return max(1, i):min(length(a), i+2)
    end
    i:i+1
end

# Coordinate transforms can leave a nominal grid-aligned station a few ulps
# to one side of a grid plane. Canonicalize that numerical ambiguity so the
# adjacent-cell source stencil cannot depend on a trigonometric round-off bit.
# Bounds are checked before this helper is called, with only an eight-ulp
# round-off allowance. Genuine sub-cell offsets retain the normal source stencil.
function _canonical_source_coordinate(a, q)
    i, _ = _cell_fraction(a,q)
    left,right=a[i],a[i+1]
    tolerance=64eps(max(abs(q),abs(left),abs(right),right-left))
    abs(q-left) <= tolerance && return left
    abs(q-right) <= tolerance && return right
    q
end

# All source-to-seed segments stay in one grid cell. Trilinear slowness along
# the segment is cubic, so two-point Gauss integration is exact for that model.
function _seed_segment(g, slow, station, point)
    d = sqrt(sum((point[a]-station[a])^2 for a in 1:3))
    coeff = Dict{Int,Float64}()
    if d > 0
        for q in (0.5 - 0.5/sqrt(3.0), 0.5 + 0.5/sqrt(3.0))
            p = ntuple(a -> clamp(muladd(q, point[a]-station[a], station[a]),
                                   min(station[a],point[a]), max(station[a],point[a])), 3)
            ids, ws, _ = _interpolation(g, p)
            for a in 1:8
                ws[a] == 0 && continue
                coeff[ids[a]] = get(coeff, ids[a], 0.0) + 0.5*d*ws[a]
            end
        end
    end
    deriv = sort!(collect(coeff); by=first)
    sum((last(p)*slow[first(p)] for p in deriv); init=0.0), deriv
end

# Accepted neighbours on opposite sides of an axis share the same spacing.
# Their smaller time is the Godunov upwind neighbour (ties are deterministic).
function _upwind_update(n, slow, times, accepted, dims, spacing, ::Val{O}=Val(1)) where O
    nx, ny, nz = dims
    i = mod(n-1, nx)+1; j = mod(div(n-1, nx), ny)+1; k = div(n-1, nx*ny)+1
    neighbours = ((i > 1 ? n-1 : 0, i < nx ? n+1 : 0),
                  (j > 1 ? n-nx : 0, j < ny ? n+nx : 0),
                  (k > 1 ? n-nx*ny : 0, k < nz ? n+nx*ny : 0))
    raw_candidates = ntuple(3) do a
        l, r = neighbours[a]
        l = l > 0 && accepted[l] ? l : 0
        r = r > 0 && accepted[r] ? r : 0
        p = l == 0 ? r : r == 0 ? l : times[l] <= times[r] ? l : r
        if p == 0
            return (Inf, 0, 0.0, 0)
        end
        t1 = times[p]
        ih2 = inv(spacing[a]^2)
        if O == 2
            coord = (i,j,k)[a]
            stride = (1,nx,nx*ny)[a]
            second = p == n-stride ? (coord > 2 ? n-2stride : 0) :
                                      (coord < dims[a]-1 ? n+2stride : 0)
            if second > 0 && accepted[second] && times[second] <= t1
                # (3T - 4T1 + T2)/(2h) = (T - Teff)/heff.
                # Teff >= T1 >= T2 enforces acceptance-order causality. The
                # second node has a negative derivative, so this is not a
                # monotone scheme in all input neighbour values.
                return (t1 + (t1-times[second])/3, p, 2.25ih2, second)
            end
        end
        (t1, p, ih2, 0)
    end
    # A three-entry sorting network keeps the inner marching update allocation
    # free. Missing neighbours sort last and never enter the quadratic.
    c1, c2, c3 = raw_candidates
    _time_less(c2[1], c2[2], c1[1], c1[2]) && ((c1,c2) = (c2,c1))
    _time_less(c3[1], c3[2], c2[1], c2[2]) && ((c2,c3) = (c3,c2))
    _time_less(c2[1], c2[2], c1[1], c1[2]) && ((c1,c2) = (c2,c1))
    candidates = (c1,c2,c3)
    count = (c1[2] != 0) + (c2[2] != 0) + (c3[2] != 0)
    if count == 0
        return O == 1 ? (Inf, (0,0,0), (0.0,0.0,0.0), 0.0) :
                        (Inf, (0,0,0,0,0,0), (0.0,0.0,0.0,0.0,0.0,0.0), 0.0)
    end
    # Shift the quadratic to the smallest neighbour time to avoid cancellation
    # when absolute travel time is large relative to one grid-cell crossing.
    base = candidates[1][1]
    aa = 0.0; bb = 0.0; cc = 0.0; t = Inf; active = 0
    for a in 1:count
        ti, _, wi, _ = candidates[a]
        d = ti-base
        aa += wi; bb += wi*d; cc += wi*d*d
        disc = bb*bb - aa*(cc-slow[n]^2)
        t = base + (bb + sqrt(max(0.0, disc)))/aa
        active = a
        (a == count || t <= candidates[a+1][1]) && break
    end
    d1 = (t-c1[1])*c1[3]
    d2 = active >= 2 ? (t-c2[1])*c2[3] : 0.0
    d3 = active >= 3 ? (t-c3[1])*c3[3] : 0.0
    den = d1+d2+d3
    den > 0 || error("noncausal or degenerate fast-marching update")
    if O == 1
        ps = (c1[2], active >= 2 ? c2[2] : 0, active >= 3 ? c3[2] : 0)
        ws = (d1/den,d2/den,d3/den)
        return t, ps, ws, slow[n]/den
    end
    p1, p2, p3 = c1[2], active >= 2 ? c2[2] : 0, active >= 3 ? c3[2] : 0
    q1, q2, q3 = c1[4], active >= 2 ? c2[4] : 0, active >= 3 ? c3[4] : 0
    w1, w2, w3 = d1/den,d2/den,d3/den
    ps = (p1,p2,p3,q1,q2,q3)
    ws = (q1 == 0 ? w1 : 4w1/3, q2 == 0 ? w2 : 4w2/3, q3 == 0 ? w3 : 4w3/3,
          q1 == 0 ? 0.0 : -w1/3, q2 == 0 ? 0.0 : -w2/3, q3 == 0 ? 0.0 : -w3/3)
    t, ps, ws, slow[n]/den
end

include("Marching.jl")
end # module FMM3D
