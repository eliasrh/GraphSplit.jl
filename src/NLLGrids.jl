"Geographic definition attached to a NonLinLoc grid; distances inside GraphSplit are metres."
struct GridProjection3D
    kind::Symbol
    latitude::Float64
    longitude::Float64
    rotation::Float64
end
Base.:(==)(a::GridProjection3D,b::GridProjection3D)=
    a.kind==b.kind&&a.latitude==b.latitude&&a.longitude==b.longitude&&a.rotation==b.rotation

struct NLLHeader3D
    file::String
    dimensions::NTuple{3,Int}
    origin_m::NTuple{3,Float64}
    spacing_m::NTuple{3,Float64}
    kind::String
    precision::DataType
    projection::GridProjection3D
    station::Union{Nothing,Tuple{String,NTuple{3,Float64}}}
end
nll_buffer_path(h::NLLHeader3D)=h.file[1:end-4]*".buf"
nll_bytes(h::NLLHeader3D)=prod(Int128.(h.dimensions))*sizeof(h.precision)

"Read metadata only; validate dimensions and buffer length before allocating a grid."
function read_nll_header(path::AbstractString; coordinate_system::String="header")
    file=endswith(lowercase(path),".hdr") ? String(path) : String(path)*".hdr"
    lines=filter(!isempty,[strip(first(split(l,'#'))) for l in eachline(file)])
    isempty(lines)&&error("Empty NonLinLoc header: $file")
    fields=split(lines[1]);length(fields) in (10,11)||error("NonLinLoc header needs ten fields, optionally followed by FLOAT or DOUBLE: $file")
    dims=Tuple(parse.(Int,fields[1:3]));all(>=(2),dims)||error("Only full 3D NonLinLoc grids are supported (each dimension >= 2); TIME2D is not supported")
    prod(Int128.(dims))<=typemax(Int32)||error("NonLinLoc grid exceeds the 32-bit node-index limit")
    origin=Tuple(1000 .* parse.(Float64,fields[4:6]));step=Tuple(1000 .* parse.(Float64,fields[7:9]))
    all(isfinite,origin)&&all(h->isfinite(h)&&h>0,step)||error("Grid origin and spacing must be finite; spacing must be positive")
    kind=uppercase(fields[10]);precision=length(fields)==10 ? "FLOAT" : uppercase(fields[11])
    precision in ("FLOAT","DOUBLE")||error("NonLinLoc precision must be FLOAT or DOUBLE")
    T=precision=="FLOAT" ? Float32 : Float64
    trans=filter(l->startswith(l,"TRANSFORM ")||startswith(l,"TRANS "),lines[2:end])
    length(trans)<=1||error("Multiple geographic transforms in $file")
    tokens=isempty(trans) ? String[] : split(only(trans))
    if coordinate_system=="local"
        (isempty(tokens)||tokens[2]=="NONE")||error("coordinate_system=local is only allowed for a missing/NONE transform; do not ignore a declared geographic projection")
        projection=GridProjection3D(:local,0.,0.,0.)
    else
        length(tokens)>=2||error("No geographic transform in $file; for a deliberately local grid set grid3d.coordinate_system=local and a manual coordinate reference")
        tokens[2]=="SIMPLE"||error("This experimental reader supports NonLinLoc SIMPLE only. $(tokens[2]) is not silently approximated; re-export in SIMPLE or supply a documented local grid.")
        if length(tokens)==8&&tokens[3]=="LatOrig"&&tokens[5]=="LongOrig"&&tokens[7]=="RotCW"
            lat,lon,rot=parse.(Float64,tokens[[4,6,8]])
        elseif length(tokens)==5&&tokens[1]=="TRANS"
            lat,lon,rot=parse.(Float64,tokens[3:5])
        else
            error("Unrecognized NonLinLoc SIMPLE transform in $file")
        end
        all(isfinite,(lat,lon,rot))&&abs(lat)<90&&abs(lon)<=180||error("Invalid SIMPLE origin/rotation")
        projection=GridProjection3D(:nll_simple,lat,lon,rot)
    end
    station=nothing
    if kind=="TIME"
        length(lines)>=2||error("TIME header is missing its station line")
        q=split(lines[2]);length(q)==4||error("TIME station line must be STA x_km y_km depth_km")
        xyz=Tuple(1000 .* parse.(Float64,q[2:4]));all(isfinite,xyz)||error("Non-finite TIME station location")
        station=(String(q[1]),xyz)
    end
    header=NLLHeader3D(abspath(file),dims,origin,step,kind,T,projection,station)
    isfile(nll_buffer_path(header))||error("Missing NonLinLoc .buf file for $file")
    filesize(nll_buffer_path(header))==nll_bytes(header)||error("NonLinLoc buffer size disagrees with its dimensions/precision: $file")
    header
end

function nll_grid(h::NLLHeader3D)
    FMM3D.Grid3D((collect(range(h.origin_m[a];step=h.spacing_m[a],length=h.dimensions[a])) for a in 1:3)...)
end
function same_nll_grid(a::NLLHeader3D,b::NLLHeader3D)
    a.dimensions==b.dimensions&&a.origin_m==b.origin_m&&a.spacing_m==b.spacing_m&&a.projection==b.projection
end

"NLL buffers have z varying fastest, then y, then x; endian choice is explicit."
function with_nll_buffer(f,h::NLLHeader3D,byte_order::String)
    native_little=ENDIAN_BOM==0x04030201
    swap=byte_order=="native" ? false : (byte_order=="little")!=native_little
    open(nll_buffer_path(h),"r") do io
        data=Mmap.mmap(io,Array{h.precision,3},reverse(h.dimensions))
        value=(i,j,k)->begin
            v=data[k,j,i]
            if swap
                v=h.precision===Float32 ? reinterpret(Float32,bswap(reinterpret(UInt32,v))) : reinterpret(Float64,bswap(reinterpret(UInt64,v)))
            end
            Float64(v)
        end
        f(value)
    end
end

function nll_velocity(value,h::NLLHeader3D)
    isfinite(value)&&value>0||error("NonLinLoc model contains a nonpositive or non-finite sample. Supply positive velocities throughout the box; use the separate surface file to mask air.")
    h.kind=="VELOCITY"&&return 1000value
    h.kind=="VELOCITY_METERS"&&return value
    h.kind=="SLOWNESS"&&return 1000/value
    h.kind=="SLOW_LEN"&&return h.spacing_m[1]/value
    error("Unsupported NonLinLoc velocity type $(h.kind); use VELOCITY, VELOCITY_METERS, SLOWNESS or SLOW_LEN")
end

function sample_nll_velocity(h::NLLHeader3D,g::FMM3D.Grid3D,method,byte_order; sampling="cell_centers")
    h.kind in ("VELOCITY","VELOCITY_METERS","SLOWNESS","SLOW_LEN")||error("Expected a velocity/slowness model, not $(h.kind)")
    h.kind=="SLOW_LEN"&&!all(d->isapprox(d,h.spacing_m[1];rtol=1e-10),h.spacing_m)&&
        error("SLOW_LEN requires equal x/y/z spacing; use VELOCITY for unequal spacing")
    sampling in ("cell_centers","nodes")||error("NLL model_sampling must be cell_centers or nodes")
    method in ("nearest","slowness_linear")||error("Unknown velocity interpolation method: $method")
    centered=sampling=="cell_centers"
    # Vel2Grid samples at +half a cell. The final plane on each axis is
    # padding, outside the Grid2Time model. TIME fields themselves are nodal.
    count=ntuple(a->h.dimensions[a]-(centered ? 1 : 0),3)
    offset=centered ? .5 : 0.
    domain=nll_grid(h);out=Array{Float64}(undef,size(g))
    with_nll_buffer(h,byte_order) do value
        for k in eachindex(g.z),j in eachindex(g.y),i in eachindex(g.x)
            xyz=(g.x[i],g.y[j],g.z[k])
            fractional=ntuple(a->begin
                FMM3D._cell_fraction((domain.x,domain.y,domain.z)[a],xyz[a])
                # Within the boundary half-cell use its sole adjacent cell
                # value. The preceding domain check forbids extension outside
                # the physical model box.
                clamp((xyz[a]-h.origin_m[a])/h.spacing_m[a]-offset,0.,count[a]-1.)
            end,3)
            if method=="nearest"
                indices=ntuple(a->floor(Int,fractional[a]+.5)+1,3)
                out[i,j,k]=nll_velocity(value(indices...),h)
            else
                lower=ntuple(a->floor(Int,fractional[a])+1,3)
                upper=ntuple(a->min(lower[a]+1,count[a]),3)
                t=ntuple(a->fractional[a]-(lower[a]-1),3)
                slow=0.
                for c in 0:1,b in 0:1,a in 0:1
                    weight=(a==0 ? 1-t[1] : t[1])*(b==0 ? 1-t[2] : t[2])*(c==0 ? 1-t[3] : t[3])
                    weight==0&&continue
                    index=(a==0 ? lower[1] : upper[1],b==0 ? lower[2] : upper[2],c==0 ? lower[3] : upper[3])
                    slow+=weight/nll_velocity(value(index...),h)
                end
                out[i,j,k]=1/slow
            end
        end
    end
    out
end

# NonLinLoc geo.h AVG_ERAD (GMT Sphere), not GraphSplit's default 6371 km.
const NLL_SIMPLE_RADIUS_M = 6_371_008.7714

"NLL SIMPLE (GMT sphere, latitude-dependent longitude scale), with chain derivatives."
function projected_xy_jacobian(p::GridProjection3D,x,y,lat0,lon0,radius)
    p.kind==:local&&return (Float64(x),Float64(y),1.,0.,0.,1.)
    lat,lon=local_xy_to_ll(x,y,lat0,lon0,radius)
    phi=lat*DEG2RAD;delta=atan(sin((lon-p.longitude)*DEG2RAD),cos((lon-p.longitude)*DEG2RAD))
    rn=NLL_SIMPLE_RADIUS_M;xx=rn*cos(phi)*delta;yy=rn*(lat-p.latitude)*DEG2RAD
    s,c=sincos(p.rotation*DEG2RAD);a=rn*cos(phi)/(radius*cos(lat0*DEG2RAD));b=-rn*sin(phi)*delta/radius;d=rn/radius
    (c*xx+s*yy,-s*xx+c*yy,c*a,c*b+s*d,-s*a,-s*b+c*d)
end

"Write nodal velocities in little-endian NLL layout; read with model_sampling=nodes. Arrays are [x,y,z], m/s."
function write_nll_velocity(path::AbstractString,g::FMM3D.Grid3D,velocity;
        latitude=0.,longitude=0.,rotation=0.,local_coordinates=false)
    size(velocity)==size(g)||error("Velocity shape does not match the grid")
    all(v->isfinite(v)&&v>0,velocity)||error("Velocities must be positive and finite")
    stem=endswith(path,".hdr") ? path[1:end-4] : String(path);mkpath(dirname(abspath(stem)))
    open(stem*".hdr","w") do io
        println(io,join(size(g)," ")," ",join([first(g.x),first(g.y),first(g.z)]./1000," ")," ",join(collect(FMM3D.grid_spacing(g))./1000," ")," VELOCITY FLOAT")
        println(io,local_coordinates ? "TRANSFORM NONE" : "TRANSFORM SIMPLE LatOrig $latitude LongOrig $longitude RotCW $rotation")
    end
    open(stem*".buf","w") do io
        for i in eachindex(g.x),j in eachindex(g.y),k in eachindex(g.z)
            write(io,htol(reinterpret(UInt32,Float32(velocity[i,j,k]/1000))))
        end
    end
    stem*".hdr"
end
