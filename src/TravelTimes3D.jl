"Fixed-model 3D times, one Float32 volume for each requested station and phase."
struct TravelTime3D <: AbstractTravelTimeModel
    geometry::Symbol
    grid::FMM3D.Grid3D
    projection::GridProjection3D
    fields::Dict{Tuple{UInt8,Float64,Float64,Float64},Array{Float32,3}}
    valid::BitArray{3}
    ref_lat::Float64
    ref_lon::Float64
    coordinate_radius_m::Float64
    metadata::Dict{String,Any}
end

function validate_grid3d_config(cfg)
    o=cfg["grid3d"]
    o["model_format"] in ("nll_velocity","nll_time")||error("grid3d.model_format must be nll_velocity or nll_time")
    o["coordinate_system"] in ("header","local")||error("grid3d.coordinate_system must be header or local")
    o["coordinate_system"]=="local"&&cfg["coordinates"]["reference"]!="manual"&&
        error("A local 3D grid requires coordinates.reference=manual so that its origin cannot change with the catalog")
    o["byte_order"] in ("little","big","native")||error("grid3d.byte_order must be little, big or native")
    o["model_interpolation"] in ("nearest","slowness_linear")||error("grid3d.model_interpolation must be nearest or slowness_linear")
    s=o["spacing_m"];length(s)==3&&all(x->isfinite(x)&&x>0,s)||error("grid3d.spacing_m needs three positive spacings in metres")
    b=o["bounds_km"];(isempty(b)||(length(b)==6&&all(isfinite,b)&&all(b[2a]>b[2a-1] for a in 1:3)))||error("grid3d.bounds_km must be empty or [xmin,xmax,ymin,ymax,zmin,zmax]")
    o["model_format"]=="nll_time"&&!isempty(b)&&error("Precomputed TIME grids retain their native extent; bounds_km must be empty")
    phases=o["phases"];!isempty(phases)&&all(p->p in ("P","S"),phases)&&length(unique(phases))==length(phases)||error("grid3d.phases must contain P, S or both, without duplicates")
    o["accuracy_order"] in (1,2)||error("grid3d.accuracy_order must be 1 or 2")
    for key in ("maximum_memory_gib","warning_memory_gib")
        isfinite(o[key])&&o[key]>0||error("grid3d.$key must be positive and finite")
    end
    o["warning_memory_gib"]<=o["maximum_memory_gib"]||error("The 3D warning threshold must not exceed the memory limit")
    !isempty(o["cache_dir"])||error("grid3d.cache_dir cannot be empty")
    cfg["coordinates"]["station_vertical"]=="depth_from_elevation"&&cfg["coordinates"]["event_vertical"]=="positive_depth"||
        error("3D grids require station_vertical=depth_from_elevation and event_vertical=positive_depth, relative to the same vertical datum")
    cfg["travel_time"]["geometry"] in ("auto","cartesian","flat")||error("3D grids use Cartesian propagation, not radial geometry")
    cfg["constraints"]["depth_bound"]["reflected_prelocation_restart"]&&error("Reflected-depth branch search is not supported with a bounded 3D grid")
    nothing
end

file_hash3d(file)=open(io->bytes2hex(sha256(io)),file,"r")

"Plan dimensions before allocating axes or model arrays. Spacing is at most the requested value."
function grid3d_plan(h,o)
    lo=collect(h.origin_m);hi=lo.+collect(h.spacing_m).*(collect(h.dimensions).-1)
    b=o["bounds_km"]
    if !isempty(b)
        l=1000 .* Float64.(b[1:2:5]);u=1000 .* Float64.(b[2:2:6])
        all(l.>=lo.-1e-7)&&all(u.<=hi.+1e-7)||error("Requested 3D bounds extend outside the input model; no velocity extrapolation is applied")
        lo=l;hi=u
    end
    dims=o["model_format"]=="nll_time" ? h.dimensions : ntuple(a->begin
        ratio=(hi[a]-lo[a])/o["spacing_m"][a]
        ratio<typemax(Int32)-1||error("Requested 3D spacing produces too many grid nodes")
        max(2,ceil(Int,ratio)+1)
    end,3)
    nodes=prod(Int128.(dims));nodes<=typemax(Int32)||error("3D grid exceeds the 32-bit node-index limit; use coarser spacing or smaller bounds")
    lo,hi,dims,nodes
end

function grid3d_memory_guard(nodes,groups,input_bytes,surface_bytes,o)
    # Retained fields + serial FMM workspace and velocity + temporary conversion,
    # maps/masks and allocator allowance. Counts mapped inputs as resident too.
    field_bytes=Int128(4)*nodes*groups
    estimated=field_bytes+Int128(80)*nodes+input_bytes+Int128(16)*surface_bytes
    gib=Float64(estimated)/2.0^30
    @printf("3D travel times: %d nodes × %d station–phase fields; %.3f GiB retained, %.3f GiB estimated build peak\n",nodes,groups,Float64(field_bytes)/2.0^30,gib)
    gib<=o["maximum_memory_gib"]||error(@sprintf("3D travel-time estimate %.3f GiB exceeds grid3d.maximum_memory_gib=%.3f. Use coarser spacing, fewer phases/stations or a smaller model. This limit excludes relocation arrays.",gib,o["maximum_memory_gib"]))
    gib>=o["warning_memory_gib"]&&@warn "3D travel times require substantial memory; leave additional RAM for observations, solver arrays and other applications" estimated_gib=gib
    Dict{String,Any}("retained_fields_gib"=>Float64(field_bytes)/2.0^30,"estimated_build_peak_gib"=>gib,
        "memory_scope"=>"travel-time fields, model input and serial FMM; excludes catalogs, observations and relocation arrays")
end

"Read a complete regular surface: model x/y in km, elevation in m, positive up."
function read_surface3d(file,g)
    isempty(file)&&return nothing
    rows=NTuple{3,Float64}[]
    for (line_number,line) in enumerate(eachline(file))
        f=split(strip(first(split(line,'#'))));isempty(f)&&continue
        length(f)==3||error("Surface line $line_number must have x_km y_km elevation_m")
        v=parse.(Float64,f);all(isfinite,v)||error("Non-finite surface value at line $line_number")
        push!(rows,(1000v[1],1000v[2],-v[3]))
    end
    isempty(rows)&&error("Empty surface file")
    xs=sort!(unique(first.(rows)));ys=sort!(unique(getindex.(rows,2)))
    length(xs)*length(ys)==length(rows)||error("Surface must be a complete rectangular grid")
    # Grid3D validates regular spacing. A two-level dummy depth axis is sufficient.
    terrain=FMM3D.Grid3D(xs,ys,[0.,1.]);z=fill(NaN,length(xs),length(ys))
    for (x,y,d) in rows
        i=searchsortedfirst(xs,x);j=searchsortedfirst(ys,y)
        isnan(z[i,j])||error("Duplicate surface coordinate");z[i,j]=d
    end
    all(isfinite,z)||error("Surface has missing cells")
    [FMM3D.surface_depth(terrain,z,x,y) for x in g.x,y in g.y]
end

function load_cached_field3d(file,dims)
    filesize(file)==4prod(Int128.(dims))||error("Incomplete 3D travel-time cache: $file")
    open(io->Mmap.mmap(io,Array{Float32,3},dims),file,"r")
end

function read_time_field3d(h,g,surface,byte_order)
    h.kind=="TIME"||error("Expected a NonLinLoc TIME grid, not $(h.kind)")
    result=Array{Float32}(undef,size(g))
    with_nll_buffer(h,byte_order) do value
        for k in eachindex(g.z),j in eachindex(g.y),i in eachindex(g.x)
            v=value(i,j,k)
            if surface!==nothing&&g.z[k]<surface[i,j]
                result[i,j,k]=Inf32
            else
                isfinite(v)&&v>=0&&isfinite(Float32(v))||error("TIME grid contains an invalid solid-node travel time; check byte_order and grid coverage")
                result[i,j,k]=v
            end
        end
    end
    result
end

function prepare_travel_time_3d(cfg,stations,catalog,state;force_build=false)
    validate_grid3d_config(cfg);o=cfg["grid3d"];format=o["model_format"]
    phases=UInt8[p=="P" ? 1 : 2 for p in o["phases"]]
    headers=Dict{Tuple{UInt8,Int},NLLHeader3D}();inputs=String[]
    for phase in phases
        for j in (format=="nll_velocity" ? (0:0) : (1:length(stations)))
            label=phase==1 ? "P" : "S"
            if j>0
                occursin(r"^[A-Za-z0-9_.-]+$",stations.id[j])||error("Station code cannot be used in an NLL filename: $(stations.id[j])")
            end
            path=format=="nll_velocity" ? o[phase==1 ? "vp_file" : "vs_file"] : o["time_root"]*".$label.$(stations.id[j]).time.hdr"
            isempty(path)&&error("Missing 3D $label model path")
            h=read_nll_header(path;coordinate_system=o["coordinate_system"]);headers[(phase,j)]=h
            push!(inputs,h.file,nll_buffer_path(h))
        end
    end
    h=first(values(headers));all(t->same_nll_grid(h,t),values(headers))||error("All 3D P/S and station grids must have identical dimensions, origin, spacing and projection")
    lo,hi,dims,nodes=grid3d_plan(h,o)
    surface_file=o["surface_file"];surface_bytes=isempty(surface_file) ? 0 : filesize(surface_file)
    resources=grid3d_memory_guard(nodes,length(phases)*length(stations),maximum(nll_bytes(t) for t in values(headers)),surface_bytes,o)
    g=FMM3D.Grid3D((collect(range(lo[a],hi[a];length=dims[a])) for a in 1:3)...)
    surface=read_surface3d(surface_file,g)
    project(x,y)=projected_xy_jacobian(h.projection,x,y,state.ref_lat,state.ref_lon,state.ref_radius_m)
    sources=[(project(stations.x_m[j],stations.y_m[j])[1:2]...,stations.z_m[j]) for j in 1:length(stations)]
    for (j,xyz) in enumerate(sources)
        FMM3D._interpolation(g,xyz)
        surface!==nothing&&xyz[3]<FMM3D.surface_depth(g,surface,xyz[1],xyz[2])-1e-7&&error("Station $(stations.id[j]) lies above the surface; check elevations and datum")
        if format=="nll_time"
            for phase in phases
                t=headers[(phase,j)];t.kind=="TIME"||error("Precomputed grids must have type TIME")
                t.station[1]==stations.id[j]||error("TIME station label disagrees with stations.txt")
                maximum(abs.(t.station[2].-xyz))<=0.01||error("TIME station coordinates disagree with stations.txt by more than 1 cm; check projection and elevation datum")
            end
        end
    end
    isempty(surface_file)||push!(inputs,surface_file)
    identity=Dict{String,Any}("format_version"=>1,"byte_order_host"=>string(ENDIAN_BOM),"options"=>o,
        "inputs"=>Dict(abspath(p)=>file_hash3d(p) for p in sort!(unique(inputs))),
        "implementation"=>Dict(p=>file_hash3d(joinpath(@__DIR__,p)) for p in ("TravelTimes3D.jl","NLLGrids.jl","forward3d/FMM3D.jl","forward3d/Marching.jl")),
        "station_names"=>stations.id,"sources"=>[collect(s) for s in sources])
    identity_text=sprint(io->TOML.print(io,identity;sorted=true));key=bytes2hex(sha256(identity_text))
    cache=joinpath(o["cache_dir"],key);manifest=joinpath(cache,"complete.toml")
    records=isfile(manifest)&&!force_build ? get(TOML.parsefile(manifest),"fields",Dict()) : Dict()
    expected=["$(Int(phase))_$(j).bin" for phase in phases for j in 1:length(stations)]
    ready=all(file->haskey(records,file)&&isfile(joinpath(cache,file))&&filesize(joinpath(cache,file))==4nodes&&file_hash3d(joinpath(cache,file))==records[file],expected)
    if !ready
        Bool(cfgget(cfg,"travel_time","auto_build";default=true))||force_build||error("No complete compatible 3D cache and travel_time.auto_build=false")
        isfile(manifest)&&!Bool(cfgget(cfg,"travel_time","auto_rebuild";default=true))&&!force_build&&error("3D cache is damaged and travel_time.auto_rebuild=false")
        mkpath(cache);isfile(manifest)&&rm(manifest);records=Dict{String,String}()
        open(io->write(io,identity_text),joinpath(cache,"inputs.toml"),"w")
        for phase in phases
            velocity=format=="nll_velocity" ? sample_nll_velocity(headers[(phase,0)],g,o["model_interpolation"],o["byte_order"]) : nothing
            for j in 1:length(stations)
                @printf("  %s %s: %s\n",stations.id[j],phase==1 ? "P" : "S",format=="nll_velocity" ? "building FMM field" : "reading TIME grid")
                times=format=="nll_velocity" ? FMM3D.march(g,velocity,sources[j];accuracy_order=Int(o["accuracy_order"]),surface=surface) :
                    read_time_field3d(headers[(phase,j)],g,surface,o["byte_order"])
                field=Float32.(times)
                # Internal caches are native endian and record the host byte order in their key.
                file="$(Int(phase))_$(j).bin";dest=joinpath(cache,file)
                open(io->write(io,field),dest*".tmp","w");mv(dest*".tmp",dest;force=true)
                records[file]=file_hash3d(dest)
                times=nothing;field=nothing;GC.gc()
            end
            velocity=nothing;GC.gc()
        end
        open(io->TOML.print(io,Dict("fields"=>records);sorted=true),manifest*".tmp","w")
        mv(manifest*".tmp",manifest;force=true)
    else
        println("Using verified 3D travel-time cache: $cache")
    end
    fields=Dict{Tuple{UInt8,Float64,Float64,Float64},Array{Float32,3}}();valid=trues(dims)
    for phase in phases,j in 1:length(stations)
        field=load_cached_field3d(joinpath(cache,"$(Int(phase))_$(j).bin"),dims)
        fields[(phase,stations.x_m[j],stations.y_m[j],stations.z_m[j])]=field
        valid .&=isfinite.(field)
    end
    metadata=merge(resources,Dict{String,Any}("type"=>"3d","geometry"=>"cartesian","model_format"=>format,
        "cache_dir"=>cache,"cache_key"=>key,"dimensions"=>collect(dims),"spacing_m"=>collect(FMM3D.grid_spacing(g)),
        "topography_mask"=>surface!==nothing,"phases"=>o["phases"],"accuracy_order"=>Int(o["accuracy_order"])))
    model=TravelTime3D(:cartesian,g,h.projection,fields,valid,state.ref_lat,state.ref_lon,state.ref_radius_m,metadata)
    for i in 1:length(state)
        valid_location3d(model,state.x[i],state.y[i],state.z[i])||error("EventID $(state.event_id[i]) is outside the usable 3D interpolation domain. Enlarge the model or refine near the surface; no event is clamped or silently excluded.")
    end
    model
end

function travel_time_gradient(m::TravelTime3D,phase::UInt8,x,y,z,sx,sy,sz)
    key=(phase,Float64(sx),Float64(sy),Float64(sz))
    haskey(m.fields,key)||error("No 3D travel-time field for this station–phase; check grid3d.phases and stations.txt")
    X,Y,a,b,c,d=projected_xy_jacobian(m.projection,x,y,m.ref_lat,m.ref_lon,m.coordinate_radius_m)
    t,gx,gy,gz=FMM3D.value_gradient(m.grid,m.fields[key],(X,Y,z))
    t,a*gx+c*gy,b*gx+d*gy,gz
end

function valid_location3d(m::TravelTime3D,x,y,z)
    X,Y=projected_xy_jacobian(m.projection,x,y,m.ref_lat,m.ref_lon,m.coordinate_radius_m)[1:2]
    try
        ids,_,_=FMM3D._interpolation(m.grid,(X,Y,z))
        all(m.valid[i] for i in ids)
    catch e
        e isa DomainError||rethrow();false
    end
end

"Shorten a full update, preserving its direction, if it would leave the valid 3D domain."
function restrict_grid3d_step!(update,state,model::TravelTime3D)
    n=length(state);factor=1.
    for attempt in 0:20
        feasible=all(valid_location3d(model,state.x[i]+factor*update[i],state.y[i]+factor*update[n+i],state.z[i]+factor*update[2n+i]) for i in 1:n)
        if feasible
            factor<1&&(update .*= factor)
            factor<1&&@warn "Shortened a relocation update at the 3D grid/terrain boundary; inspect domain coverage before interpreting this solution" factor=factor
            return factor
        end
        factor/=2
    end
    error("Relocation cannot take a feasible step inside the 3D grid. Enlarge the model or refine the terrain grid; no boundary-clamped result is returned.")
end
