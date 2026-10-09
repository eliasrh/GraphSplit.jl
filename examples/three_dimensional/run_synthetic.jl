#!/usr/bin/env julia
# Reproducible small test; no packages beyond the Julia standard library.
using Random, Printf, TOML, LinearAlgebra
include(joinpath(@__DIR__,"..","..","src","GraphSplit.jl"))
using .GraphSplit
const F = GraphSplit.FMM3D

terrain_elevation(x,y) = 220.0 + 180.0*exp(-((x-1500)^2+(y+500)^2)/5000^2) + .012x
background(z) = 5000.0 + .12z
velocity(x,y,z) = background(z)*(1 + .16*exp(-((x-1800)^2+(y+1200)^2)/2600^2-((z-2500)/4000)^2) - .12*exp(-((x+2300)^2+(y-1200)^2)/2300^2-((z-3200)/4000)^2))

function write_catalog(file,xyz,lat0,lon0)
    open(file,"w") do io
        for i in axes(xyz,1)
            lat,lon=GraphSplit.local_xy_to_ll(xyz[i,1],xyz[i,2],lat0,lon0,6371000.)
            @printf(io,"2026 1 1 0 0 0.000000 %.10f %.10f %.7f 1.0 %d\n",lat,lon,xyz[i,3]/1000,i)
        end
    end
end
function error_metrics(xyz,truth)
    err=xyz.-truth;h=hypot.(err[:,1],err[:,2]);z=abs.(err[:,3])
    Dict("horizontal_rmse_m"=>sqrt(sum(abs2,h)/length(h)),"depth_rmse_m"=>sqrt(sum(abs2,z)/length(z)),
        "spatial_rmse_m"=>sqrt(sum(abs2,err)/size(err,1)),"events"=>size(err,1))
end

function main(out=joinpath(@__DIR__,"generated"))
    out=abspath(out);mkpath(out);rng=MersenneTwister(104729);n=100;ns=8;noise=.008;lat0=64.;lon0=-20.
    truth=hcat(4000 .* rand(rng,n).-2000,4000 .* rand(rng,n).-2000,1800 .+3000 .*rand(rng,n))
    initial=truth.+randn(rng,n,3).*reshape([600.,600.,900.],1,3);initial[:,3]=max.(700.,initial[:,3]);initial[1,:]=truth[1,:]
    write_catalog(joinpath(out,"truth.txt"),truth,lat0,lon0);write_catalog(joinpath(out,"catalog.txt"),initial,lat0,lon0)
    sx=Float64[-4800,-4800,-2400,2400,4800,4800,2400,-2400]
    sy=Float64[-2400,2400,4800,4800,2400,-2400,-4800,-4800]
    # Nodes shared by both grids; a 0.1 m buried sensor avoids rounding above ground.
    sz=-terrain_elevation.(sx,sy) .+ 0.1
    open(joinpath(out,"stations.txt"),"w") do io
        for j in 1:ns
            lat,lon=GraphSplit.local_xy_to_ll(sx[j],sy[j],lat0,lon0,6371000.)
            @printf(io,"ST%02d %.10f %.10f %.8f\n",j,lat,lon,-sz[j])
        end
    end
    # Input is a nodal NLL model; its native 150 m spacing also defines the
    # independent truth calculation. Relocation uses a 300 m grid.
    fine=F.Grid3D(-6000.:150.:6000.,-6000.:150.:6000.,-600.:150.:7500.)
    vp=[velocity(x,y,z) for x in fine.x,y in fine.y,z in fine.z];vs=vp./sqrt(3.)
    GraphSplit.write_nll_velocity(joinpath(out,"vp"),fine,vp;local_coordinates=true)
    GraphSplit.write_nll_velocity(joinpath(out,"vs"),fine,vs;local_coordinates=true)
    open(joinpath(out,"surface.txt"),"w") do io
        for x in fine.x,y in fine.y;@printf(io,"%.8f %.8f %.8f\n",x/1000,y/1000,terrain_elevation(x,y));end
    end
    surface=[-terrain_elevation(x,y) for x in fine.x,y in fine.y]
    mkpath(joinpath(out,"theta"));mkpath(joinpath(out,"thetastd"))
    println("Generating noisy travel-time potentials on the independent 150 m grid")
    truth_started=time()
    for (p,vel) in (("P",vp),("S",vs)),j in 1:ns
        t=F.march(fine,vel,(sx[j],sy[j],sz[j]);surface=surface)
        observed=[F.value_gradient(fine,t,Tuple(truth[i,:]))[1]+noise*randn(rng) for i in 1:n];observed.-=observed[1]
        name=@sprintf("ST%02d",j)
        open(joinpath(out,"theta","theta_$(name)_$p.txt"),"w") do io
            for i in 1:n;@printf(io,"%d %.10f 1\n",i,observed[i]);end
        end
        open(joinpath(out,"thetastd","std_theta_$(name)_$p.txt"),"w") do io
            for i in 1:n;@printf(io,"%d %.10f 1 %d\n",i,noise,n-1);end
        end
        t=nothing;GC.gc()
    end
    truth_seconds=time()-truth_started
    open(joinpath(out,"vm.txt"),"w") do io
        for z in (-1000.,8000.);@printf(io,"%.3f %.8f %.8f\n",z/1000,background(z)/1000,background(z)/sqrt(3.)/1000);end
    end
    cfg=GraphSplit.default_config()
    cfg["coordinates"]["reference"]="manual";cfg["coordinates"]["reference_latitude"]=lat0;cfg["coordinates"]["reference_longitude"]=lon0
    cfg["gauge"]["mode"]="pin";cfg["gauge"]["pin_event_ids"]=[1];cfg["gauge"]["pin_reference_catalog"]="truth.txt"
    cfg["graph"]["neighbors"]=20;cfg["graph"]["maximum_degree"]=40;cfg["graph"]["mutual"]=false
    cfg["observations"]["minimum_theta_degree"]=0
    for stage in ("prelocation","relocation")
        cfg[stage]["max_outer_iterations"]=40;cfg[stage]["damping_lambda"]=1e-6
        cfg[stage]["max_event_step_m"]=300.;cfg[stage]["stop_stall_iterations"]=5
    end
    cfg["travel_time"]["geometry"]="cartesian";cfg["travel_time"]["clamp_to_grid"]=false
    cfg["lookup"]["horizontal_step_m"]=100.;cfg["lookup"]["depth_step_m"]=100.;cfg["lookup"]["station_depth_step_m"]=100.
    cfg["lookup"]["maximum_distance_km"]=15.;cfg["lookup"]["minimum_depth_km"]=-.6;cfg["lookup"]["maximum_depth_km"]=7.5
    cfg["lookup"]["depth_margin_km"]=0.
    cfg["grid3d"]["coordinate_system"]="local";cfg["grid3d"]["vp_file"]="vp.hdr";cfg["grid3d"]["vs_file"]="vs.hdr"
    cfg["grid3d"]["spacing_m"]=[300.,300.,300.];cfg["grid3d"]["model_interpolation"]="slowness_linear";cfg["grid3d"]["surface_file"]="surface.txt"
    results=Dict{String,Any}("experiment"=>Dict("seed"=>104729,"events"=>n,"stations"=>ns,"potential_noise_std_s"=>noise,
        "truth_spacing_m"=>150.,"relocation_3d_spacing_m"=>300.,"truth_generation_seconds"=>truth_seconds,
        "reference"=>"Event 1 spatially pinned to truth in both runs; origin-time updates have zero mean"),"starting"=>error_metrics(initial,truth))
    for (label,kind) in (("one_dimensional","lookup"),("three_dimensional","3d"))
        c=deepcopy(cfg);c["travel_time"]["type"]=kind;c["io"]["output_dir"]="output_"*label
        path=joinpath(out,label*".toml");open(io->TOML.print(io,c;sorted=true),path,"w")
        start=time();result=GraphSplit.run(GraphSplit.load_config(path));elapsed=time()-start
        xyz=hcat(result.dd_state.x,result.dd_state.y,result.dd_state.z);metrics=error_metrics(xyz,truth);metrics["elapsed_seconds_including_tables"]=elapsed
        metrics["dd_observations"]=length(result.dd_observations);metrics["outer_iterations"]=result.dd_stats.iterations
        results[label]=metrics
        open(joinpath(out,label*"_errors.csv"),"w") do io
            println(io,"event_id,truth_x_m,truth_y_m,truth_z_m,initial_x_m,initial_y_m,initial_z_m,relocated_x_m,relocated_y_m,relocated_z_m")
            for i in 1:n;println(io,join([i;truth[i,:];initial[i,:];xyz[i,:]],','));end
        end
        GC.gc()
    end
    open(io->TOML.print(io,results;sorted=true),joinpath(out,"metrics.toml"),"w")
    println("\nKnown-model synthetic results (RMSE, metres):")
    for name in ("starting","one_dimensional","three_dimensional")
        m=results[name];@printf("%-20s H %8.2f   Z %8.2f   3D %8.2f\n",name,m["horizontal_rmse_m"],m["depth_rmse_m"],m["spatial_rmse_m"])
    end
    results
end
if abspath(PROGRAM_FILE)==@__FILE__
    main(isempty(ARGS) ? joinpath(@__DIR__,"generated") : ARGS[1])
end
