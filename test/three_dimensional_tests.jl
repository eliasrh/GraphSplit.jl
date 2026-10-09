const F3 = GraphSplit.FMM3D

function tiny3d_fixture(dir;surface=false)
    g=F3.Grid3D(-2000.:250.:2000.,-2000.:250.:2000.,-500.:250.:3000.)
    vp=fill(5000.,size(g));vs=vp./sqrt(3.)
    for (name,v) in (("vp",vp),("vs",vs))
        GraphSplit.write_nll_velocity(joinpath(dir,name),g,v;local_coordinates=true)
    end
    write(joinpath(dir,"catalog.txt"),"2026 1 1 0 0 0 0.002 0.003 1.2 1 1\n2026 1 1 0 0 0 -0.003 0.005 1.5 1 2\n")
    write(joinpath(dir,"stations.txt"),"A 0.001 0.002 0\nB -0.004 -0.005 100\n")
    cfg=GraphSplit.default_config();cfg["travel_time"]["type"]="3d";cfg["coordinates"]["reference"]="manual"
    cfg["grid3d"]["coordinate_system"]="local";cfg["grid3d"]["spacing_m"]=[250.,250.,250.]
    for name in ("vp","vs");cfg["grid3d"][name*"_file"]=joinpath(dir,name*".hdr");end
    cfg["grid3d"]["cache_dir"]=joinpath(dir,"cache")
    catalog=GraphSplit.read_catalog(joinpath(dir,"catalog.txt"),cfg);stations=GraphSplit.read_stations(joinpath(dir,"stations.txt"))
    state=GraphSplit.attach_coordinates!(stations,catalog,cfg)
    cfg,catalog,stations,state,g
end

@testset "Experimental fixed 3D model" begin
    @testset "FMM homogeneous accuracy and refinement" begin
        errors=Float64[]
        for spacing in (200.,100.,50.)
            g=F3.Grid3D(0:spacing:2000.,0:spacing:2000.,0:spacing:2000.)
            t=F3.march(g,fill(5000.,size(g)),(0.,0.,0.))
            points=((1131.,1377.,1621.),(731.,997.,1871.),(1789.,1457.,613.))
            push!(errors,maximum(abs(F3.value_gradient(g,t,p)[1]-norm(collect(p))/5000) for p in points))
            @test all(isfinite,t)
            @test t[1,1,1]==0
        end
        @test errors[3]<errors[2]<errors[1]
        @test errors[3]<.0015
        g=F3.Grid3D(0:100.:1000.,0:100.:1000.,0:100.:1000.)
        one=F3.march(g,fill(5000.,size(g)),(0.,0.,0.);accuracy_order=1)
        two=F3.march(g,fill(5000.,size(g)),(0.,0.,0.);accuracy_order=2)
        @test abs(two[end,end,end]-sqrt(3.)/5)<abs(one[end,end,end]-sqrt(3.)/5)
        @test_throws DomainError F3.value_gradient(g,two,(1001.,500.,500.))
    end

    @testset "Station elevation, sloping terrain and air exclusion" begin
        g=F3.Grid3D(-1000.:100.:1000.,-1000.:100.:1000.,-500.:100.:2000.)
        surface=[-200.0 + 0.1*x + 0.05*y for x in g.x, y in g.y]
        source=(33.,-71.,-200.0 + 0.1*33 - 0.05*71)
        t=F3.march(g,fill(5000.,size(g)),source;surface=surface)
        @test all(isfinite(t[i,j,k])==(g.z[k]>=surface[i,j]) for k in eachindex(g.z),j in eachindex(g.y),i in eachindex(g.x))
        @test_throws ErrorException F3.march(g,fill(5000.,size(g)),(source[1],source[2],source[3]-1.);surface=surface)
        @test_throws DomainError F3.value_gradient(g,t,(0.,0.,-300.))
        p=(313.,271.,803.);analytic=norm(collect(p).-collect(source))/5000
        @test abs(F3.value_gradient(g,t,p)[1]-analytic)<.007
        @test !F3.seed_below_surface(g,surface,source,(source[1],source[2],source[3]-100.))
    end

    @testset "NLL data order, byte order, velocity units and projection" begin
        mktempdir() do dir
            g=F3.Grid3D(0:100.:200.,0:200.:600.,-100.:100.:300.)
            v=[4000.0 + 100*i + 10*j + k for i in 1:3,j in 1:4,k in 1:5]
            path=GraphSplit.write_nll_velocity(joinpath(dir,"vp"),g,v;latitude=64.,longitude=-20.,rotation=27.)
            h=GraphSplit.read_nll_header(path)
            @test GraphSplit.sample_nll_velocity(h,g,"nearest","little")≈v rtol=1e-7
            # Re-encode independently in big-endian SLOWNESS DOUBLE.
            write(path,replace(read(path,String),"VELOCITY FLOAT"=>"SLOWNESS DOUBLE"))
            open(joinpath(dir,"vp.buf"),"w") do io
                for i in 1:3,j in 1:4,k in 1:5;write(io,hton(reinterpret(UInt64,1000/v[i,j,k])));end
            end
            h=GraphSplit.read_nll_header(path)
            @test GraphSplit.sample_nll_velocity(h,g,"slowness_linear","big")≈v
            @test_throws ErrorException GraphSplit.read_nll_header(path;coordinate_system="local")
            p=h.projection;x=10013.;y=17801.;args=(64.,-20.,6371000.)
            xy=GraphSplit.projected_xy_jacobian(p,x,y,args...)
            e=.02
            for (axis,derivative) in ((1,(xy[3],xy[5])),(2,(xy[4],xy[6])))
                a=GraphSplit.projected_xy_jacobian(p,x+(axis==1 ? e : 0),y+(axis==2 ? e : 0),args...)
                b=GraphSplit.projected_xy_jacobian(p,x-(axis==1 ? e : 0),y-(axis==2 ? e : 0),args...)
                @test collect(derivative)≈(collect(a[1:2])-collect(b[1:2]))/(2e) rtol=1e-7 atol=1e-9
            end
            # Independent NLL SIMPLE formula in geographic coordinates, including rotation.
            lat,lon=GraphSplit.local_xy_to_ll(x,y,args...);xx=6371000*cosd(lat)*deg2rad(lon+20);yy=6371000*deg2rad(lat-64)
            @test xy[1]≈cosd(27)*xx+sind(27)*yy
            @test xy[2]≈-sind(27)*xx+cosd(27)*yy
            write(path,replace(read(path,String),"SIMPLE"=>"LAMBERT"))
            @test_throws ErrorException GraphSplit.read_nll_header(path)
        end
    end

    @testset "Provider, cache, spatial derivatives and domain guard" begin
        mktempdir() do dir
            cfg,cat,sta,state,g=tiny3d_fixture(dir)
            @test GraphSplit.validate_config(cfg)===cfg
            limited=deepcopy(cfg);limited["grid3d"]["maximum_memory_gib"]=1e-6;limited["grid3d"]["warning_memory_gib"]=1e-7
            @test_throws ErrorException GraphSplit.prepare_travel_time(limited,sta,cat,state)
            @test !isdir(cfg["grid3d"]["cache_dir"])
            m=GraphSplit.prepare_travel_time(cfg,sta,cat,state)
            @test length(m.fields)==4
            q=[333.,277.,1333.];source=(sta.x_m[1],sta.y_m[1],sta.z_m[1])
            v=GraphSplit.travel_time_gradient(m,UInt8(1),q...,source...)
            for a in 1:3
                plus=copy(q);minus=copy(q);plus[a]+=.1;minus[a]-=.1
                derivative=(GraphSplit.travel_time_gradient(m,UInt8(1),plus...,source...)[1]-GraphSplit.travel_time_gradient(m,UInt8(1),minus...,source...)[1])/.2
                @test derivative≈v[a+1] atol=1e-12
            end
            reread=GraphSplit.prepare_travel_time(cfg,sta,cat,state)
            @test reread.metadata["cache_key"]==m.metadata["cache_key"]
            @test GraphSplit.travel_time_gradient(reread,UInt8(1),q...,source...)==v
            update=zeros(8);update[1]=4000.
            factor=GraphSplit.restrict_grid3d_step!(update,state,m)
            @test 0<factor<1
            @test GraphSplit.valid_location3d(m,state.x[1]+update[1],state.y[1],state.z[1])
            @test !GraphSplit.valid_location3d(m,9999.,0.,1000.)
            # A changed model must not reuse old times.
            GraphSplit.write_nll_velocity(joinpath(dir,"vp"),g,fill(6000.,size(g));local_coordinates=true)
            changed=GraphSplit.prepare_travel_time(cfg,sta,cat,state)
            @test changed.metadata["cache_key"]!=m.metadata["cache_key"]
            @test GraphSplit.travel_time_gradient(changed,UInt8(1),q...,source...)[1]<v[1]
            # A truncated field is rebuilt, or explicitly refused when rebuilding is disabled.
            file=joinpath(changed.metadata["cache_dir"],"1_1.bin");write(file,"broken")
            refuse=deepcopy(cfg);refuse["travel_time"]["auto_rebuild"]=false
            @test_throws ErrorException GraphSplit.prepare_travel_time(refuse,sta,cat,state)
            repaired=GraphSplit.prepare_travel_time(cfg,sta,cat,state)
            @test isfinite(GraphSplit.travel_time_gradient(repaired,UInt8(1),q...,source...)[1])
        end
    end
    @testset "Layer refraction and convergence" begin
        errors=Float64[]
        for step in (100.,50.)
            g=F3.Grid3D(0:step:8000.,0:step:500.,0:step:3000.)
            v=[z<1000. ? 4000. : 6000. for x in g.x,y in g.y,z in g.z]
            t=F3.march(g,v,(0.,0.,0.))
            expected=6000/6000+2000*sqrt(1/4000^2-1/6000^2)
            push!(errors,abs(F3.value_gradient(g,t,(6000.,0.,0.))[1]-expected))
        end
        @test errors[2]<errors[1]
        @test errors[2]<.02
    end

    @testset "Precomputed NLL TIME round trip and source checks" begin
        mktempdir() do dir
            cfg,cat,sta,state,g=tiny3d_fixture(dir)
            cfg["grid3d"]["phases"]=["P"]
            m=GraphSplit.prepare_travel_time(cfg,sta,cat,state)
            root=joinpath(dir,"precomputed")
            for j in 1:length(sta)
                stem=root*".P.$(sta.id[j]).time"
                field=m.fields[(UInt8(1),sta.x_m[j],sta.y_m[j],sta.z_m[j])]
                open(stem*".hdr","w") do io
                    println(io,join(size(g),' ')," ",join([g.x[1],g.y[1],g.z[1]]./1000,' ')," ",join(collect(F3.grid_spacing(g))./1000,' ')," TIME FLOAT")
                    println(io,sta.id[j]," ",join([sta.x_m[j],sta.y_m[j],sta.z_m[j]]./1000,' '))
                    println(io,"TRANSFORM NONE")
                end
                open(stem*".buf","w") do io
                    for i in eachindex(g.x),k in eachindex(g.y),l in eachindex(g.z)
                        write(io,htol(reinterpret(UInt32,field[i,k,l])))
                    end
                end
            end
            cfg["grid3d"]["model_format"]="nll_time";cfg["grid3d"]["time_root"]=root
            cfg["grid3d"]["spacing_m"]=[1000.,1000.,1000.] # imported TIME spacing is preserved
            imported=GraphSplit.prepare_travel_time(cfg,sta,cat,state)
            @test size(imported.grid)==size(m.grid)
            @test imported.fields==m.fields
            q=(333.,277.,1333.);source=(sta.x_m[1],sta.y_m[1],sta.z_m[1])
            @test GraphSplit.travel_time_gradient(imported,UInt8(1),q...,source...)==GraphSplit.travel_time_gradient(m,UInt8(1),q...,source...)
            sta.z_m[1]+=1.
            @test_throws ErrorException GraphSplit.prepare_travel_time(cfg,sta,cat,state)
        end
    end

    @testset "Small complete 3D relocation writes outputs" begin
        mktempdir() do dir
            cfg,cat,sta,state,g=tiny3d_fixture(dir)
            m=GraphSplit.prepare_travel_time(cfg,sta,cat,state)
            mkpath(joinpath(dir,"theta"))
            for phase in (UInt8(1),UInt8(2)),j in 1:length(sta)
                times=[GraphSplit.travel_time_gradient(m,phase,state.x[i],state.y[i],state.z[i],sta.x_m[j],sta.y_m[j],sta.z_m[j])[1] for i in 1:2]
                label=phase==1 ? "P" : "S"
                write(joinpath(dir,"theta","theta_$(sta.id[j])_$label.txt"),"1 0.0 1\n2 $(times[2]-times[1]) 1\n")
            end
            for (key,value) in (("catalog_file","catalog.txt"),("stations_file","stations.txt"),("theta_dir","theta"),("thetastd_dir","thetastd"),("output_dir","output"))
                cfg["io"][key]=joinpath(dir,value)
            end
            cfg["observations"]["minimum_observations_per_pair"]=1;cfg["observations"]["minimum_theta_degree"]=0
            cfg["graph"]["neighbors"]=1;cfg["graph"]["maximum_degree"]=1
            cfg["gauge"]["mode"]="pin";cfg["gauge"]["pin_event_ids"]=[1]
            for stage in ("prelocation","relocation");cfg[stage]["max_outer_iterations"]=5;cfg[stage]["verbose"]=false;end
            result=GraphSplit.run(cfg)
            @test maximum(abs.(result.dd_state.z.-state.z))<1e-5
            @test result.dd_stats.converged
            @test result.dd_stats.rms_s[end]<1e-10
            @test TOML.parsefile(joinpath(dir,"output","run_summary.toml"))["travel_time"]["type"]=="3d"
            @test isfile(joinpath(dir,"output","catalog_dd.txt"))
        end
    end

end
