"Indexed heap: each trial node has one slot, so its storage cannot exceed N nodes."
mutable struct TrialHeap
    nodes::Vector{Int32}
    position::Vector{Int32}
    count::Int
    times::Array{Float64,3}
end
TrialHeap(times) = TrialHeap(zeros(Int32,length(times)),zeros(Int32,length(times)),0,times)
less(h,a,b) = _time_less(h.times[a],a,h.times[b],b)
function lower!(h::TrialHeap, node::Int)
    k=Int(h.position[node])
    if k==0
        h.count+=1;k=h.count;h.nodes[k]=node
    end
    while k>1
        parent=k>>1;p=Int(h.nodes[parent])
        less(h,node,p)||break
        h.nodes[k]=p;h.position[p]=k;k=parent
    end
    h.nodes[k]=node;h.position[node]=k
end
function pop_trial!(h::TrialHeap)
    node=Int(h.nodes[1]);last=Int(h.nodes[h.count]);h.count-=1;h.position[node]=0
    if h.count>0
        k=1
        while 2k<=h.count
            c=2k
            c<h.count&&less(h,h.nodes[c+1],h.nodes[c])&&(c+=1)
            child=Int(h.nodes[c]);less(h,child,last)||break
            h.nodes[k]=child;h.position[child]=k;k=c
        end
        h.nodes[k]=last;h.position[last]=k
    end
    node
end

"Depth of a bilinear land surface, in the same positive-down datum as the grid."
function surface_depth(g::Grid3D, surface::Matrix{Float64}, x,y)
    i,u=_cell_fraction(g.x,x);j,v=_cell_fraction(g.y,y)
    (1-u)*(1-v)*surface[i,j]+u*(1-v)*surface[i+1,j]+(1-u)*v*surface[i,j+1]+u*v*surface[i+1,j+1]
end

"Check a short source-to-node segment within a single horizontal grid cell."
function seed_below_surface(g,surface,source,point)
    surface===nothing&&return true
    clearance(t)=begin
        x=(1-t)*source[1]+t*point[1];y=(1-t)*source[2]+t*point[2]
        z=(1-t)*source[3]+t*point[3]
        z-surface_depth(g,surface,x,y)
    end
    # Bilinear height on a straight segment is quadratic. Its minimum is
    # determined by the endpoints and, when convex, its interior vertex.
    c0=clearance(0.);c1=clearance(1.);cm=clearance(.5)
    a=2*(c1+c0-2cm);b=c1-c0-a
    min(c0,c1)>=-1e-7||return false
    a>0&&0 < -b/(2a) < 1 ? clearance(-b/(2a))>=-1e-7 : true
end

"Relocation-only Cartesian FMM; stores times, without model-derivative arrays."
function march(g::Grid3D, velocity::Array{Float64,3}, station;
        accuracy_order::Int=2, surface::Union{Nothing,Matrix{Float64}}=nothing)
    accuracy_order in (1,2)||error("FMM accuracy_order must be 1 or 2")
    length(g)<=typemax(Int32)||error("3D grid exceeds the 32-bit node-index limit")
    size(velocity)==size(g)||error("Velocity size does not match the 3D grid")
    all(v->isfinite(v)&&v>0,velocity)||error("FMM velocities must be positive and finite")
    source=ntuple(a->Float64(station[a]),3)
    _interpolation(g,source)
    source=ntuple(a->_canonical_source_coordinate((g.x,g.y,g.z)[a],source[a]),3)
    if surface!==nothing
        size(surface)==size(g)[1:2]||error("Surface size does not match the grid")
        source[3]>=surface_depth(g,surface,source[1],source[2])-1e-7||
            error("Station lies above the supplied topography. Correct its elevation or the surface; stations are not moved automatically.")
    end
    nx,ny,nz=size(g);n=length(g)
    active=surface===nothing ? trues(n) : vec(BitArray([z>=surface[i,j] for i in 1:nx,j in 1:ny,z in g.z]))
    slow=vec(1 ./ velocity)
    # Only boundary-cell quadrature uses the extension above the surface.
    # No air node is ever accepted by marching. Avoid mixing an arbitrary
    # supplied air velocity into a short receiver-to-rock seed segment.
    if surface!==nothing
        for j in 1:ny,i in 1:nx
            k=findfirst(z->z>=surface[i,j],g.z)
            k===nothing&&error("Topography leaves a grid column without solid nodes")
            for l in 1:k-1;slow[i+nx*(j-1)+nx*ny*(l-1)]=slow[i+nx*(j-1)+nx*ny*(k-1)];end
        end
    end
    times=fill(Inf,size(g));accepted=falses(n);heap=TrialHeap(times)
    iz,_=_cell_fraction(g.z,source[3])
    # One extra depth plane provides below-ground seeds for an off-node
    # receiver on sloping terrain. Every segment is checked against the DEM.
    zseeds=surface===nothing ? _source_axis_nodes(g.z,source[3]) : max(1,iz-1):min(nz,iz+2)
    for k in zseeds,j in _source_axis_nodes(g.y,source[2]),i in _source_axis_nodes(g.x,source[1])
        m=i+nx*(j-1)+nx*ny*(k-1);active[m]||continue
        point=(g.x[i],g.y[j],g.z[k]);seed_below_surface(g,surface,source,point)||continue
        # Extended vertical seeds can cross a horizontal cell boundary;
        # integrate each vertical-grid interval separately.
        cuts=Float64[0.,1.]
        if point[3]!=source[3]
            for z in g.z
                t=(z-source[3])/(point[3]-source[3]);0<t<1&&push!(cuts,t)
            end
        end
        sort!(cuts);value=0.
        for l in 1:length(cuts)-1
            a=ntuple(q->source[q]+cuts[l]*(point[q]-source[q]),3)
            b=ntuple(q->source[q]+cuts[l+1]*(point[q]-source[q]),3)
            value+=first(_seed_segment(g,slow,a,b))
        end
        times[m]=value;lower!(heap,m)
    end
    heap.count>0||error("No solid grid node can be connected to the station; refine the grid or check topography")
    spacing=grid_spacing(g);tag=Val(accuracy_order)
    while heap.count>0
        m=pop_trial!(heap);accepted[m]=true
        i=mod(m-1,nx)+1;j=mod(div(m-1,nx),ny)+1;k=div(m-1,nx*ny)+1
        neighbors=(i>1 ? m-1 : 0,i<nx ? m+1 : 0,j>1 ? m-nx : 0,
                   j<ny ? m+nx : 0,k>1 ? m-nx*ny : 0,k<nz ? m+nx*ny : 0)
        for q in neighbors
            (q==0||accepted[q]||!active[q])&&continue
            value=first(_upwind_update(q,slow,times,accepted,size(g),spacing,tag))
            value<times[q]||continue
            times[q]=value;lower!(heap,q)
        end
    end
    all(accepted[active])||error("The solid grid is disconnected from this station; enlarge the depth extent or inspect the surface")
    times
end

"Time and the exact derivative of its trilinear interpolant (seconds/metre)."
function value_gradient(g::Grid3D, times::AbstractArray{<:Real,3}, xyz)
    ids,ws,(u,v,w)=_interpolation(g,xyz);hx,hy,hz=grid_spacing(g)
    t=0.;gx=0.;gy=0.;gz=0.
    for n in 1:8
        a=(n-1)&1;b=((n-1)>>1)&1;c=((n-1)>>2)&1
        ax=a==0 ? 1-u : u;ay=b==0 ? 1-v : v;az=c==0 ? 1-w : w
        tn=Float64(times[ids[n]])
        isfinite(tn)||throw(DomainError(xyz,"Travel-time interpolation touches air or an unreachable node. Use a finer grid or a deeper event; no extrapolation is applied."))
        t+=ws[n]*tn;gx+=(a==0 ? -1 : 1)*ay*az*tn/hx
        gy+=(b==0 ? -1 : 1)*ax*az*tn/hy;gz+=(c==0 ? -1 : 1)*ax*ay*tn/hz
    end
    t,gx,gy,gz
end
