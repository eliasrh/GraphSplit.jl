struct KDTree
    points::Matrix{Float64}   # dimensions x points
    point::Vector{Int32}
    left::Vector{Int32}
    right::Vector{Int32}
    axis::Vector{UInt8}
    root::Int32
end

@inline function kd_less(points::Matrix{Float64}, axis::Int, first::Int, second::Int)
    a, b = points[axis, first], points[axis, second]
    return a < b || (a == b && first < second)
end

@inline function swap_indices!(values::Vector{Int}, first::Int, second::Int)
    values[first], values[second] = values[second], values[first]
end

function kd_partition!(permutation::Vector{Int}, lo::Int, hi::Int, pivot_position::Int,
        points::Matrix{Float64}, axis::Int)
    pivot_index = permutation[pivot_position]
    swap_indices!(permutation, pivot_position, hi)
    store = lo
    for position in lo:(hi - 1)
        if kd_less(points, axis, permutation[position], pivot_index)
            swap_indices!(permutation, store, position)
            store += 1
        end
    end
    swap_indices!(permutation, store, hi)
    return store
end

"Select one axis median in place without sorting an entire recursive subtree."
function kd_select!(permutation::Vector{Int}, lo::Int, hi::Int, target::Int,
        points::Matrix{Float64}, axis::Int)
    while lo < hi
        pivot = (lo + hi) >>> 1
        pivot = kd_partition!(permutation, lo, hi, pivot, points, axis)
        pivot == target && return
        if target < pivot
            hi = pivot - 1
        else
            lo = pivot + 1
        end
    end
end

function KDTree(points_by_row::Matrix{Float64})
    n, dimensions = size(points_by_row)
    n > 0 || error("Cannot build a k-d tree with no points")
    points = permutedims(points_by_row)
    permutation = collect(1:n)
    node_point, left, right, axes = Int32[], Int32[], Int32[], UInt8[]
    function build!(lo::Int, hi::Int, depth::Int)::Int32
        lo > hi && return Int32(0)
        axis = mod(depth, dimensions) + 1
        middle = (lo + hi) >>> 1
        kd_select!(permutation, lo, hi, middle, points, axis)
        push!(node_point, Int32(permutation[middle]))
        push!(left, Int32(0)); push!(right, Int32(0)); push!(axes, UInt8(axis))
        node = Int32(length(node_point))
        left[node] = build!(lo, middle - 1, depth + 1)
        right[node] = build!(middle + 1, hi, depth + 1)
        return node
    end
    root = build!(1, n, 0)
    return KDTree(points, node_point, left, right, axes, root)
end

mutable struct NeighborHeap
    index::Vector{Int32}
    distance2::Vector{Float64}
    count::Int
end

NeighborHeap(k::Int) = NeighborHeap(fill(Int32(0), k), fill(Inf, k), 0)

function heap_swap!(heap::NeighborHeap, a::Int, b::Int)
    heap.index[a], heap.index[b] = heap.index[b], heap.index[a]
    heap.distance2[a], heap.distance2[b] = heap.distance2[b], heap.distance2[a]
end

function heap_sift_up!(heap::NeighborHeap, position::Int)
    while position > 1
        parent = position >>> 1
        heap.distance2[parent] >= heap.distance2[position] && break
        heap_swap!(heap, parent, position)
        position = parent
    end
end

function heap_sift_down!(heap::NeighborHeap, position::Int)
    while true
        child = position << 1
        child > heap.count && break
        if child < heap.count && heap.distance2[child + 1] > heap.distance2[child]
            child += 1
        end
        heap.distance2[position] >= heap.distance2[child] && break
        heap_swap!(heap, position, child)
        position = child
    end
end

function heap_consider!(heap::NeighborHeap, index::Int32, distance2::Float64)
    capacity = length(heap.index)
    if heap.count < capacity
        heap.count += 1
        heap.index[heap.count] = index
        heap.distance2[heap.count] = distance2
        heap_sift_up!(heap, heap.count)
    elseif distance2 < heap.distance2[1]
        heap.index[1] = index
        heap.distance2[1] = distance2
        heap_sift_down!(heap, 1)
    end
end

function knn(tree::KDTree, query_index::Int, k::Int, radius_m::Float64=Inf)
    k = min(k, size(tree.points, 2) - 1)
    k <= 0 && return Int32[], Float64[]
    heap = NeighborHeap(k)
    radius2 = radius_m * radius_m
    function visit(node::Int32)
        node == 0 && return
        point_index = tree.point[node]
        axis = Int(tree.axis[node])
        difference = tree.points[axis, query_index] - tree.points[axis, point_index]
        near = difference <= 0.0 ? tree.left[node] : tree.right[node]
        far = difference <= 0.0 ? tree.right[node] : tree.left[node]
        visit(near)
        if point_index != query_index
            distance2 = 0.0
            @inbounds for dim in axes(tree.points, 1)
                d = tree.points[dim, query_index] - tree.points[dim, point_index]
                distance2 += d * d
            end
            distance2 <= radius2 && heap_consider!(heap, point_index, distance2)
        end
        threshold = heap.count < k ? radius2 : min(radius2, heap.distance2[1])
        difference * difference <= threshold && visit(far)
    end
    visit(tree.root)
    order = sortperm(view(heap.distance2, 1:heap.count))
    return heap.index[order], sqrt.(heap.distance2[order])
end

edge_key(i::Integer, j::Integer) = (UInt64(min(i, j)) << 32) | UInt64(max(i, j))
directed_key(i::Integer, j::Integer) = (UInt64(i) << 32) | UInt64(j)

function metric_points(state::State, cfg::AbstractDict)
    metric = lowercase(String(cfgget(cfg, "graph", "metric"; default="xyz_scaled")))
    if metric == "xy"
        return hcat(state.x, state.y)
    elseif metric == "xyz"
        return hcat(state.x, state.y, state.z)
    end
    scale = Float64(cfgget(cfg, "graph", "depth_scale"; default=0.5))
    return hcat(state.x, state.y, scale .* state.z)
end

function directed_neighbors(points::Matrix{Float64}, k::Int, radius::Float64)
    tree = KDTree(points)
    source, target, distance = Int32[], Int32[], Float64[]
    for i in axes(points, 1)
        neighbors, distances = knn(tree, i, k, radius)
        append!(source, fill(Int32(i), length(neighbors)))
        append!(target, neighbors)
        append!(distance, distances)
    end
    return source, target, distance
end

function undirected_candidates(points::Matrix{Float64}, k::Int, radius::Float64, mutual::Bool)
    source, target, distance = directed_neighbors(points, k, radius)
    directed = mutual ? Set(directed_key(source[e], target[e]) for e in eachindex(source)) : Set{UInt64}()
    by_key = Dict{UInt64,Tuple{Int32,Int32,Float64}}()
    for e in eachindex(source)
        i, j = source[e], target[e]
        mutual && !(directed_key(j, i) in directed) && continue
        key = edge_key(i, j)
        if !haskey(by_key, key) || distance[e] < by_key[key][3]
            by_key[key] = (min(i, j), max(i, j), distance[e])
        end
    end
    candidates = collect(values(by_key))
    sort!(candidates; by=edge -> edge[3])
    return candidates
end

mutable struct UnionFind
    parent::Vector{Int32}
    size::Vector{Int32}
end
UnionFind(n::Int) = UnionFind(Int32.(1:n), ones(Int32, n))

function uf_find!(uf::UnionFind, x::Int32)
    root = x
    while uf.parent[root] != root
        root = uf.parent[root]
    end
    while uf.parent[x] != x
        next = uf.parent[x]
        uf.parent[x] = root
        x = next
    end
    return root
end

function uf_union!(uf::UnionFind, a::Int32, b::Int32)
    ra, rb = uf_find!(uf, a), uf_find!(uf, b)
    ra == rb && return false
    if uf.size[ra] < uf.size[rb]
        ra, rb = rb, ra
    end
    uf.parent[rb] = ra
    uf.size[ra] += uf.size[rb]
    return true
end

function graph_components(n::Int, edges_i::Vector{Int32}, edges_j::Vector{Int32})
    uf = UnionFind(n)
    for e in eachindex(edges_i)
        uf_union!(uf, edges_i[e], edges_j[e])
    end
    roots = Dict{Int32,Int32}()
    labels = Vector{Int32}(undef, n)
    sizes = Dict{Int32,Int32}()
    for i in 1:n
        root = uf_find!(uf, Int32(i))
        label = get!(roots, root, Int32(length(roots) + 1))
        labels[i] = label
        sizes[label] = get(sizes, label, Int32(0)) + 1
    end
    event_size = Int32[get(sizes, label, Int32(1)) for label in labels]
    return labels, event_size, length(roots)
end

function approximate_fiedler(n::Int, edges_i, edges_j, comp_id, degree, iterations::Int)
    rng = MersenneTwister(0x4753504c)
    x = randn(rng, n)
    ncomp = maximum(comp_id)
    function project_components!(v)
        sums = zeros(ncomp); counts = zeros(Int, ncomp)
        for i in 1:n
            c = comp_id[i]; sums[c] += v[i]; counts[c] += 1
        end
        for i in 1:n
            v[i] -= sums[comp_id[i]] / counts[comp_id[i]]
        end
    end
    project_components!(x)
    normx = norm(x); normx > 0.0 && (x ./= normx)
    alpha = 0.9 / max(2.0, 2.0 * maximum(degree))
    lap = zeros(n)
    for _ in 1:iterations
        fill!(lap, 0.0)
        for e in eachindex(edges_i)
            i, j = edges_i[e], edges_j[e]
            difference = x[i] - x[j]
            lap[i] += difference; lap[j] -= difference
        end
        x .-= alpha .* lap
        project_components!(x)
        normx = norm(x)
        normx > 0.0 && (x ./= normx)
    end
    return x
end

function edge_distances(state::State, points::Matrix{Float64}, i::Int32, j::Int32)
    dmetric2 = 0.0
    for dim in axes(points, 2)
        d = points[i, dim] - points[j, dim]
        dmetric2 += d * d
    end
    dxy = hypot(state.x[i] - state.x[j], state.y[i] - state.y[j])
    dxyz = hypot(dxy, state.z[i] - state.z[j])
    return sqrt(dmetric2), dxy, dxyz
end

function quantile_sorted(values::Vector{Float64}, probability::Float64)
    isempty(values) && return NaN
    ordered = sort(values)
    position = clamp(1.0 + probability * (length(ordered) - 1), 1.0, length(ordered))
    lo, hi = floor(Int, position), ceil(Int, position)
    return lo == hi ? ordered[lo] : ordered[lo] + (position - lo) * (ordered[hi] - ordered[lo])
end

"Build GraphSplit's sparse local graph with exact dependency-free kNN."
function build_event_graph(state::State, cfg::AbstractDict)
    n = length(state)
    n <= 1 && return EventGraph(Int32[], Int32[], Float64[], Float64[], Float64[], zeros(Int32, n),
        Int32.(1:n), ones(Int32, n), n, falses(0))
    points = metric_points(state, cfg)
    k = min(Int(cfgget(cfg, "graph", "neighbors"; default=20)), n - 1)
    maximum_degree_cfg = Int(cfgget(cfg, "graph", "maximum_degree"; default=30))
    maximum_degree = maximum_degree_cfg <= 0 ? typemax(Int32) : maximum_degree_cfg
    max_distance_km = Float64(cfgget(cfg, "graph", "maximum_distance_km"; default=0.0))
    radius = max_distance_km > 0.0 ? 1000.0 * max_distance_km : Inf
    mutual = Bool(cfgget(cfg, "graph", "mutual"; default=true))
    candidates = undirected_candidates(points, k, radius, mutual)
    isempty(candidates) && error("Event graph has no edges; relax mutual/radius graph settings")

    selected = Tuple{Int32,Int32,Float64}[]
    selected_keys = Set{UInt64}()
    degree = zeros(Int32, n)
    for edge in candidates
        i, j, _ = edge
        if degree[i] < maximum_degree && degree[j] < maximum_degree
            push!(selected, edge); push!(selected_keys, edge_key(i, j))
            degree[i] += 1; degree[j] += 1
        end
    end
    minimum_degree = Int(cfgget(cfg, "graph", "minimum_degree"; default=0))
    if minimum_degree > 0
        for edge in candidates
            i, j, _ = edge
            edge_key(i, j) in selected_keys && continue
            (degree[i] < minimum_degree || degree[j] < minimum_degree) || continue
            degree[i] < maximum_degree && degree[j] < maximum_degree || continue
            push!(selected, edge); push!(selected_keys, edge_key(i, j))
            degree[i] += 1; degree[j] += 1
        end
    end

    if Bool(cfgget(cfg, "graph", "ensure_connected"; default=false))
        uf = UnionFind(n)
        for (i, j, _) in selected
            uf_union!(uf, i, j)
        end
        bridge_cfg = Float64(cfgget(cfg, "graph", "maximum_bridge_distance_km"; default=0.0))
        bridge_max = bridge_cfg > 0.0 ? 1000.0 * bridge_cfg : Inf
        for edge in candidates
            edge[3] <= bridge_max || continue
            if uf_union!(uf, edge[1], edge[2])
                edge_key(edge[1], edge[2]) in selected_keys || begin
                    push!(selected, edge); push!(selected_keys, edge_key(edge[1], edge[2]))
                    degree[edge[1]] += 1; degree[edge[2]] += 1
                end
            end
        end
    end

    augmented = falses(length(selected))
    aug_cfg = cfgget(cfg, "graph", "augmentation"; default=Dict{String,Any}())
    if Bool(get(aug_cfg, "enabled", false))
        ci = Int(get(aug_cfg, "candidate_neighbors", max(2k, 40)))
        broad = undirected_candidates(points, min(ci, n - 1), radius, mutual)
        base_i = Int32[edge[1] for edge in selected]; base_j = Int32[edge[2] for edge in selected]
        comp_id, _, _ = graph_components(n, base_i, base_j)
        fiedler = approximate_fiedler(n, base_i, base_j, comp_id, degree, Int(get(aug_cfg, "power_iterations", 30)))
        add_count = Int(get(aug_cfg, "add_edges", 0))
        add_per_event = Float64(get(aug_cfg, "add_edges_per_event", 0.0))
        add_count <= 0 && (add_count = ceil(Int, 0.5n * add_per_event))
        maximum_added = Int(get(aug_cfg, "maximum_added_per_event", 3))
        added_degree = zeros(Int32, n)
        scored = Tuple{Float64,Tuple{Int32,Int32,Float64}}[]
        for edge in broad
            i, j, distance = edge
            edge_key(i, j) in selected_keys && continue
            score = comp_id[i] != comp_id[j] ? 1.0e12 / max(distance, 1.0) : (fiedler[i] - fiedler[j])^2 / (1.0 + distance)
            push!(scored, (score, edge))
        end
        sort!(scored; by=item -> item[1], rev=true)
        accepted = 0
        for (_, edge) in scored
            accepted >= add_count && break
            i, j, _ = edge
            added_degree[i] < maximum_added && added_degree[j] < maximum_added || continue
            push!(selected, edge); push!(selected_keys, edge_key(i, j)); push!(augmented, true)
            added_degree[i] += 1; added_degree[j] += 1; degree[i] += 1; degree[j] += 1
            accepted += 1
        end
    end

    edges_i = Int32[edge[1] for edge in selected]
    edges_j = Int32[edge[2] for edge in selected]
    dmetric, dxy, dxyz = Float64[], Float64[], Float64[]
    for (i, j, _) in selected
        dm, dh, d3 = edge_distances(state, points, i, j)
        push!(dmetric, dm); push!(dxy, dh); push!(dxyz, d3)
    end
    comp_id, comp_size, ncomp = graph_components(n, edges_i, edges_j)
    @printf("Built event graph: N=%d, edges=%d, k=%d, components=%d\n", n, length(edges_i), k, ncomp)
    @printf("  metric distance: p50=%.1f m, p90=%.1f m, max=%.1f m; degree min/median/max=%d/%.0f/%d\n",
        quantile_sorted(dmetric, 0.5), quantile_sorted(dmetric, 0.9), maximum(dmetric), minimum(degree), quantile_sorted(Float64.(degree), 0.5), maximum(degree))
    return EventGraph(edges_i, edges_j, dmetric, dxy, dxyz, degree, comp_id, comp_size, ncomp, augmented)
end
