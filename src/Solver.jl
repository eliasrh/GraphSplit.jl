struct LinearizedSystem
    i::Vector{Int32}
    j::Vector{Int32}
    gix::Vector{Float64}
    giy::Vector{Float64}
    giz::Vector{Float64}
    gjx::Vector{Float64}
    gjy::Vector{Float64}
    gjz::Vector{Float64}
    sqrt_weight::Vector{Float64}
    n_events::Int
    gauge_mode::Symbol
    constraint_weight::Float64
    reference_velocity_ms::Float64
end

function copy_state(state::State)
    return State(copy(state.event_id), copy(state.x), copy(state.y), copy(state.z), copy(state.t0),
        state.ref_lat, state.ref_lon, state.ref_radius_m)
end

function predict_and_gradients(state::State, stations::Stations, obs::Observations,
        travel_time::AbstractTravelTimeModel)
    nobs = length(obs)
    prediction = Vector{Float64}(undef, nobs)
    gix, giy, giz = zeros(nobs), zeros(nobs), zeros(nobs)
    gjx, gjy, gjz = zeros(nobs), zeros(nobs), zeros(nobs)
    for row in 1:nobs
        i, j, station = obs.i[row], obs.j[row], obs.station[row]
        ti, ix, iy, iz = travel_time_gradient(travel_time, obs.phase[row], state.x[i], state.y[i], state.z[i],
            stations.x_m[station], stations.y_m[station], stations.z_m[station])
        tj, jx, jy, jz = travel_time_gradient(travel_time, obs.phase[row], state.x[j], state.y[j], state.z[j],
            stations.x_m[station], stations.y_m[station], stations.z_m[station])
        isfinite(ti) && isfinite(tj) || error("Travel-time lookup returned a non-finite value for observation $row")
        prediction[row] = state.t0[i] - state.t0[j] + ti - tj
        gix[row], giy[row], giz[row] = ix, iy, iz
        gjx[row], gjy[row], gjz[row] = jx, jy, jz
    end
    return prediction, gix, giy, giz, gjx, gjy, gjz
end

huber_weight(value::Float64, threshold::Float64) = abs(value) <= threshold ? 1.0 : threshold / abs(value)

function gauge_definition(cfg::AbstractDict)
    mode = lowercase(String(cfgget(cfg, "gauge", "mode"; default="zero_mean")))
    if mode == "zero_mean"
        value = lowercase(String(cfgget(cfg, "gauge", "zero_mean"; default="origin_time")))
        value in ("origin_time", "t0_only") && return :origin_time
        value in ("xyzt0", "xyzt0_unscaled", "t0_xyz") && return :xyzt0
        value in ("xyzt0_scaled", "t0_xyz_scaled") && return :xyzt0_scaled
        error("Unknown gauge.zero_mean mode: $value")
    elseif mode == "pin"
        fields = lowercase(String(cfgget(cfg, "gauge", "pin_fields"; default="xyz")))
        return occursin("t0", fields) || occursin('t', fields) ? :none : :origin_time
    end
    error("gauge.mode must be zero_mean or pin")
end

gauge_rows(system::LinearizedSystem) = system.gauge_mode == :none ? 0 : system.gauge_mode == :origin_time ? 1 : 4

function apply_A(system::LinearizedSystem, model_update::Vector{Float64})
    n, nobs = system.n_events, length(system.i)
    output = zeros(nobs + gauge_rows(system))
    for row in 1:nobs
        i, j = system.i[row], system.j[row]
        left = system.gix[row] * model_update[i] + system.giy[row] * model_update[n + i] +
            system.giz[row] * model_update[2n + i] + model_update[3n + i]
        right = system.gjx[row] * model_update[j] + system.gjy[row] * model_update[n + j] +
            system.gjz[row] * model_update[2n + j] + model_update[3n + j]
        output[row] = system.sqrt_weight[row] * (left - right)
    end
    rows = gauge_rows(system)
    rows == 0 && return output
    scale = sqrt(system.constraint_weight)
    if system.gauge_mode == :origin_time
        output[nobs + 1] = scale * sum(@view model_update[3n + 1:4n]) / n
    else
        spatial_scale = system.gauge_mode == :xyzt0_scaled ? system.reference_velocity_ms : 1.0
        output[nobs + 1] = scale * sum(@view model_update[1:n]) / (n * spatial_scale)
        output[nobs + 2] = scale * sum(@view model_update[n + 1:2n]) / (n * spatial_scale)
        output[nobs + 3] = scale * sum(@view model_update[2n + 1:3n]) / (n * spatial_scale)
        output[nobs + 4] = scale * sum(@view model_update[3n + 1:4n]) / n
    end
    return output
end

function apply_At(system::LinearizedSystem, data::Vector{Float64})
    n, nobs = system.n_events, length(system.i)
    length(data) == nobs + gauge_rows(system) || error("Linear operator data length mismatch")
    output = zeros(4n)
    for row in 1:nobs
        i, j = system.i[row], system.j[row]
        q = system.sqrt_weight[row] * data[row]
        output[i] += q * system.gix[row]; output[n + i] += q * system.giy[row]
        output[2n + i] += q * system.giz[row]; output[3n + i] += q
        output[j] -= q * system.gjx[row]; output[n + j] -= q * system.gjy[row]
        output[2n + j] -= q * system.gjz[row]; output[3n + j] -= q
    end
    rows = gauge_rows(system)
    rows == 0 && return output
    scale = sqrt(system.constraint_weight) / n
    if system.gauge_mode == :origin_time
        value = scale * data[nobs + 1]
        @views output[3n + 1:4n] .+= value
    else
        spatial_scale = system.gauge_mode == :xyzt0_scaled ? system.reference_velocity_ms : 1.0
        @views output[1:n] .+= scale * data[nobs + 1] / spatial_scale
        @views output[n + 1:2n] .+= scale * data[nobs + 2] / spatial_scale
        @views output[2n + 1:3n] .+= scale * data[nobs + 3] / spatial_scale
        @views output[3n + 1:4n] .+= scale * data[nobs + 4]
    end
    return output
end

function pinned_active_mask(state::State, cfg::AbstractDict)
    active = trues(4length(state))
    lowercase(String(cfgget(cfg, "gauge", "mode"; default="zero_mean"))) == "pin" || return active
    pin_ids = Int64.(cfgget(cfg, "gauge", "pin_event_ids"; default=Int[]))
    isempty(pin_ids) && error("gauge.mode=pin requires at least one gauge.pin_event_ids entry")
    event_rows = Dict(id => row for (row, id) in enumerate(state.event_id))
    missing = setdiff(pin_ids, state.event_id)
    isempty(missing) || error("Pinned event IDs were not found in the catalog: $(join(missing, ','))")
    fields = lowercase(String(cfgget(cfg, "gauge", "pin_fields"; default="xyz")))
    n = length(state)
    for id in pin_ids
        row = event_rows[id]
        occursin('x', fields) && (active[row] = false)
        occursin('y', fields) && (active[n + row] = false)
        occursin('z', fields) && (active[2n + row] = false)
        (occursin("t0", fields) || fields == "t") && (active[3n + row] = false)
    end
    return active
end

function block_jacobi(system::LinearizedSystem, damping::Float64, active::BitVector)
    n = system.n_events
    blocks = zeros(4, 4, n)
    for row in eachindex(system.i)
        weight = system.sqrt_weight[row]^2
        vi = (system.gix[row], system.giy[row], system.giz[row], 1.0)
        vj = (system.gjx[row], system.gjy[row], system.gjz[row], 1.0)
        for a in 1:4, b in 1:4
            blocks[a, b, system.i[row]] += weight * vi[a] * vi[b]
            blocks[a, b, system.j[row]] += weight * vj[a] * vj[b]
        end
    end
    inverses = similar(blocks)
    for event in 1:n
        block = Matrix(@view blocks[:, :, event])
        for d in 1:4
            block[d, d] += damping
            global_index = (d - 1) * n + event
            if !active[global_index]
                block[d, :] .= 0.0; block[:, d] .= 0.0; block[d, d] = 1.0
            end
        end
        inverses[:, :, event] .= inv(Symmetric(block))
    end
    return inverses
end

function diagonal_jacobi(system::LinearizedSystem, damping::Float64, active::BitVector)
    n = system.n_events
    diagonal = fill(damping, 4n)
    for row in eachindex(system.i)
        weight = system.sqrt_weight[row]^2
        i, j = system.i[row], system.j[row]
        diagonal[i] += weight * system.gix[row]^2; diagonal[n + i] += weight * system.giy[row]^2
        diagonal[2n + i] += weight * system.giz[row]^2; diagonal[3n + i] += weight
        diagonal[j] += weight * system.gjx[row]^2; diagonal[n + j] += weight * system.gjy[row]^2
        diagonal[2n + j] += weight * system.gjz[row]^2; diagonal[3n + j] += weight
    end
    diagonal[.!active] .= 1.0
    return 1.0 ./ max.(diagonal, eps(Float64))
end

function preconditioner(system::LinearizedSystem, damping::Float64, active::BitVector, name::String)
    n = system.n_events
    if name == "block_jacobi"
        inverse_blocks = block_jacobi(system, damping, active)
        return function (residual::Vector{Float64})
            result = similar(residual)
            for event in 1:n
                for a in 1:4
                    value = 0.0
                    for b in 1:4
                        value += inverse_blocks[a, b, event] * residual[(b - 1) * n + event]
                    end
                    result[(a - 1) * n + event] = value
                end
            end
            return result
        end
    elseif name == "diagonal"
        inverse_diagonal = diagonal_jacobi(system, damping, active)
        return residual -> inverse_diagonal .* residual
    end
    return residual -> copy(residual)
end

function pcg(operator, right_hand_side::Vector{Float64}, apply_preconditioner;
        tolerance::Float64=1.0e-5, maximum_iterations::Int=250)
    x = zeros(length(right_hand_side))
    residual = copy(right_hand_side)
    z = apply_preconditioner(residual)
    direction = copy(z)
    rz = dot(residual, z)
    norm_rhs = max(norm(right_hand_side), eps(Float64))
    norm(residual) <= tolerance * norm_rhs && return x, 0, true
    for iteration in 1:maximum_iterations
        product = operator(direction)
        denominator = dot(direction, product)
        denominator > 0.0 && isfinite(denominator) || return x, iteration, false
        alpha = rz / denominator
        x .+= alpha .* direction
        residual .-= alpha .* product
        norm(residual) <= tolerance * norm_rhs && return x, iteration, true
        z = apply_preconditioner(residual)
        rz_new = dot(residual, z)
        beta = rz_new / rz
        direction .= z .+ beta .* direction
        rz = rz_new
    end
    return x, maximum_iterations, false
end

function direct_solve(operator, right_hand_side::Vector{Float64}, active::BitVector)
    free = findall(active)
    length(free) <= 4000 || error("The dependency-free direct solver is limited to 4000 free parameters; use PCG")
    matrix = Matrix{Float64}(undef, length(free), length(free))
    basis = zeros(length(right_hand_side))
    for (column, index) in enumerate(free)
        fill!(basis, 0.0); basis[index] = 1.0
        matrix[:, column] .= operator(basis)[free]
    end
    result = zeros(length(right_hand_side))
    result[free] = Symmetric(0.5 .* (matrix .+ matrix')) \ right_hand_side[free]
    return result, 0, true
end

function apply_step_limits!(update::Vector{Float64}, n::Int, solver_cfg::AbstractDict)
    maximum_spatial = Float64(get(solver_cfg, "max_event_step_m", 0.0))
    maximum_time = Float64(get(solver_cfg, "max_origin_step_s", 0.0))
    if maximum_spatial > 0.0
        for event in 1:n
            length3 = sqrt(update[event]^2 + update[n + event]^2 + update[2n + event]^2)
            if length3 > maximum_spatial
                scale = maximum_spatial / length3
                update[event] *= scale; update[n + event] *= scale; update[2n + event] *= scale
            end
        end
    end
    if maximum_time > 0.0
        time_view = view(update, 3n + 1:4n)
        time_view .= clamp.(time_view, -maximum_time, maximum_time)
    end
    return update
end

"Robust matrix-free Gauss-Newton/IRLS relocation for either stage."
function solve_relocation(initial::State, stations::Stations, obs::Observations,
        travel_time::AbstractTravelTimeModel, solver_cfg::AbstractDict, cfg::AbstractDict)
    isempty(obs) && return copy_state(initial), SolveStats(0, Float64[], Float64[], Float64[], Int[], true)
    state = copy_state(initial)
    n = length(state)
    active = pinned_active_mask(state, cfg)
    maximum_outer = Int(get(solver_cfg, "max_outer_iterations", 20))
    minimum_outer = Int(get(solver_cfg, "min_outer_iterations", 3))
    minimum_sigma = Float64(get(solver_cfg, "min_sigma_s", 0.002))
    huber_k = Float64(get(solver_cfg, "huber_k", 1.345))
    damping = max(Float64(get(solver_cfg, "damping_lambda", 0.0)), 1.0e-8)
    step_damping = Float64(get(solver_cfg, "step_damping", 1.0))
    inner_tolerance = Float64(get(solver_cfg, "inner_tolerance", 1.0e-5))
    inner_maximum = Int(get(solver_cfg, "inner_max_iterations", 250))
    preconditioner_name = lowercase(String(get(solver_cfg, "preconditioner", "block_jacobi")))
    linear_solver = lowercase(String(get(solver_cfg, "linear_solver", "pcg")))
    rms_history, spatial_history, time_history = Float64[], Float64[], Float64[]
    inner_history = Int[]
    stall = 0
    converged = false
    for iteration in 1:maximum_outer
        prediction, gix, giy, giz, gjx, gjy, gjz = predict_and_gradients(state, stations, obs, travel_time)
        residual = obs.dt .- prediction
        sigma = max.(obs.sigma, minimum_sigma)
        robust = [huber_weight(residual[row] / sigma[row], huber_k) for row in eachindex(residual)]
        weight = robust ./ (sigma .^ 2)
        sqrt_weight = sqrt.(weight)
        system = LinearizedSystem(obs.i, obs.j, gix, giy, giz, gjx, gjy, gjz, sqrt_weight, n,
            gauge_definition(cfg), Float64(cfgget(cfg, "gauge", "constraint_weight"; default=10.0)),
            Float64(cfgget(cfg, "gauge", "reference_velocity_ms"; default=5000.0)))
        data = vcat(sqrt_weight .* residual, zeros(gauge_rows(system)))
        rhs = apply_At(system, data)
        rhs[.!active] .= 0.0
        operator = function (vector::Vector{Float64})
            work = copy(vector); work[.!active] .= 0.0
            result = apply_At(system, apply_A(system, work)) .+ damping .* work
            result[.!active] .= vector[.!active]
            return result
        end
        apply_prec = preconditioner(system, damping, active, preconditioner_name)
        if linear_solver == "direct"
            update, inner_iterations, inner_converged = direct_solve(operator, rhs, active)
        else
            update, inner_iterations, inner_converged = pcg(operator, rhs, apply_prec;
                tolerance=inner_tolerance, maximum_iterations=inner_maximum)
        end
        if !inner_converged
            @warn @sprintf("Inner linear solve did not converge at outer iteration %d after %d inner iterations",
                iteration, inner_iterations)
        end
        update .*= step_damping
        apply_step_limits!(update, n, solver_cfg)
        update[.!active] .= 0.0
        dx = view(update, 1:n)
        dy = view(update, n + 1:2n)
        dz = view(update, 2n + 1:3n)
        dt = view(update, 3n + 1:4n)
        state.x .+= dx; state.y .+= dy; state.z .+= dz; state.t0 .+= dt
        if gauge_definition(cfg) == :origin_time
            state.t0 .-= sum(state.t0) / n
        end
        robust_rms = sqrt(sum(weight .* residual .^ 2) / max(sum(weight), eps(Float64)))
        spatial_rms = sqrt(sum(dx .^ 2 .+ dy .^ 2 .+ dz .^ 2) / n)
        time_rms = sqrt(sum(dt .^ 2) / n)
        push!(rms_history, robust_rms); push!(spatial_history, spatial_rms); push!(time_history, time_rms); push!(inner_history, inner_iterations)
        Bool(get(solver_cfg, "verbose", true)) && @printf("  iter %3d: robust RMS %.6f s; step RMS %.3f m, %.3f ms; PCG %d\n",
            iteration, robust_rms, spatial_rms, 1000 * time_rms, inner_iterations)
        if iteration >= minimum_outer
            small_step = spatial_rms <= Float64(get(solver_cfg, "stop_step_rms_m", 0.1)) &&
                time_rms <= Float64(get(solver_cfg, "stop_step_rms_s", 1.0e-4))
            if small_step
                converged = true
                break
            end
            if length(rms_history) > 1
                improvement = rms_history[end - 1] - rms_history[end]
                abs(improvement) <= Float64(get(solver_cfg, "stop_rms_improvement_s", 1.0e-5)) ? (stall += 1) : (stall = 0)
                maximum_stall = Int(get(solver_cfg, "stop_stall_iterations", 5))
                if maximum_stall > 0 && stall >= maximum_stall
                    converged = true
                    break
                end
            end
        end
    end
    return state, SolveStats(length(rms_history), rms_history, spatial_history, time_history, inner_history, converged)
end
