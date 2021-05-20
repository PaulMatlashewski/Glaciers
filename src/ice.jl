using Plots

mutable struct Ice{T}
    σ::Vector{T}
    u::Vector{T}
    h::Vector{T}
    xm::T
    u0::Vector{T}
    h0::Vector{T}
    xm0::T
    hm::T
    n::T
    ϵ::T
    δ::T
    β::T
    C::T
    Δx::T
    Δt::T
    N::Int
    tsteps::Int
    colors::Vector{Int}
    store_trace::Bool
    u_trace::Matrix{T}
    h_trace::Matrix{T}
    xm_trace::Vector{T}
    c_trace::Matrix{Int}
end

function Ice(u0, h0, xm0, params; store_trace=false)
    u = copy(u0)
    h = copy(h0)
    xm = xm0
    σ = collect(range(0, 1, length=params.N+1))
    colors = zeros(Int, length(h))
    if store_trace
        u_trace = zeros(length(u), params.tsteps)
        h_trace = zeros(length(h), params.tsteps)
        xm_trace = zeros(params.tsteps)
        c_trace = zeros(Int, length(h), params.tsteps)
    else
        u_trace = zeros(0, 0)
        h_trace = zeros(0, 0)
        xm_trace = zeros(0)
        c_trace = zeros(Int, 0)
    end
    return Ice(σ, u, h, xm, u0, h0, xm0,
               params.hm, params.n, params.ϵ, params.δ, params.β,
               params.C, params.Δx, params.Δt, params.N, params.tsteps,
               colors, store_trace, u_trace, h_trace, xm_trace, c_trace)
end

function animate_height(ice::Ice{T}; stretch_coord=false, fps=60) where {T}
    N = ice.N
    duration = 1/fps
    ylims = (minimum(ice.h_trace) - 0.01, maximum(ice.h_trace) + 0.01)
    x = zeros(N, ice.tsteps)
    # Stretch coordinates or fixed coordinates
    if stretch_coord
        xlims = (-0.01, 1.01)
        xs = 0.5 * (ice.σ[1:end-1] .+ ice.σ[2:end])
        for i in 1:ice.tsteps
            x[:, i] .= xs
        end
    else
        xlims = (-0.01, maximum(ice.xm_trace) + 0.01)
        for i in 1:ice.tsteps
            xs = collect(range(0, ice.xm_trace[i], length=N+1))
            x[:, i] .= 0.5 * (xs[1:end-1] .+ xs[2:end])
        end
    end
    # Color markers by flux category
    cmap = Dict(
        0 => RGBA(0.0000, 0.0000, 0.0000, 1.0), # NA
        1 => RGBA(0.0000, 0.6056, 0.9787, 1.0), # ++
        2 => RGBA(0.8889, 0.4356, 0.2781, 1.0), # --
        3 => RGBA(0.2422, 0.6433, 0.3044, 1.0), # +-
        4 => RGBA(0.7644, 0.4441, 0.8243, 1.0)  # -+
    )
    c_rgba = [[cmap[color] for color in data] for data in eachcol(ice.c_trace[2:end-1, :])]
    gr(show=true)
    p = scatter(-1*ones(1, 4), -1*ones(1, 4), ylims=ylims, xlims=xlims, label=["++" "--" "+-" "-+"])
    scatter!(x[:, 1], ice.h_trace[2:end-1, 1], label=nothing, color=c_rgba[1])
    for i in 1:ice.tsteps
        t = time()
        p.series_list[5].plotattributes[:x] .= x[:, i]
        p.series_list[5].plotattributes[:y] .= ice.h_trace[2:end-1, i]
        p.series_list[5].plotattributes[:markercolor] .= c_rgba[i]
        display(p)
        Δt = time() - t
        if Δt < duration
            sleep(duration - Δt)
        end
    end
end

function animate_velocity(ice::Ice{T}; stretch_coord=false, fps=60) where {T}
    N = ice.N
    duration = 1/fps
    ylims = (minimum(ice.u_trace) - 0.01, maximum(ice.u_trace) + 0.01)
    x = zeros(N + 1, ice.tsteps)
    # Stretch coordinates or fixed coordinates
    if stretch_coord
        xlims = (-0.01, 1.01)
        for i in 1:ice.tsteps
            x[:, i] .= ice.σ
        end
    else
        xlims = (-0.01, maximum(ice.xm_trace) + 0.01)
        for i in 1:ice.tsteps
            xs = collect(range(0, ice.xm_trace[i], length=N+1))
            x[:, i] .= xs
        end
    end
    gr(show=true)
    p = scatter(ylims=ylims, xlims=xlims)
    scatter!(x[:, 1], ice.u_trace[:, 1], label=nothing)
    for i in 1:ice.tsteps
        t = time()
        p.series_list[1].plotattributes[:x] .= x[:, i]
        p.series_list[1].plotattributes[:y] .= ice.u_trace[:, i]
        display(p)
        Δt = time() - t
        if Δt < duration
            sleep(duration - Δt)
        end
    end
end
