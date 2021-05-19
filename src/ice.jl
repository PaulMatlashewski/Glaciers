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
end

function Ice(u0, h0, xm0, params)
    u = copy(u0)
    h = copy(h0)
    xm = xm0
    σ = collect(range(0, 1, length=params.N+1))
    return Ice(σ, u, h, xm, u0, h0, xm0,
               params.hm, params.n, params.ϵ, params.δ,
               params.β, params.C, params.Δx, params.Δt, params.N)
end

function animate_ice(x, h, fps)
    n = size(x)[1]
    duration = 1/fps
    ylims = (minimum(h) - 0.01, maximum(h) + 0.01)
    xlims = (-0.01, maximum(x) + 0.01)
    gr(show=true)
    for i in 1:n
        t = time()
        p = scatter(x[i, :], h[i, :], ylims=ylims, xlims=xlims, label=nothing)
        display(p)
        Δt = time() - t
        if Δt < duration
            sleep(duration - Δt)
        end
    end
end
