using NLsolve
using Plots

struct FVGrid{T}
    u::Vector{T}
    u0::Vector{T}
    v::Vector{T}
    Δx::T
    Δt::T
    N::Int
end

struct FDGrid{T}
    u::Vector{T}
    u0::Vector{T}
    v::Vector{T}
    Δx::T
    Δt::T
    N::Int
end

# Implicit time steps
function implicit_solve(F, iters)
    function f!(G, x)
        F.u[2:end-1] .= x
        for i in 2:F.N+1
            G[i-1] = F.u[i] - F.u0[i] + div(F, i)
        end
    end
    gr(show=true)
    for i in 2:iters
        println("Step $(i) / $(iters) Mass: $(sum(F.u))")
        nlsolve(f!, F.u0[2:end-1])
        F.u0 .= F.u
        p = scatter(F.u[2:end-1], ylims=(-1.2, 2.2), label=nothing)
        display(p)
    end
end

# Explicit time steps
function explicit_solve(F, iters)
    gr(show=true)
    for i in 2:iters
        println("Step $(i) / $(iters) Mass: $(sum(F.u))")
        for j in 2:F.N+1
            F.u[j] = F.u0[j] - div(F, j)
        end
        F.u0 .= F.u
        p = scatter(F.u[2:end-1], ylims=(-1.2, 2.2), label=nothing)
        display(p)
    end
end

###### Finite Volume ######
# TODO: Need to return flux instead of velocity times field
function flux(F, i)
    uL = F.u[i]
    uR = F.u[i + 1]
    vL = F.v[i]
    vR = F.v[i + 1]
    # Upwind
    vL > 0 && vR > 0 && return vL * uL
    vL < 0 && vR < 0 && return vR * uR
    # Rarefaction
    vL <= 0 && vR >= 0 && return 0.0
    # Shock
    ṡ = 0.5 * (vR*uR - vL*uL) / (uR - uL) # Rankine-Hugonoit
    ṡ > 0 && return vL * uL
    ṡ < 0 && return vR * uR
    return 0.0
end

function div(F::FVGrid{T}, i) where {T}
    return (flux(F, i) - flux(F, i-1)) * F.Δt / F.Δx
end

###### Finite Difference ######

function div(F::FDGrid{T}, i) where {T}
    uL = F.u[i - 1]
    uC = F.u[i]
    uR = F.u[i + 1]
    vL = F.v[i - 1]
    vC = F.v[i]
    vR = F.v[i + 1]
    v̂L = 0.5 * (vC + vL)
    v̂R = 0.5 * (vC + vR)
    # Upwind
    v̂L > 0 && v̂R > 0 && return (uC * vC - uL * vL) * F.Δt / F.Δx
    v̂L < 0 && v̂R < 0 && return (uR * vR - uC * vC) * F.Δt / F.Δx
    # Rarefaction
    v̂L <= 0 && v̂R >= 0 && return 0.0
    # Shock
    ṡ = 0.5 * (vR*uR - vL*uL) / (uR - uL) # Rankine-Hugonoit
    ṡ > 0 && return (uC * vC - uL * vL) * F.Δt / F.Δx
    ṡ < 0 && return (uR * vR - uC * vC) * F.Δt / F.Δx
    return 0.0
end

function hyperbolic_test()
    # Test parameters
    N = 100
    iters = 300
    α = 0.05

    # Initial data (CLAWPACK example)
    # function n_wave(x)
    #     if x > -pi && x < pi
    #         return (cos(x) + 1) * (2sin(3x) + cos(2x) + 0.2)
    #     else
    #         return 0.0
    #     end
    # end
    # u = n_wave.(range(-8, 8, length=N+2))
    # Initial data
    u = zeros(N + 2)
    c = Int(round(N/2))
    w = Int(round(N/8))
    # u[(c-w):(c+w)] .= -1.0
    u[1:c] .= -1
    u[c+1:end] .= 2
    v = u
    Δx = 1 / N
    Δt = α * Δx / maximum(abs.(v))
    
    # Solve
    F = FDGrid(u, copy(u), v, Δx, Δt, N)
    implicit_solve(F, iters)
end
