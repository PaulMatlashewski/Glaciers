using NLsolve
using Plots
using LinearAlgebra

struct Grid{T}
    u::Vector{T}
    u0::Vector{T}
    h::Vector{T}
    xm::T
    hm::T
    n::T
    ϵ::T
    δ::T
    β::T
    C::T
    Δx::T
    N::Int
end

function Grid(u, h, xm, params)
    return Grid(copy(u), copy(u), h, xm, params.hm, params.n, params.ϵ,
               params.δ, params.β, params.C, params.Δx, params.N)
end

Δ⁺(u, i) = u[i + 1] - u[i]
Δ⁻(u, i) = u[i] - u[i - 1]

function σ1(ice, i)
    u = ice.u
    h = ice.h
    xm = ice.xm
    n = ice.n
    ϵ = ice.ϵ
    Δx = ice.Δx

    du_dx1 = Δ⁻(u, i) / (Δx * xm)
    du_dx2 = Δ⁺(u, i) / (Δx * xm)
    ν1 = sign(du_dx1) * abs(du_dx1)^(1/n)
    ν2 = sign(du_dx2) * abs(du_dx2)^(1/n)
    h̄1 = (h[i] + h[i - 1]) / 2
    h̄2 = (h[i] + h[i + 1]) / 2
    return 4ϵ/(Δx * xm) * (h̄2*ν2 - h̄1*ν1)
end

function σ2(ice, i)
    u = ice.u
    h = ice.h
    xm = ice.xm
    n = ice.n
    δ = ice.δ
    β = ice.β
    C = ice.C
    Δx = ice.Δx

    dh_dx = (h[i] - h[i - 1]) / (Δx * xm)
    a = 1 - δ * dh_dx
    b = β * sign(u[i]) * abs(u[i])^(1/n)
    c = C * sign(u[i]) * abs(u[i])^(1/n)
    return h[i] * (a - b) - c
end

function elliptic(ice, i)
    return σ1(ice, i) + σ2(ice, i)
end

function bc(ice)
    u = ice.u
    h = ice.h
    xm = ice.xm
    n = ice.n
    ϵ = ice.ϵ
    δ = ice.δ
    β = ice.β
    C = ice.C
    Δx = ice.Δx

    du_dx = (u[end] - u[end - 1]) / (Δx * xm)
    dh_dx = (h[end] - h[end - 1]) / (Δx * xm)
    h̄ = (h[end] + h[end - 1]) / 2
    ν = sign(du_dx) * abs(du_dx)^(1/n)
    a = (h[end]^2 / 2 - 4ϵ*h̄*ν) / (Δx * xm)
    b = 1 - δ * dh_dx
    c = β * sign(u[end]) * abs(u[end])^(1/n)
    d = C * sign(u[end]) * abs(u[end])^(1/n)
    return a + h[end] * (b - c) - d
end

function solve(ice)
    function f!(F, x)
        ice.u[2:end] .= x
        for i = 2:(ice.N)
            F[i-1] = elliptic(ice, i)
        end
        F[end] = bc(ice)
    end
    return nlsolve(f!, ice.u0[2:end], store_trace=true)
end

function elliptic_test()
    # Test parameters
    N = 30
    params = (hm=0.01, n=3.0, ϵ=0.01, δ=0.05, β=0.1, C=0.2, Δx=1/N, N=N)

    # Initial data
    u = [-16*(x - 0.5)^2 + 4 for x in range(0, 1, length=N+1)]
    u[1] = 0.0
    h = [-(1.0 - params.hm) * x^4 + 1.0 for x in range(0, 1, length=N+1)]
    xm = 1.0

    # Solve
    ice = Grid(u, h, xm, params)
    sol = solve(ice)
    F = zeros(N)
    for i = 2:(ice.N)
        F[i-1] = elliptic(ice, i)
    end
    F[end] = bc(ice)
    return sol, ice, F
end
