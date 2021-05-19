using NLsolve
using LinearAlgebra

# Elliptic problem for velocity
∂h_∂σ(ice, i)  = (ice.h[i + 1] - ice.h[i]) / (ice.Δx * ice.xm)
∂u_∂σ⁺(ice, i) = (ice.u[i + 1] - ice.u[i]) / (ice.Δx * ice.xm)
∂u_∂σ⁻(ice, i) = (ice.u[i] - ice.u[i - 1]) / (ice.Δx * ice.xm)
qu⁺(ice, i) = ice.h[i] * sign(∂u_∂σ⁺(ice, i)) * abs(∂u_∂σ⁺(ice, i))
qu⁻(ice, i) = ice.h[i - 1] * sign(∂u_∂σ⁻(ice, i)) * abs(∂u_∂σ⁻(ice, i))

function elliptic(ice::Ice{T}, i) where {T}
    a = 4ice.ϵ * (qu⁺(ice, i) - qu⁻(ice, i)) / (ice.Δx * ice.xm)
    b = 1 - ice.δ * ∂h_∂σ(ice, i)
    c = ice.β * sign(ice.u[i]) * abs(ice.u[i])^(1/ice.n)
    d = ice.C * sign(ice.u[i]) * abs(ice.u[i])^(1/ice.n)
    h = (ice.h[i + 1] + ice.h[i]) / 2
    return a + h * (b - c) - d
end

function bc(ice::Ice{T}) where {T}
    a = -4ice.ϵ * qu⁻(ice, ice.N+1) / (ice.Δx * ice.xm)
    b = 1 - ice.δ * ∂h_∂σ(ice, ice.N+1)
    c = ice.β * sign(ice.u[end]) * abs(ice.u[end])^(1/ice.n)
    d = ice.C * sign(ice.u[end]) * abs(ice.u[end])^(1/ice.n)
    h = ice.hm
    return a + h * (b - c) - d
end

# Hyperbolic problem for height and moving boundary
dxₘ_dt(ice) = (ice.xm - ice.xm0) / ice.Δt

function flux(ice::Ice{T}, i) where {T}
    ẋm = dxₘ_dt(ice)
    v = ice.u[i] - ice.σ[i]*ẋm
    v >= 0 ? (return v * ice.h[i]) : (return v * ice.h[i + 1])
end

function div_flux(ice::Ice{T}, i) where {T}
    return ice.Δt*(flux(ice, i) - flux(ice, i-1)) / (ice.Δx * ice.xm)
end

function source(ice::Ice{T}, i) where {T}
    return ice.Δt*(dxₘ_dt(ice) * ice.h[i]) / ice.xm
end

function solve_velocity(ice::Ice{T}) where {T}
    function f!(F, x)
        ice.u[2:end] .= x
        for i = 2:(ice.N)
            F[i-1] = elliptic(ice, i)
        end
        F[end] = bc(ice)
    end
    return nlsolve(f!, ice.u0[2:end])
end

function solve_height(ice::Ice{T}) where {T}
    function f!(F, x)
        ice.h[2:end-1] .= x[1:end-1]
        ice.xm = x[end]
        ice.h[1] = ice.h[2]
        for i in 2:ice.N+1
            F[i-1] = ice.h[i] - ice.h0[i] - div_flux(ice, i) + source(ice, i)
        end
        F[end] = ice.h[ice.N+1] - 2ice.hm
    end
    sol = nlsolve(f!, vcat(ice.h0[2:end-1], [ice.xm0]))
    return sol
end
