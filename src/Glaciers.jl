module Glaciers

include("ice.jl")
include("solver.jl")

function ice_solve()
    # Test parameters
    N = 100; hm = 0.01; n = 3.0; ϵ = 0.05; tsteps = 700
    δ = 0.05; β = 0.1; C = 0.1; Δx = 1/N; Δt = Δx
    params = (hm=hm, n=n, ϵ=ϵ, δ=δ, β=β, C=C, Δx=Δx, Δt=Δt, N=N, tsteps = 700)

    # Initial data
    u0 = [-16*(x - 0.5)^2 + 4 for x in range(0, 1, length=N+1)]
    u0[1] = 0.0
    h0 = [-(0.2 - 2hm) * x^4 + 0.2 for x in range(0, 1, length=N)]
    pushfirst!(h0, h0[1])
    push!(h0, 0.0)
    xm0 = 1.0
    ice = Ice(u0, h0, xm0, params; store_trace=true)

    # Solve for initial velocity
    solve_velocity(ice)
    ice.u0 .= ice.u
    ice.Δt = ice.Δx / maximum(ice.u)

    # Solve problem
    solve!(ice)
    return ice
end

export ice_solve, animate_height, animate_velocity

end # module
