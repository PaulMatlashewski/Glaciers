module Glaciers

include("ice.jl")
include("solver.jl")

function ice_test()
    # Test parameters
    N = 50; hm = 0.01; n = 3.0; ϵ = 1.0
    δ = 0.5; β = 0.1; C = 0.1; Δx = 1/N; Δt = Δx
    params = (hm=hm, n=n, ϵ=ϵ, δ=δ, β=β, C=C, Δx=Δx, Δt=Δt, N=N)

    # Initial data
    u0 = [-16*(x - 0.5)^2 + 4 for x in range(0, 1, length=N+1)]
    u0[1] = 0.0
    h0 = [-(5.0 - 2hm) * x^4 + 5.0 for x in range(0, 1, length=N)]
    pushfirst!(h0, h0[1])
    push!(h0, 0.0)
    xm0 = 1.0
    ice = Ice(u0, h0, xm0, params)

    # Solve for initial velocity
    solve_velocity(ice)
    ice.u0 .= ice.u

    iters = 100
    h_data = zeros(iters, N + 1)
    u_data = zeros(iters, N + 1)
    x_data = zeros(iters, N + 1)
    for i in 1:iters
        print("Step $(i) / $(iters)")
        # Solve for velocity given h
        sol_v = solve_velocity(ice)
        # Solve for h and xm given velocity
        sol_h = solve_height(ice)
        # Update previous step
        ice.u0 .= ice.u
        ice.h0 .= ice.h
        ice.xm0 = ice.xm
        # Record results
        h_data[i, :] .= ice.h[1:end-1]
        u_data[i, :] .= ice.u
        x_data[i, :] .= collect(range(0, ice.xm, length=N+1))
        println("    Velocity Converged: $(sol_v.f_converged), Height Converged: $(sol_h.f_converged)")
    end
    return h_data, u_data, x_data, ice
end

function ice_solve()
    # Test parameters
    N = 100; hm = 0.01; n = 3.0; ϵ = 1.0
    δ = 0.5; β = 0.1; C = 0.1; Δx = 1/N; Δt = Δx
    params = (hm=hm, n=n, ϵ=ϵ, δ=δ, β=β, C=C, Δx=Δx, Δt=Δt, N=N)

    # Initial data
    u0 = [-16*(x - 0.5)^2 + 4 for x in range(0, 1, length=N+1)]
    u0[1] = 0.0
    h0 = [-(5.0 - 2hm) * x^4 + 5.0 for x in range(0, 1, length=N)]
    pushfirst!(h0, h0[1])
    push!(h0, 0.0)
    xm0 = 1.0
    ice = Ice(u0, h0, xm0, params)

    # Solve for initial velocity
    solve_velocity(ice)
    ice.u0 .= ice.u

    iters = 200
    h_data = zeros(iters, N + 1)
    u_data = zeros(iters, N + 1)
    x_data = zeros(iters, N + 1)
    for i in 1:iters
        print("Step $(i) / $(iters)")
        sol = solve(ice)
        # Update previous step
        ice.u0 .= ice.u
        ice.h0 .= ice.h
        ice.xm0 = ice.xm
        # Record results
        h_data[i, :] .= ice.h[1:end-1]
        u_data[i, :] .= ice.u
        x_data[i, :] .= collect(range(0, ice.xm, length=N+1))
        println("    Converged: $(sol.f_converged), Iterations: $(sol.iterations), Inf-norm: $(sol.residual_norm)")
    end
    return h_data, u_data, x_data
end

export ice_test, ice_solve, animate_ice

end # module
