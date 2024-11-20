function compute_time_step!(p::cm_Problem, Δt) 
    (; setup, K, M, aⁿ, aⁿ⁺¹) = p 
    (; ch, g, J) = setup
    # .nzval assures that structural zeros are NOT dropped (-> needed to apply constraints)

    dim = 3
    g .= (M ./ Δt .+ K ./ 2) * aⁿ
    J.nzval .= (M.nzval ./ Δt + K.nzval ./ 2)
    apply!(J, g, ch) 
    aⁿ⁺¹ .= J \ g 
    aⁿ .= aⁿ⁺¹

    return cm_Problem{dim}(setup, K, M, aⁿ, aⁿ⁺¹), -1
end


function solve_load_case(problem::RVEProblem{dim}, load::LoadCase{dim};  Δt=0.25, t_total=1) where {dim}
    setup = prepare_setup(problem, load)
	K, M= assemble_rve_system(setup)
    aⁿ = zeros(ndofs(setup.dh))
	aⁿ⁺¹ = zeros(ndofs(setup.dh))
    @show ndofs(setup.dh)

	cm_problem = cm_Problem{dim}(setup, K, M, aⁿ, aⁿ⁺¹)
    nsteps = ceil(Int, t_total/Δt)
	res = Vector{Tuple{Float64,Int,Vector{Float64}}}(undef, nsteps+1)
	res[1] = (0.0, 0, deepcopy(aⁿ))


    for i in 1:nsteps
        @show norm(aⁿ)
        cm_problem, niter = compute_time_step!(cm_problem, Δt)
        res[i+1] = (i*Δt, niter, deepcopy(aⁿ))
    end

    return res, setup.dh
end;

function solve_RVE()
    dim_rve = 3
    #Material properties
    G = 4.0e5
    K = 6.67e5
    α = 1.0
    η = 1.0
    μ_ref = 1.0
    c_m = 1.0
    θ_ref = 1.0
    c_ref = 1.0
    R = 8.31446261815324
    δ(i,j) = i == j ? 1.0 : 0.0

    P = Material_pre(G, K, α, R, θ_ref, c_m, c_ref, η, μ_ref, dim_rve)
    
    M = Material_pre(G, K, α, R, θ_ref, c_m, c_ref, η, μ_ref, dim_rve)
	

    #RVE problem setups
    Load = LoadCase_pre(dim_rve, 1.0, 1.0, 1.0, 1.0)

    #Grid
    d = 0.1 
    ϕ = 0.1 #packing fraction
    meshsize = 0.1
    
    rve_grid= generate_rve_grid(; ϕ=ϕ, d=d, meshsize=meshsize, dx=(1.0,1.0, 1.0))

    rve_problem = RVEProblem{dim_rve}(rve_grid, P, M)

    res, dh = solve_load_case(rve_problem, Load)

    plot_grid(rve_grid)

    plot_result(res, dh, "μ")




end
