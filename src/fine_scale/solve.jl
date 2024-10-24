function solve_load_case(problem::RVEProblem{dim}, load::LoadCase{dim};  Δt=0.025, t_total=1.0) where {dim}
    setup = prepare_setup(problem, load)
	K = allocate_matrix(setup.dh, setup.ch)
    f = zeros(ndofs(setup.dh))
    a = zeros(ndofs(setup.dh))
	a_old = copy(a)

	pv_problem = pv_Problem{dim}(setup, K, f, a, a_old)

	pvd = paraview_collection("porous_media")
    step = 0

    for t in 0:Δt:t_total
        if t>0
            update!(setup.ch, t)
            apply!(a, setup.ch)
            doassemble_K_f!(pv_problem, Δt)
            apply_zero!(K, f, setup.ch)
			
            Δa = -K\f
            apply_zero!(Δa, setup.ch)
            a .+= Δa
            copyto!(a_old, a)
        end
        step += 1
        VTKGridFile("porous_media_$step", setup.dh) do vtk
            write_solution(vtk, setup.dh, a)
            pvd[t] = vtk
        end
    end

    vtk_save(pvd);
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
	#@show P

    #RVE problem setups
    Load = LoadCase_pre(dim_rve, 1.0, 1.0, 1.0, 1.0)

    #Grid
    d = 0.1 
    ϕ = 0.1 #packing fraction
    meshsize = 0.05
    
    rve_grid= generate_rve_grid(; ϕ=ϕ, d=d, meshsize=meshsize, dx=(1.0,1.0,1.0))

    rve_problem = RVEProblem{dim_rve}(rve_grid, P, M)

    solve_load_case(rve_problem, Load)

    
end
