function solve_load_case(problem::RVEProblem{dim}, load::LoadCase{dim}) where {dim}
    setup = prepare_setup(problem, load)
	K = create_sparsity_pattern(setup.dh, setup.ch)
    f = zeros(ndofs(setup.dh))
    a = zeros(ndofs(setup.dh))
	a_old = copy(a)
	

	pv_problem = pv_Problem{dim}(setup, K, f, a, a_old)

	pvd = paraview_collection("porous_media")
    step = 0
    for t in 0:Δt:t_total
        if t>0
            update!(ch, t)
            apply!(a, ch)
            doassemble_K_f!(pv_problem, Δt)
            apply_zero!(K, f, ch)
            Δa = -K\f
            apply_zero!(Δa, ch)
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

#=function solve_dns_problem(grid::Grid{dim}, materials::NamedTuple, Δt::Real; steps=0.0:10.0:1000.0) where {dim}
	setup = prepare_dns_setup(grid, materials)
	(; dh, ch) = setup
	M, K, b = assemble_dns_system(setup)

	A = (Δt .* K) + M
    rhsdata = get_rhs_data(ch, A)
    apply!(A, ch)

    uₙ = zeros(ndofs(dh))
    
	t = 0.0
	res = Dict{Float64, Vector}()
	for step in steps
		while t < step
			t += Δt
			update!(ch, t)
			b .= M*uₙ
			apply_rhs!(rhsdata, b, ch)
			uₙ .= A\b
		end
		merge!(res, Dict([t => deepcopy(uₙ)]))
	end
    return res, setup
end=#

function solve_RVE()
    dim_rve = 3 
    #Material properties
    G = 80e3
    K = 160e3
    κ = 0.05
    α = 0.9
    β = 1/15e3

    P = Material_pre(G, K, κ, α, β, dim_rve)
    M = Material_pre(G, K, κ, α, β, dim_rve)

    #RVE problem setups
    Load = LoadCase(1.0, 0.0, dim_rve)

    #Grid
    d = 0.1 
    ϕ = 0.1 #packing fraction
    meshsize = 0.05
    
    rve_grid= generate_rve_grid(; ϕ=ϕ, d=d, meshsize=meshsize, dx=(1.0,1.0))

    rve_problem = RVEProblem{dim_rve}(rve_grid, P, M)

    solve_load_case(rve_problem, Load)

    
end


#u, setup, grid, problem = solve_RVE_prepare()

#Δt = 1

#@time res, (grid1, dh1, ch1, cv1) = solve_macro_problem(grid, problem, Δt, steps=0.0:10.0:1000.0)

#t, u = select_state(res, 5.0)

#vtk_grid("Macro_potential_$t", dh1) do vtk
#    vtk_point_data(vtk, dh1, u, "u")
#end