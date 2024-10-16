function solve_macro_problem(grid::Grid{dim}, problem::RVEProblem{dim},  Δt::Real; steps=0.0:10.0:1000.0) where {dim}
    dh, ch, cv, uₙ = prepare_macro_setup(grid)
    K = create_sparsity_pattern(dh, ch)
    M = create_sparsity_pattern(dh, ch)
    assemble_K!(K, dh, cv, problem)
    assemble_M!(M, dh, cv, problem)

    A = (Δt .* K) + M
    rhsdata = get_rhs_data(ch, A)
    apply!(A, ch)

    b  = zeros(ndofs(dh))
    uₙ = zeros(ndofs(dh))
    
	t = 0.0
    res = Dict{Float64,Vector}()
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
    return res, (grid=grid, dh=dh, ch=ch, cv=cv)
end
