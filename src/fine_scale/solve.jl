function solve_load_case(problem::RVEProblem{dim}, load::LoadCase{dim}) where {dim}
    setup = prepare_setup(problem, load)
    K, f  = assemble_rve_system(setup)
	u = K\f
	u = u[1:ndofs(setup.dh)]
	apply!(u, setup.ch)
	return u, setup
end

function solve_dns_problem(grid::Grid{dim}, materials::NamedTuple, Δt::Real; steps=0.0:10.0:1000.0) where {dim}
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
end