δ(i,j) = i == j ? 1 : 0

function add_bc!(ch::ConstraintHandler, grid::Grid, load::LoadCase{3})
	(; μ̄, ζ̄) = load
	∂Ω = vcat( collect_periodic_facets(grid, "left", "right"),
               collect_periodic_facets(grid, "bottom", "top"),
			   collect_periodic_facets(grid, "front", "back") )

	displacement_bc = (x, t) -> μ̄ .+ ζ̄ .* x
	add!(ch, PeriodicDirichlet(:u, ∂Ω, displacement_bc))
	add!(ch, PeriodicDirichlet(:p, ∂Ω, (x,t) -> 0.0))
	return ch
end

function prepare_base_setup(grid::Grid{dim}) where {dim}
	if dim == 3
        bshape = RefTetrahedron
    end
	setP, setM = get_phase_cell_sets(grid)
		
	ip  = (u = Lagrange{bshape,1}()^(dim), p = Lagrange{bshape,1}())

	dh = DofHandler(grid)
	add!(dh, :u, ip.u)
	add!(dh, :p, ip.p)
	#push!(dh, :u, 2, ip)
	#push!(dh, :p, 1, ip)
	close!(dh)

	qr  = QuadratureRule{bshape}(4)
	cv  = (u = CellValues(qr, ip.u), p = CellValues(qr, ip.p))
	nbf = (u = getnbasefunctions(cv.u), p = getnbasefunctions(cv.p))
	@show nbf
	@show ndofs(dh)

	ch = ConstraintHandler(dh)
	return dh, ch, cv, nbf, setP, setM
end
	
function prepare_setup(problem::RVEProblem{dim}, Load::LoadCase{dim}) where {dim}
    (; grid, P, M) = problem

	dh, ch, cv, nbf, setP, setM = prepare_base_setup(grid)
	add_bc!(ch, grid, Load)
	close!(ch)
	update!(ch)

	setupP = iso_pv_ElementSetup{dim}(dh, cv, nbf, P)
	setupM = iso_pv_ElementSetup{dim}(dh, cv, nbf, M)

	sets   = (P=setP,   M=setM)
	setups = (P=setupP, M=setupM)
	return FESetup_base{dim}(grid, dh, ch, Load, sets, setups) 
end
