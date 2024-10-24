δ(i,j) = i == j ? 1 : 0

function add_bc!(ch::ConstraintHandler, grid::Grid, load::LoadCase{3})
	(; μ̄, ζ̄) = load
	∂Ω = vcat( collect_periodic_facets(grid, "left", "right"),
               collect_periodic_facets(grid, "bottom", "top"),
			   collect_periodic_facets(grid, "front", "back") )

	
	add!(ch, PeriodicDirichlet(:μ, ∂Ω, (x,t) -> μ̄ + ζ̄⋅x))
	add!(ch, PeriodicDirichlet(:u, ∂Ω, [1,2,3]))

	return ch
end

function prepare_base_setup(grid::Grid{dim}) where {dim}
	
    bshape = RefTetrahedron    #RefTetrahedron

	setP, setM = get_phase_cell_sets(grid)
		
	ip  = (u = Lagrange{bshape,1}()^(dim), μ = Lagrange{bshape,1}(), c = μ = Lagrange{bshape,1}())

	dh = DofHandler(grid)
	add!(dh, :u, ip.u)
	add!(dh, :μ, ip.μ)
	add!(dh, :c, ip.c)

	close!(dh)

	qr  = QuadratureRule{bshape}(4)
	cv  = (u = CellValues(qr, ip.u), μ = CellValues(qr, ip.μ), c = CellValues(qr, ip.c))
	nbf = (u = getnbasefunctions(cv.u), μ = getnbasefunctions(cv.μ), c = getnbasefunctions(cv.c))
	

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
