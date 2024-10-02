δ(i,j) = i == j ? 1 : 0

function add_bc!(ch::ConstraintHandler, grid::Grid, load::LoadCase{2})
	(; μ̄, ζ̄) = load
	∂Ω = vcat( collect_periodic_faces(grid, "left", "right"),
               collect_periodic_faces(grid, "bottom", "top") )
	add!(ch, PeriodicDirichlet(:u, ∂Ω, (x,t) -> μ̄ + ζ̄⋅x))
	add!(ch, PeriodicDirichlet(:p, ∂Ω, (x,t) -> 0.0))
	return ch
end

function prepare_base_setup(grid::Grid{dim}) where {dim}
	if dim == 2
        bshape = RefTetrahedron
    end
	setP, setM = get_phase_cell_sets(grid)
		
	ip  = (u = Lagrange{2,bshape,1}()^(dim), p = Lagrange{2,bshape,1}())

	dh = DofHandler(grid)
	add!(dh, :u, ip.u)
	add!(dh, :p, ip.p)
	#push!(dh, :u, 2, ip)
	#push!(dh, :p, 1, ip)
	close!(dh)

	qr  = QuadratureRule{dim, bshape}(4)
	cv  = (u = CellVectorValues(qr, ip.u), p = CellScalarValues(qr, ip.p))
	nbf = (u = getnbasefunctions(cv.u), p = getnbasefunctions(cv.p))

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
	return FESetup{dim}(grid, dh, ch, Load, sets, setups) 
end
