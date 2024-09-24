δ(i,j) = i == j ? 1 : 0

function add_bc!(ch::ConstraintHandler, grid::Grid, load::LoadCase{2})
	(; μ̄, ζ̄) = load
	∂Ω = vcat( collect_periodic_faces(grid, "left", "right"),
               collect_periodic_faces(grid, "bottom", "top") )
	add!(ch, PeriodicDirichlet(:u, ∂Ω, (x,t) -> μ̄ + ζ̄⋅x))
	return ch
end

function prepare_base_setup(grid::Grid{dim}) where {dim}
	if dim == 2
        bshape = RefTetrahedron
    end
	setP, setM = get_phase_cell_sets(grid)
		
	ip  = Lagrange{2,bshape,1}()

	dh = DofHandler(grid)
	push!(dh, :u, 1, ip)
	close!(dh)

	qr  = QuadratureRule{dim, bshape}(4)
	cv  = CellScalarValues(qr, ip)
	nbf = Ferrite.getnbasefunctions(cv)

	ch = ConstraintHandler(dh)
	return dh, ch, cv, nbf, setP, setM
end
	
function prepare_setup(problem::RVEProblem{dim}, Load::LoadCase{dim}) where {dim}
    (; grid, P, M) = problem

	dh, ch, cv, nbf, setP, setM = prepare_base_setup(grid)
	add_bc!(ch, grid, Load)
	close!(ch)
	update!(ch)

	setupP = BulkElementSetup{dim}(dh, cv, nbf, P)
	setupM = BulkElementSetup{dim}(dh, cv, nbf, M)

	sets   = (P=setP,   M=setM)
	setups = (P=setupP, M=setupM)
	return FESetup{dim}(grid, dh, ch, Load, sets, setups) 
end

function prepare_dns_setup(grid::Grid{dim}, materials::NamedTuple) where {dim}
	dh, ch, cv, nbf, setP, setM = prepare_base_setup(grid)
	
	∂Ωₗ = getfaceset(grid, "left") 
    ∂Ωᵣ = getfaceset(grid, "right")
	add!(ch, Dirichlet(:u, ∂Ωₗ, (x,t) -> ramp(t, 10.0, 1.0)))
    add!(ch, Dirichlet(:u, ∂Ωᵣ, (x,t) -> 0.0))
	close!(ch)
	update!(ch)

	setupP = BulkElementSetup{dim}(dh, cv, nbf, materials.P, 0.0, zero(Vec{dim}), zero(Tensor{2,dim}))
	setupM = BulkElementSetup{dim}(dh, cv, nbf, materials.M, 0.0, zero(Vec{dim}), zero(Tensor{2,dim}))

	sets   = (P=setP,   M=setM)
	setups = (P=setupP, M=setupM)
	return FESetup{dim}(grid, dh, ch, nothing, sets, setups) 
end