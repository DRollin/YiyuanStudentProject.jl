"""
    TODO
"""
function add_bc!(ch::ConstraintHandler, grid::Grid{3}, load::LoadCase{3})
	(; ε̄, μ̄, ζ̄) = load
	∂Ω = vcat( collect_periodic_facets(grid, "left", "right"),
               collect_periodic_facets(grid, "bottom", "top"),
			   collect_periodic_facets(grid, "front", "back") )
	add!(ch, PeriodicDirichlet(:μ, ∂Ω, (x,t) -> ζ̄⋅x))
	add!(ch, PeriodicDirichlet(:u, ∂Ω, (x,t) -> ε̄⋅x)) #0.0, [1,2,3])) # TODO: why zero and not ε̄ ⋅ x ?
	#add!(ch, PeriodicDirichlet(:c, ∂Ω, (x,t) -> [1])) 
	return ch
end

_get_ref_shape(::Val{1}) = RefLine
_get_ref_shape(::Val{2}) = RefTriangle
_get_ref_shape(::Val{3}) = RefTetrahedron

"""
    TODO
"""
function prepare_setup(rve::RVE{dim}) where {dim}
    (; grid, P, M) = rve
	
	refshape   = _get_ref_shape(Val(dim))
	setP, setM = getcellset(grid, "particles"), getcellset(grid, "matrix")
		
	ip = (u = Lagrange{refshape,1}()^dim,
	      μ = Lagrange{refshape,1}(),
		  c = Lagrange{refshape,1}())

	dh = DofHandler(grid)
	add!(dh, :u, ip.u)
	add!(dh, :μ, ip.μ)
	add!(dh, :c, ip.c)
	close!(dh)

	qr  = QuadratureRule{refshape}(3)
	cv  = (u = CellValues(qr, ip.u),
	       μ = CellValues(qr, ip.μ),
		   c = CellValues(qr, ip.c))
	nbf = (u = getnbasefunctions(cv.u),
	       μ = getnbasefunctions(cv.μ),
		   c = getnbasefunctions(cv.c))

	ch = ConstraintHandler(dh)
	add_bc!(ch, grid, LoadCase(dim))
	close!(ch)

	Kₑ = zeros(sum(nbf), sum(nbf))
	Mₑ = deepcopy(Kₑ)
	submatrices = (
		Kₑuu = @view(Kₑ[dof_range(dh, :u), dof_range(dh, :u)]),
		Kₑuμ = @view(Kₑ[dof_range(dh, :u), dof_range(dh, :μ)]),
		Kₑμu = @view(Kₑ[dof_range(dh, :μ), dof_range(dh, :u)]),
		Kₑuc = @view(Kₑ[dof_range(dh, :u), dof_range(dh, :c)]),
		Kₑcu = @view(Kₑ[dof_range(dh, :c), dof_range(dh, :u)]),
		Kₑμμ = @view(Kₑ[dof_range(dh, :μ), dof_range(dh, :μ)]),
		Kₑμc = @view(Kₑ[dof_range(dh, :μ), dof_range(dh, :c)]),
		Kₑcμ = @view(Kₑ[dof_range(dh, :c), dof_range(dh, :μ)]),
		Kₑcc = @view(Kₑ[dof_range(dh, :c), dof_range(dh, :c)]),
		Mₑuu = @view(Mₑ[dof_range(dh, :u), dof_range(dh, :u)]),
		Mₑuμ = @view(Mₑ[dof_range(dh, :u), dof_range(dh, :μ)]),
		Mₑμu = @view(Mₑ[dof_range(dh, :μ), dof_range(dh, :u)]),
		Mₑuc = @view(Mₑ[dof_range(dh, :u), dof_range(dh, :c)]),
		Mₑcu = @view(Mₑ[dof_range(dh, :c), dof_range(dh, :u)]),
		Mₑμμ = @view(Mₑ[dof_range(dh, :μ), dof_range(dh, :μ)]),
		Mₑμc = @view(Mₑ[dof_range(dh, :μ), dof_range(dh, :c)]),
		Mₑcμ = @view(Mₑ[dof_range(dh, :c), dof_range(dh, :μ)]),
		Mₑcc = @view(Mₑ[dof_range(dh, :c), dof_range(dh, :c)]),
	)
 
	setups = (P = PhaseSetup{dim}(dh, setP, cv, nbf, P, Kₑ, Mₑ, submatrices), 
	          M = PhaseSetup{dim}(dh, setM, cv, nbf, M, Kₑ, Mₑ, submatrices))

	K = allocate_matrix(dh, ch)
	M = allocate_matrix(dh, ch)
	J = allocate_matrix(dh, ch)
	g    = zeros(ndofs(dh))
	aⁿ   = zeros(ndofs(dh))
    aⁿ⁺¹ = zeros(ndofs(dh))

	return RVESetup{dim}(grid, dh, setups, K, M, J, g, aⁿ, aⁿ⁺¹) 
end
