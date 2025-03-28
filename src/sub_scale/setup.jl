"""
    add_bc!(ch::ConstraintHandler, grid::Grid{3}, load::LoadCase{3})

Add a Dirichlet boundary condition for the unknown fields respectively on the `∂Ω` part of the boundary. `∂Ω` is defined by the collection of facets from the `grid`. 
Boundary condition values are given by object`LoadCase`. 
"""
function add_bc!(ch::ConstraintHandler, grid::Grid{3}, load::LoadCase{3})
	(; ε̄, μ̄, ζ̄) = load
	
	∂Ω = union(getfacetset.([grid], ["left", "right", "bottom", "top", "back", "front"])...)
	add!(ch, Dirichlet(:μ, ∂Ω, (x,t) -> μ̄ + ζ̄⋅x))
	add!(ch, Dirichlet(:u, ∂Ω, (x,t) -> ε̄⋅x))

	return ch
end

_get_ref_shape(::Val{1}) = RefLine
_get_ref_shape(::Val{2}) = RefTriangle
_get_ref_shape(::Val{3}) = RefTetrahedron

function get_volume(grid, cellvalues)
    V  = 0.0 
    for c in CellIterator(grid)
        reinit!(cellvalues, c)
        for qp in 1:getnquadpoints(cellvalues)
            dΩ = getdetJdV(cellvalues, qp)
            V += dΩ
        end
    end
    return V
end

"""
    prepare_setup(rve::RVE{dim}) where {dim}

Return a `RVESetup` object for the element-wise assembly and time stepping. 
	
With the object `RVE` as argument, the `RVESetup` is created by constructing `PhaseSetup` objects for material phase particles `P` and matrix `M`, initializing 
necessary global matrices and vectors, and the total volume of RVE.

Three unknown fields deformation `u`, chemical potential `μ` and ion concentration `c` are prescribed for all cells over both material phases.

"""
function prepare_setup(rve::RVE{dim}) where {dim}
    @info "Preparing RVE setup"
	(; grid, P, M) = rve
	
	P_material = P
	M_material = M

	
	refshape   = _get_ref_shape(Val(dim))
	Ωᴾ, Ωᴹ = getcellset(grid, "particles"), getcellset(grid, "matrix")
		
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
	fₑ = zeros(sum(nbf))
	subarrays = (
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
		fₑu  = @view(fₑ[dof_range(dh, :u)]),
		fₑμ  = @view(fₑ[dof_range(dh, :μ)]),
		fₑc  = @view(fₑ[dof_range(dh, :c)]),
	)
 
	setups = (P = PhaseSetup{dim}(dh, Ωᴾ, cv, nbf, P, Kₑ, Mₑ, fₑ, subarrays), 
	          M = PhaseSetup{dim}(dh, Ωᴹ, cv, nbf, M, Kₑ, Mₑ, fₑ, subarrays))

	K = allocate_matrix(dh, ch)
	M = deepcopy(K)
	J = deepcopy(K)
	f = zeros(ndofs(dh))
	g    = deepcopy(f)
	aⁿ   = deepcopy(f)
    aⁿ⁺¹ = deepcopy(f)

	apply_analytical!(aⁿ, dh, :c, (x -> P_material.cʳᵉᶠ), Ωᴾ)
	apply_analytical!(aⁿ, dh, :c, (x -> M_material.cʳᵉᶠ), Ωᴹ)
	apply_analytical!(aⁿ, dh, :μ, (x -> P_material.μʳᵉᶠ), Ωᴾ)
	apply_analytical!(aⁿ, dh, :μ, (x -> M_material.μʳᵉᶠ), Ωᴹ)

	Vʳᵛᵉ = get_volume(grid, cv.u)
	setup = RVESetup{dim}(grid, dh, setups, K, M, f, J, g, aⁿ, aⁿ⁺¹, Vʳᵛᵉ) 
	@info "RVE Setup prepared"
	return setup
end


