"""
    add_bc!(ch::ConstraintHandler, grid::Grid{3}, load::LoadCase{3})

Create a periodic boundary condition on `u` and `μ` unknown fields respectively on the `∂Ω` part of the boundary. `∂Ω` is defined by the collection of periodic facets from the `grid`. 
"""
function add_bc!(ch::ConstraintHandler, grid::Grid{3}, load::LoadCase{3})
	(; ε̄, μ̄, ζ̄) = load
	∂Ω = vcat( collect_periodic_facets(grid, "left", "right"),
               collect_periodic_facets(grid, "bottom", "top"),
			   collect_periodic_facets(grid, "front", "back") )
	add!(ch, PeriodicDirichlet(:μ, ∂Ω, (x,t) -> ζ̄⋅x))
	add!(ch, PeriodicDirichlet(:u, ∂Ω, (x,t) -> ε̄⋅x))

	centernode_idx = argmin((idx_node) -> norm(idx_node[2].x), enumerate(grid.nodes))[1]

	centernode = OrderedSet{Int}([centernode_idx])

	#centernode =  OrderedSet{Int}([ argmin(n -> norm(n.x), grid.nodes) ])
	add!(ch, Dirichlet(:u, centernode, (x,t) -> zero(Vec{3})))
	#add!(ch, Dirichlet(:μ, centernode, (x,t) -> μ̄))
	return ch
end

_get_ref_shape(::Val{1}) = RefLine
_get_ref_shape(::Val{2}) = RefTriangle
_get_ref_shape(::Val{3}) = RefTetrahedron

"""
    prepare_setup(rve::RVE{dim}) where {dim}

Return a `RVESetup` struct for the elementweise assembly and time stepping. 
	
# Arguments:
-`rve::RVE{dim}`: An `RVE` object containing the following fields:

- `grid`: 		The grid (mesh) for RVE.
- `P`: 			Phase data for "particles" in the material.
- `M`: 			Phase data for the "matrix" in the material.

# Implementation Details:
The interpolation `ip` is defined by passing the corresponding `refshape` using the function `_get_ref_shape(Val(dim))`. 

After the 'DofHandler' is created based on the `grid` from the argument `rve`, unknown fields can be added to it using `add!`
By calling `close!` to finalize the construction. 

Cellvallues are created for each discrete fields respectively passing quadrature rule with 3 integration points.

Number of base function for each unknown fields is defined using `getnbasefunctions`.

Calling function `add_bc!` to initialize the boundary conditions and associate this to the dofhandler.

Local stiffness and mass matrices are initialized. subarrays for locating the corresponding field interaction in elementweise by calling the macro `@view`.

A named tuple with the Struct `PhaseSetup` is then prepared for each `P` stands for particals and `M` for matrix.

Gloable stiffness, mass, and jacobian matrices is initialized using `allocate_matrix` for a sparse matrix pattern.

initialize the gloable residual vector and the solution vectors for both current and next time step.
"""
function prepare_setup(rve::RVE{dim}) where {dim}
    @info "Preparing setup"
	(; grid, P, M) = rve
	
	P_material = P
	M_material = M

	
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
	Cₑ = zeros(nbf.μ)
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
 
	setups = (P = PhaseSetup{dim}(dh, setP, cv, nbf, P, Kₑ, Mₑ, Cₑ, fₑ, subarrays), 
	          M = PhaseSetup{dim}(dh, setM, cv, nbf, M, Kₑ, Mₑ, Cₑ, fₑ, subarrays))

		# System matrix with Lagrange multiplier for <μ> = μ̄
	dofsμ = Set{Int}([ dof for cell in 1:getncells(grid) for dof in celldofs(dh, cell)[dof_range(dh, :μ)] ])
	dofsμ = sort!(collect(dofsμ))
	C = sparse(ones(Int, length(dofsμ)), dofsμ, zeros(length(dofsμ)), 1, ndofs(dh))
	K = allocate_matrix(dh, ch)
	K = hcat( vcat(K, C), copy(hcat(C, spzeros(1,1))') )

	M = deepcopy(K)
	J = deepcopy(K)
	f = zeros(ndofs(dh) + 1)
	g    = deepcopy(f)
	aⁿ   = deepcopy(f)
    aⁿ⁺¹ = deepcopy(f)

	apply_analytical!(aⁿ, dh, :c, (x -> P_material.cʳᵉᶠ), setP)
	apply_analytical!(aⁿ, dh, :c, (x -> M_material.cʳᵉᶠ), setM)
	apply_analytical!(aⁿ, dh, :μ, (x -> P_material.μʳᵉᶠ), setP)
	apply_analytical!(aⁿ, dh, :μ, (x -> M_material.μʳᵉᶠ), setM)

	Vʳᵛᵉ = get_volume(grid, cv.u)
	
	setup = RVESetup{dim}(grid, dh, setups, K, M, f, J, g, aⁿ, aⁿ⁺¹, Vʳᵛᵉ) 
	@info "Setup prepared"
	return setup
end


function get_volume(grid::Grid, cv::CellValues)
    V = 0.0
    for cc in CellIterator(grid)
		reinit!(cv, cc)
		for qp in 1:getnquadpoints(cv)
        	V +=  getdetJdV(cv, qp)
		end
    end
    return V
end