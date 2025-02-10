"""
   add_macro_bc!(ch::ConstraintHandler, grid::Grid{3}, Δt)

Add a time dependent dirichlet boundary condition for the macro scale unknown fields respectively on the possible facet `∂Ω1` and `∂Ω2`.

Apply time dependent boundary condition on both unknown fields `u` and `μ`.

"""
function add_macro_bc!(ch::ConstraintHandler, grid::Grid{3}, Δt, bc::MacroBCParams)
	∂Ω1 = getfacetset(grid, bc.facetset_1) 
    ∂Ω2 = getfacetset(grid, bc.facetset_2)
	ramp(t) = t < 10*Δt ? t/(10*Δt) : 1.0
    add!(ch, Dirichlet(:u, ∂Ω1, (x,t) -> ramp(t)*bc.u_1 ))
    add!(ch, Dirichlet(:u, ∂Ω2, (x,t) -> ramp(t)*bc.u_2))
    add!(ch, Dirichlet(:μ, ∂Ω1, (x,t) -> ramp(t)*bc.μ_1))
    add!(ch, Dirichlet(:μ, ∂Ω2, (x,t) -> ramp(t)*bc.μ_2)) 
	return ch
end


"""
    prepare_macro_setup(grid::Grid{dim}, rvesetup::RVESetup{dim}, Δt) where {dim}

Return a `SolveSetup` object for the element-wise assembly and solving the time dependent macro scale problem. 
	
# Arguments:
- `grid`:		prescribed macro scale grid with no distinguished phases,
- `rvesetup`:	[RVESetup](@ref "RVE{dim}"),
- `Δt`:			The time step size.

# Implementation Details:
The interpolation `ip` is defined by passing the corresponding `refshape` using the function `_get_ref_shape(Val(dim))`. 

After the 'DofHandler' is created based on the `grid`, unknown fields can be added to it using `add!`
By calling `close!` to finalize the construction. 

Cellvallues are created for each discrete fields respectively passing quadrature rule with 3 integration points.

Number of base function for each unknown fields is defined using `getnbasefunctions`.

Calling function `add_bc!` to initialize the boundary conditions and associate this to the ``Dofhandler``.

Local stiffness and solution vector are initialized. subarrays for locating the corresponding field interaction in element-wise by calling the macro `@view`.


Global stiffness matrix, global right-hand side and global solution vector are initialized. For the stiffness matrix using `allocate_matrix` for a sparse matrix pattern.

Create the necessary data storage ``GaussPointData`` for each gauss point in macro scale with the initialized solution vector from `RVESetup`, `0.0` for `c` concentration and `c₂` the gradient of concentration.

A ``AssemblySetup`` for elementweise assembly and ``SolveSetup`` for the final solving of the time dependent problem are prepared.
"""
function prepare_macro_setup(grid::Grid{dim,C}, rvesetup::RVESetup{dim}, rve::RVE{dim}, Δt, bc::MacroBCParams) where {dim,C}
	@info "Preparing Macro setup"
    refshape = getrefshape(C)
		
	ip = (u = Lagrange{refshape,1}()^dim,
	      μ = Lagrange{refshape,1}())

    dh = DofHandler(grid)
    add!(dh, :u, ip.u)
	add!(dh, :μ, ip.μ)  
	
	close!(dh)

    qr  = QuadratureRule{refshape}(3)

	cv  = (u = CellValues(qr, ip.u),
	       μ = CellValues(qr, ip.μ))
	nbf = (u = getnbasefunctions(cv.u),
	       μ = getnbasefunctions(cv.μ))


    ch = ConstraintHandler(dh)
    add_macro_bc!(ch, grid, Δt, bc)
	close!(ch)

    Kₑ = zeros(sum(nbf), sum(nbf))      
	#fₑ = zeros(sum(nbf))
	aₑ = zeros(sum(nbf))

    subarrays = (
		Kₑuu = @view(Kₑ[dof_range(dh, :u), dof_range(dh, :u)]),
		Kₑuμ = @view(Kₑ[dof_range(dh, :u), dof_range(dh, :μ)]),
		Kₑμu = @view(Kₑ[dof_range(dh, :μ), dof_range(dh, :u)]),
		Kₑμμ = @view(Kₑ[dof_range(dh, :μ), dof_range(dh, :μ)]),
	)


	K = allocate_matrix(dh, ch)
	f = zeros(ndofs(dh))
	aⁿ = deepcopy(f)

	apply_analytical!(aⁿ, dh, :μ, (x -> rve.M.μʳᵉᶠ), 1:getncells(grid))

    data = [[ GaussPointData(rvesetup.aⁿ, 0.0, zero(Tensor{1,dim}))
                for j in 1:getnquadpoints(cv.u)]
                for i in 1:length(grid.cells)]

    assemblysetup = AssemblySetup{dim}(dh, cv, nbf, Kₑ, aₑ, subarrays, data, rvesetup)
    setup = SolveSetup{dim}(grid, dh, ch, assemblysetup, K, f, aⁿ) 

	@info "Macro setup prepared"
	return setup
end
