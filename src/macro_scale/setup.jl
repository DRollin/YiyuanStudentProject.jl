"""
	add_macro_bc!(ch::ConstraintHandler, grid::Grid{3}, Δt, bc::MacroBCParams)

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
    prepare_macro_setup(grid::Grid{dim,C}, rvesetup::RVESetup{dim}, rve::RVE{dim}, Δt, bc::MacroBCParams) where {dim,C}

Return a `SolveSetup` object for the element-wise assembly and solving the time dependent macro scale problem. 

A initialized buffer `data` is created with the object `GaussPointData` at each integration point. 

The object `PhaseSetup` is created based on the given macro scale grid `grid`, the boundary condition `MacroBCParams` with the buffer `data`.

With necessary parameters, initialized global matrix and vectors, and `PhaseSetup` the `SolveSetup` is computed.

Two unknown fields deformation `u`, chemical potential `μ` are prescribed for all cells.

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
