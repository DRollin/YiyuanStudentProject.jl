"""
TODO

"""
function add_macro_bc!(ch::ConstraintHandler, grid::Grid{3}, Δt)
	∂Ωₗ = getfacetset(grid, "left") 
    ∂Ωᵣ = getfacetset(grid, "right")
	ramp(t) = t < 10*Δt ? t/(10*Δt) : 1.0
    add!(ch, Dirichlet(:u, ∂Ωₗ, (x,t) -> (ramp(t)*1.0, 0.0, 0.0) ))
    add!(ch, Dirichlet(:u, ∂Ωᵣ, (x,t) -> (0.0, 0.0, 0.0)))
    add!(ch, Dirichlet(:μ, ∂Ωₗ, (x,t) -> ramp(t)*1.0))
    add!(ch, Dirichlet(:μ, ∂Ωᵣ, (x,t) -> 0.0)) 
	return ch
end

"""
TODO

"""
function prepare_macro_setup(grid::Grid{dim}, rvesetup::RVESetup{dim}, Δt) where {dim}
	@info "Preparing Macro setup"
    refshape   = _get_ref_shape(Val(dim))
		
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
    add_macro_bc!(ch, grid, Δt)
	close!(ch)

    Kₑ = zeros(sum(nbf), sum(nbf))      
	fₑ = zeros(sum(nbf))
	aₑ = deepcopy(fₑ)

    subarrays = (
		Kₑuu = @view(Kₑ[dof_range(dh, :u), dof_range(dh, :u)]),
		Kₑuμ = @view(Kₑ[dof_range(dh, :u), dof_range(dh, :μ)]),
		Kₑμu = @view(Kₑ[dof_range(dh, :μ), dof_range(dh, :u)]),
		Kₑμμ = @view(Kₑ[dof_range(dh, :μ), dof_range(dh, :μ)]),
	)


	K = allocate_matrix(dh, ch)
	f = zeros(ndofs(dh))
	aⁿ = deepcopy(f)

	
    data = [[ GaussPointData(rvesetup.aⁿ,
							 (P = GaussPointPhaseData(0.0, zero(Tensor{1,dim})),
							  M = GaussPointPhaseData(0.0, zero(Tensor{1,dim})))
							 )
                for j in 1:getnquadpoints(cv.u)]
                for i in 1:length(grid.cells)]

    setups = AssemblySetup{dim}(dh, cv, nbf, Kₑ, subarrays, data, aₑ)

    setup = SolveSetup{dim}(grid, dh, ch, setups, K, f, aⁿ) 

	@info "Macro setup prepared"
	return setup
end