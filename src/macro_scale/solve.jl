
"""
    solve_macro_problem(grid::Grid{dim}, rvesetup::RVESetup{dim}, rve::RVE{dim}, bc::MacroBCParams; Δt=0.1, t_total=1) where {dim}

Return results throughout the time stepping process for visualizing the fine scale problem results.  
    
The macro scale problem is prepared, assembled and solved over a time series using the given `grid`, `RVESetup`, `RVE`, `MacroBCParams`, `t_total`, `Δt` . By default `t_total` is set to `1` and `Δt` to `0.1`.

Results on both macro scale and the reference RVE are stored in a `NamedTuple` with fields total time `t`  and solution vector `a` respectively.

The time dependent boundary condition is updated to the global stiffness matrix and force vector as well as the solution vector.

Using solver from ``Pardiso.jl`` for eventually ill conditioned system matrices.

"""
function solve_macro_problem(grid::Grid{dim}, rvesetup::RVESetup{dim}, rve::RVE{dim}, bc::MacroBCParams; Δt=0.1, t_total=1) where {dim}
    setup = prepare_macro_setup(grid, rvesetup, rve, Δt, bc)
    (; grid, ch, K, f, aⁿ, assemblysetup) = setup
    (; gpdata) = assemblysetup
    
    nsteps = ceil(Int, t_total/Δt)
	res = (t = Vector{Float64}(undef, nsteps+1),
           a = Vector{Vector{Float64}}(undef, nsteps+1))

    res_rve = (t = Vector{Float64}(undef, nsteps+1),
               a = Vector{Vector{Float64}}(undef, nsteps+1))

    res.t[1] = 0.0
    res.a[1] = deepcopy(aⁿ)

    res_rve.t[1] = 0.0
    res_rve.a[1] = deepcopy(gpdata[1][1].aᵣᵥₑⁿ)
    
    @withprogress name="Time stepping" begin
        for i in 1:nsteps
            tⁿ⁺¹ = i*Δt
            update!(ch, tⁿ⁺¹)

            assemble_macro_K!(setup, Δt)

            apply!(K, f, ch)

            # Solver for ill-conditioned matrices
            prob = LinearSolve.LinearProblem(K, f)
            sol  = LinearSolve.solve(prob, LinearSolve.MKLPardisoFactorize())
            aⁿ .= sol.u
            cacheval = sol.cache.cacheval
            Pardiso.set_phase!(cacheval, Pardiso.RELEASE_ALL)
            Pardiso.pardiso(cacheval)

            #aⁿ .= K \ f

            apply!(aⁿ, ch)

            res.t[i+1] = tⁿ⁺¹
            res.a[i+1] = deepcopy(aⁿ)
            res_rve.t[i+1] = tⁿ⁺¹
            res_rve.a[i+1] = deepcopy(gpdata[1][1].aᵣᵥₑⁿ)
            @logprogress i/nsteps
        end  
    end 
    
    return res, res_rve, setup
end
