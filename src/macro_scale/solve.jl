
"""
    solve_macro_problem(grid::Grid{dim}, rvesetup::RVESetup{dim}; Δt=0.1, t_total=1) where {dim}

return solution storage and the macro scale setup for result visualizing. 
Compute the solution storage in a `NamedTuple` with fields `t` total time and `a` solution vector series contain `u` displacement and `μ` chemical potantial
for the whole time series with a certain time step width. 

# Arguments
- `grid`:       Grid{dim}: The macro scale grid,
- `rvesetup`:   RVESetup{dim}: The object containing setups for solving rve problem,
- `Δt`:         The time step size,
- `t-total`:    the total computational time.


# Implementation Details
Initialize the solution storage. Do time stepping with updating the boundary condition.

The updated boundary condition is then applied to the global stiffness matrix and force vector as well as the solution vector.

Using solver from Pardiso.jl for eventually ill conditioned system matrices. Update the boundary condition on the result vector.

Copy the new solution vector in the storage with the corresponding time step.

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
