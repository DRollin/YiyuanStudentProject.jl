"""
    compute_time_step!(setup::RVESetup{dim}, load::LoadCase{dim}, Δt) 

Compute the solution for the next time step in a Representative Volume Element (RVE) simulation using an implicit time integration scheme.

# Arguments
- `setup::RVESetup{dim}`:   The setup object containing parameters for the RVE simulation
- `load::LoadCase{dim}`:    The load case specifying boundary conditions and external loads for the current time step.
- `Δt::Real`:               The time step size.


# Implementation Details
This function uses Crank-Nicolson Method to compute the solution in an implicit time stepping scheme for solving 
the linear system in the next time step:

The function applies boundary conditions to the system matrix and force vector using a constraint handler before solving 
for the next time step solutions.

"""
function compute_time_step!(setup::RVESetup{dim}, load::LoadCase{dim}, Δt) where{dim}
    @info "Compute RVE time step"
    (; grid, dh, K, M, f, g, J, aⁿ, aⁿ⁺¹) = setup
    ch = ConstraintHandler(dh)
	add_bc!(ch, grid, load)
	close!(ch)
        # .nzval assures that structural zeros are NOT dropped (-> needed to apply constraints)
    g .= Δt .* f .+ (M .- K .* Δt ./ 2) * aⁿ
    J.nzval .= (M.nzval .+ K.nzval .* Δt ./ 2 .+ 1e-10)
    apply!(J, g, ch) 
    aⁿ⁺¹ .= J \ g
    apply!(aⁿ⁺¹, ch) 
    
    @info "RVE time step computed"
    return setup, aⁿ⁺¹
end


"""
    solve_time_series(rve::RVE{dim}, load::LoadCase{dim};  Δt=0.25, t_total=1) where {dim}

For visualizing the fine scale problem results.
Compute the results in a `NamedTuple` with fields `t` total time and `a` solution vector series contain `u` displacement, `μ` chemical potantial, 
and `c` concentration for the whole time series with a certain time step width. 

# Arguments
- `rve`:    Object for solving function `prepare_setup`,
- `load`:   Object `LoadCase` with a 

"""
function solve_time_series(rve::RVE{dim}, load::LoadCase{dim};  Δt=0.1, t_total=1) where {dim}
    setup = prepare_setup(rve)
    (; aⁿ, aⁿ⁺¹) = setup
    assemble_K_M_f!(setup)
    

    nsteps = ceil(Int, t_total/Δt)
	res = (t = Vector{Float64}(undef, nsteps+1),
           a = Vector{Vector{Float64}}(undef, nsteps+1))
    
    
	res.t[1] = 0.0
    res.a[1] = deepcopy(aⁿ)
    @withprogress name="Time stepping" begin
        for i in 1:nsteps
            compute_time_step!(setup, load, Δt)
            aⁿ .= aⁿ⁺¹
            res.t[i+1] = i*Δt
            res.a[i+1] = deepcopy(aⁿ)
            @logprogress i/nsteps
        end
    end

    return res, setup
end;
