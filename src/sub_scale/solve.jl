"""
    compute_time_step!(setup::RVESetup{dim}, load::LoadCase{dim}, Δt) 

Return the object `RVESetup` with updated next time step solution vector `aⁿ⁺¹`.

Using Crank-Nicolson Method to compute the solution `aⁿ⁺¹`  in an implicit time stepping scheme for solving 
the linear system in the next time step.

New boundary condition values `LoadCase` are applied before for the next time step solutions.
"""
function compute_time_step!(setup::RVESetup{dim}, load::LoadCase{dim}, Δt) where{dim}
    @info "Compute RVE time step"
    (; grid, dh, K, M, f, g, J, aⁿ, aⁿ⁺¹) = setup
    ch = ConstraintHandler(dh)
	add_bc!(ch, grid, load)
	close!(ch)
        # .nzval assures that structural zeros are NOT dropped (-> needed to apply constraints)
    g .= Δt .* f .+ (M .- K .* Δt ./ 2) * aⁿ
    J.nzval .= M.nzval .+ K.nzval .* Δt ./ 2
    apply!(J, g, ch) 
    aⁿ⁺¹ .= J \ g
    apply!(aⁿ⁺¹, ch) 
    
    @info "RVE time step computed"
    return setup, aⁿ⁺¹
end


"""
    solve_time_series(rve::RVE{dim}, load::LoadCase{dim};  Δt=0.1, t_total=1) where {dim}

Return results throughout the time stepping process for visualizing the fine scale problem results.  
    
The RVE problem is prepared, assembled and solved over a time series using the given `RVE`, `LoadCase`, `t_total`, `Δt` . By default `t_total` is set to `1` and `Δt` to `0.1`.

Results are stored in a `NamedTuple` with fields total time `t`  and solution vector `a`.

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
