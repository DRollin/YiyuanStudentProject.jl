"""
    TODO
"""
function compute_time_step!(setup::RVESetup{dim}, load::LoadCase{dim}, Δt) where{dim}
    (; grid, dh, K, M, g, J, aⁿ, aⁿ⁺¹) = setup
    ch = ConstraintHandler(dh)
	add_bc!(ch, grid, load)
	close!(ch)
        # .nzval assures that structural zeros are NOT dropped (-> needed to apply constraints)
    g .= (M .- K .* Δt ./ 2) * aⁿ # TODO: I changed the sign of K!!!
    J.nzval .= (M.nzval .+ K.nzval .* Δt ./ 2)
    apply!(J, g, ch) 
    aⁿ⁺¹ .= J \ g 
    return setup
end

# TODO: Maybe better solve_time_series ?
"""
    TODO
"""
function solve_load_case(rve::RVE{dim}, load::LoadCase{dim};  Δt=0.25, t_total=1) where {dim}
    setup = prepare_setup(rve)
    (; aⁿ, aⁿ⁺¹) = setup
    assemble_K!(setup)
    assemble_M!(setup)

    nsteps = ceil(Int, t_total/Δt)
	res = (t = Vector{Float64}(undef, nsteps+1),
           a = Vector{Vector{Float64}}(undef, nsteps+1))
    
    Vector{Tuple{Float64,Int,Vector{Float64}}}(undef, nsteps+1) # TODO: maybe rather a named Tuple of vectors? -> easier to use later
	res.t[1] = 0.0
    res.a[1] = deepcopy(aⁿ)

    for i in 1:nsteps
        compute_time_step!(setup, load, Δt)
        aⁿ .= aⁿ⁺¹
        res.t[i+1] = i*Δt
        res.a[i+1] = deepcopy(aⁿ)
    end

    return res, setup
end;
