function solve_macro_problem(grid::Grid{dim}, rvesetup::RVESetup{dim}; Δt=0.1, t_total=1) where {dim}
    setup = prepare_macro_setup(grid, rvesetup, Δt)

    (; grid, ch, K, f, aⁿ) = setup
    
    nsteps = ceil(Int, t_total/Δt)
	res = (t = Vector{Float64}(undef, nsteps+1),
           a = Vector{Vector{Float64}}(undef, nsteps+1))

    res.t[1] = 0.0
    res.a[1] = deepcopy(aⁿ)
    
    @withprogress name="Time stepping" begin
        for i in 1:nsteps
            update!(ch, i*Δt)

            K, setup.Assemblysetup.data = assemble_macro_K!(setup, rvesetup, Δt)
            
            apply!(K, f, ch)
            
            aⁿ .= K \ f
            apply!(aⁿ, ch)
            
            res.t[i+1] = i*Δt
            res.a[i+1] = deepcopy(aⁿ)
            @logprogress i/nsteps
        end  
    end 
    
    return res, setup
end
