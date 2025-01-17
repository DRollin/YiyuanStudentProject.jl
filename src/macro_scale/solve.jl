function solve_macro_problem(grid::Grid{dim}, rvesetup::RVESetup{dim}; Δt=0.1, t_total=1) where {dim}
    setup = prepare_macro_setup(grid, rvesetup, Δt)
    (; grid, ch, K, f, aⁿ, assemblysetup) = setup
    #(; gpdata) = assemblysetup
    
    nsteps = ceil(Int, t_total/Δt)
	res = (t = Vector{Float64}(undef, nsteps+1),
           a = Vector{Vector{Float64}}(undef, nsteps+1))
    (; t, a) = res
    t[1] = 0.0
    a[1] = deepcopy(aⁿ)
    
    @withprogress name="Time stepping" begin
        for i in 1:nsteps
            tⁿ⁺¹ = i*Δt
            update!(ch, tⁿ⁺¹)

            #K, setup.Assemblysetup.data = 
            assemble_macro_K!(setup, Δt)

            #fdata = get_rhs_data(ch, K)
            apply!(K, f, ch)
            #apply_rhs!(fdata, f, ch) 
            aⁿ .= K \ f
            apply!(aⁿ, ch)
            #@show aⁿ
            t[i+1] = tⁿ⁺¹
            a[i+1] = deepcopy(aⁿ)
            @logprogress i/nsteps
        end  
    end 
    
    return res, setup
end
