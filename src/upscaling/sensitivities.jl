function compute_effective_response(problem::RVEProblem{dim}, load::LoadCase{dim}) where {dim}
    u, setup = solve_load_case(problem, load)
    (; sets, setups) = setup
    c̄, c̄₂, j̄ = average_bulk_quantities(u, sets[:P], setups[:P]) .+ average_bulk_quantities(u, sets[:M], setups[:M])
    return EffectiveResponse{dim,Float64}(c̄, c̄₂, j̄)
end

function compute_sensitivities(problem::RVEProblem{dim}) where {dim}
    load  = LoadCase{dim}()
    l = compute_effective_response(problem, load)
    #@show l

    load  = LoadCase{dim}(μ̄=1.0)
    ū = compute_effective_response(problem, load) - l
    #@show ū

    load = LoadCase{dim}(ζ̄₁=1.0)
    ∇ū₁ = compute_effective_response(problem, load) - l

    load = LoadCase{dim}(ζ̄₂=1.0)
    ∇ū₂ = compute_effective_response(problem, load) - l
    
    ∇ū  =  (∇ū₁, ∇ū₂)

    res = RVEResponses(l, ū, ∇ū)
    return Sensitivities(res)

end

function compute_sensitivities_parallel(problem::RVEProblem{dim}) where {dim}
    load  = LoadCase{dim}()
    load_mu   = LoadCase{dim}(μ̄=1.0)
    load_zeta1 = LoadCase{dim}(ζ̄₁=1.0)
    load_zeta2 = LoadCase{dim}(ζ̄₂=1.0)


    EffR_pre =Any[]

    Load_Case = [load, load_mu, load_zeta1, load_zeta2]

    Threads.@threads for i in 1:4
        load_case_i = Load_Case[i]
        EffR_pre_i = compute_effective_response(problem, load_case_i)
        push!(EffR_pre, EffR_pre_i)
    end

    @show EffR_pre

    l = EffR_pre[1]
    ū = EffR_pre[2] - l
    ∇ū₁ = EffR_pre[3] - l
    ∇ū₂ = EffR_pre[4] - l
    ∇ū  =  (∇ū₁, ∇ū₂)

    res = RVEResponses(l, ū, ∇ū)
    return Sensitivities(res)

end

function get_macro_field(x̂::Union{SensitivitiesForScalar{dim,T},SensitivitiesForVector{dim,T}}; 
                         ū::Real=0.0, ∇ū::Tensor{1,dim}=zero(Tensor{1,dim}), l::Bool=false) where {dim,T}
    l ? x̄ = x̂.l : x̄ = zero(T)
    return x̄ .+ x̂.ū*ū .+ x̂.∇ū⋅∇ū
end

get_j̄(sens;  kwargs...) = get_macro_field(sens.j̄̂;  kwargs...)
get_c̄(sens;  kwargs...) = get_macro_field(sens.c̄̂;  kwargs...)
get_c̄₂(sens; kwargs...) = get_macro_field(sens.c̄̂₂; kwargs...)