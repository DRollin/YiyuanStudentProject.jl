"""
TODO

"""
function assemble_macro_K!(setup::SolveSetup, problem::RVESetup, Δt) where {dim}
    (; Assemblysetup, K, aⁿ) = setup
    (; data) = Assemblysetup
    assembler = start_assemble(K)
    
    assemble_macro_K!(assembler, Assemblysetup, problem, aⁿ, Δt)
    
    return K, data
end


function assemble_macro_K!(assembler, setup::AssemblySetup{dim}, problem::RVESetup, aⁿ, Δt) where {dim}
    (; dh, cv, Kₑ, aₑ) = setup
    @info "Assembling macro system"
    for cc in CellIterator(dh)
        reinit!(cv.u, cc)
        reinit!(cv.μ, cc)
        fill!(Kₑ, 0)
        aₑ .= aⁿ[celldofs(cc)]
        assemble_macro_element!(setup, problem, Δt, cellid(cc))
        assemble!(assembler, celldofs(cc), Kₑ)
    end
    @info "Macro system assembled"
    return assembler
end

"""
TODO

"""
function assemble_macro_element!(setup::AssemblySetup{dim}, problem::RVESetup, Δt, cellid) where {dim}
    (; dh, cv, nbf, Kₑ, subarrays, data, aₑ) = setup
    (; Kₑuu, Kₑμμ) = subarrays

    μₑ = @view(aₑ[dof_range(dh, :μ)])
    uₑ = @view(aₑ[dof_range(dh, :u)])
    
    for qp in 1:getnquadpoints(cv.u)
        dΩ = getdetJdV(cv.u, qp)

        μ̄ = function_value(cv.μ, qp, μₑ)
        ζ̄ = function_gradient(cv.μ, qp, μₑ)
        ε̄ = function_symmetric_gradient(cv.u, qp, uₑ)

        load = LoadCase{dim}(ε̄, μ̄, ζ̄)

       
        σ̄, c̄̇, c̄̇₂, j̄ = compute_effective_response(problem, load, data[cellid][qp], Δt)
        @show σ̄

        for i in 1:nbf.u
            δNϵi = shape_symmetric_gradient(cv.u, qp, i)
            for j in 1:nbf.u
                Kₑuu[i,j] += (δNϵi ⊡ σ̄ ) * dΩ
            end
        end

        for i in 1:nbf.μ
            δN∇μi = shape_gradient(cv.μ, qp, i)
            δNμi = shape_value(cv.μ, qp, i)

            for j in 1:nbf.μ
                Kₑμμ[i,j] += (δNμi * c̄̇   -  δN∇μi ⋅ (c̄̇₂ - j̄) ) * dΩ
            end
        end
    end
    return Kₑ
end

