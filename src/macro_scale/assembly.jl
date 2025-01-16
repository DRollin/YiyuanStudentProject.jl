
function assemble_macro_K!(setup::SolveSetup, problem::RVESetup, Δt) where {dim}
    (; Assemblysetup, K) = setup
    (; data) = Assemblysetup
    assembler = start_assemble(K)
    
    assemble_macro_K!(assembler, Assemblysetup, problem, Δt)
    
    return K, data
end


function assemble_macro_K!(assembler, setup::AssemblySetup{dim}, problem::RVESetup, Δt) where {dim}
    (; dh, cv, Kₑ) = setup
    @info "Assembling macro system"
    for cc in CellIterator(dh)
        reinit!(cv.u, cc)
        reinit!(cv.μ, cc)
        fill!(Kₑ, 0)
        aₑ .= aₙ[celldofs(cc)]
        assemble_macro_element!(setup, problem, Δt, cellid(cc))
        assemble!(assembler, celldofs(cc), Kₑ)
    end
    @info "Macro System assembled"
    return assembler
end


function assemble_macro_element!(setup::AssemblySetup{dim}, problem::RVESetup, Δt, cellid) where {dim}
    (; cv, nbf, Kₑ, subarrays, data) = setup
    (; Kₑuu, Kₑμμ, μₑ, uₑ) = subarrays
    
    for qp in 1:getnquadpoints(cv.u)
        dΩ = getdetJdV(cv.u, qp)

        μ̄ = function_value(cv.μ, qp, μₑ)
        ζ̄ = function_gradient(cv.μ, qp, μₑ)
        ε̄ = function_symmetric_gradient(cv.u, qp, uₑ)

        load = LoadCase{dim}(ε̄, μ̄, ζ̄)

       
        σ̄, c̄̇, c̄̇₂, j̄, dataₙ₊₁ = compute_effective_response(problem, load, data[cellid][qp], Δt)


        for i in 1:nbf.u
            δNϵi = shape_gradient(cv.u, qp, i)
            for j in 1:nbf.u
                Nϵj = shape_value(cv.u, qp, j)  #??
                Kₑuu[i,j] += (δNϵi ⊡ σ̄  ) * dΩ
            end
        end

        for i in 1:nbf.μ
            δN∇μi = shape_gradient(cv.μ, qp, i)
            δNμi = shape_value(cv.μ, qp, i)

            for j in 1:nbf.μ
                N∇μj = shape_gradient(cv.μ, qp, j)#??
                Nμj = shape_value(cv.μ, qp, j)

                Kₑμμ[i,j] += (δNμi * c̄̇   -  δN∇μi ⋅ (c̄̇₂ - j̄) ) * dΩ
            end
        end
    end
    return Kₑ
end

