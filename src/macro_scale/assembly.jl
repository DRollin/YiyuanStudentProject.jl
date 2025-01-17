
function assemble_macro_K!(setup::SolveSetup{dim}, Δt) where {dim}
    (; assemblysetup, K, aⁿ) = setup
    assembler = start_assemble(K)
    assemble_macro_K!(assembler, assemblysetup, aⁿ, Δt)
    return setup
end

function assemble_macro_K!(assembler, setup::AssemblySetup{dim}, aⁿ, Δt) where {dim}
    (; dh, cv, Kₑ, aₑ, gpdata) = setup
    @info "Assembling macro system"
    for cc in CellIterator(dh)
        reinit!(cv.u, cc)
        reinit!(cv.μ, cc)
        fill!(Kₑ, 0)
        aₑ .= aⁿ[celldofs(cc)]
        assemble_macro_element!(setup, Δt, gpdata[cellid(cc)])
        assemble!(assembler, celldofs(cc), Kₑ)
    end
    @info "Macro System assembled"
    return assembler
end

function assemble_macro_element!(setup::AssemblySetup{dim}, Δt, gpdata::Vector{GaussPointData{dim}}) where {dim}
    (; cv, nbf, Kₑ, subarrays, rvesetup) = setup
    (; Kₑuu, Kₑμμ, μₑ, uₑ) = subarrays
    
    for qp in 1:getnquadpoints(cv.u)
        dΩ = getdetJdV(cv.u, qp)

        μ̄ = function_value(cv.μ, qp, μₑ)
        ζ̄ = function_gradient(cv.μ, qp, μₑ)
        ε̄ = function_symmetric_gradient(cv.u, qp, uₑ)

        load = LoadCase{dim}(ε̄, μ̄, ζ̄)

        σ̄, ċ, ċ₂, j̄ = compute_effective_response!(gpdata[qp], rvesetup, load, Δt)

        for i in 1:nbf.u
            δNϵi = shape_gradient(cv.u, qp, i)
            for j in 1:nbf.u
                Kₑuu[i,j] += (δNϵi ⊡ σ̄  ) * dΩ
            end
        end

        for i in 1:nbf.μ
            δN∇μi = shape_gradient(cv.μ, qp, i)
            δNμi = shape_value(cv.μ, qp, i)
            for j in 1:nbf.μ
                Kₑμμ[i,j] += (δNμi * ċ - δN∇μi ⋅ (ċ₂ - j̄) ) * dΩ
            end
        end
    end
    return Kₑ
end

