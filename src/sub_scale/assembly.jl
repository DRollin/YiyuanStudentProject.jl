"""
    assemble_K_M_f!(setup::RVESetup)

Return the assembled stiffness matrix `K`, mass matrix `M` and right hand side vector `f`.

The assembly of all cells across all material phases is performed using the given input `RVESetup`.
"""
function assemble_K_M_f!(setup::RVESetup)
	(; phasesetups, K, M, f) = setup
    assembler_Kf = start_assemble(K, f)
    assembler_M  = start_assemble(M)
    for phase in phasesetups
        assemble_K_M_f!(assembler_Kf, assembler_M, phase)
    end
    return K, M, f
end

function assemble_K_M_f!(assembler_Kf, assembler_M, setup::PhaseSetup{dim}) where {dim}
    (; dh, cells, cv, Kₑ, Mₑ, fₑ) = setup
    @info "Assembling RVE system"
    for cc in CellIterator(dh, cells)
        reinit!(cv.u, cc)
        reinit!(cv.c, cc)
        reinit!(cv.μ, cc)
        fill!(Kₑ, 0)
        fill!(Mₑ, 0)
        fill!(fₑ, 0)
        assemble_element!(setup)
        assemble!(assembler_Kf, celldofs(cc), Kₑ, fₑ)
        assemble!(assembler_M, celldofs(cc), Mₑ)
    end
    @info "RVE System assembled"
    return assembler_Kf, assembler_M
end

"""
    assemble_element!(setup::PhaseSetup)

Return the assembled local stiffness and mass matrices (`Kₑ`, `Mₑ`), the local right hand side vector (`fₑ`).  

Different contributions of `Kₑ`, `Mₑ` and `fₑ` are computed with stored material parameters and necessary assembly setups provided by the input `PhaseSetup`.
"""
function assemble_element!(setup::PhaseSetup)
    (; cv, nbf, material, Kₑ, Mₑ, fₑ, subarrays) = setup
    (; E, αᶜʰ, k, cʳᵉᶠ, M, μʳᵉᶠ) = material
    (; Kₑuu, Kₑuc, Kₑcu, Kₑcc, Kₑcμ, Kₑμμ, Mₑμc, fₑu, fₑc) = subarrays

   
    for qp in 1:getnquadpoints(cv.u)
        dΩ = getdetJdV(cv.u, qp)
        #μ = function_value(cv.μ, qp, aₑμ)
        for i in 1:nbf.u
            δNϵi = shape_symmetric_gradient(cv.u, qp, i)
            
            for j in 1:nbf.u
                Nϵj = shape_symmetric_gradient(cv.u, qp, j)
                Kₑuu[i, j]  += (δNϵi ⊡ E ⊡ Nϵj) * dΩ
            end	

            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Kₑuc[i, j] -= (δNϵi ⊡ E ⊡ αᶜʰ * Ncj) * dΩ
            end

            fₑu[i] -= (δNϵi ⊡ E ⊡ αᶜʰ * cʳᵉᶠ) * dΩ
        end

        for i in 1:nbf.c
            δNci = shape_value(cv.c, qp, i)

            for j in 1:nbf.u
                Nϵj = shape_symmetric_gradient(cv.u, qp, j)
                Kₑcu[i, j] += (δNci * (Nϵj ⊡ E ⊡ αᶜʰ)) * dΩ
            end

            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Kₑcc[i, j] -= (δNci * (k + αᶜʰ ⊡ E ⊡ αᶜʰ) * Ncj) * dΩ
            end

            for j in 1:nbf.μ
                Nμj = shape_value(cv.μ, qp, j)
                Kₑcμ[i, j] += (δNci * Nμj) * dΩ
            end

            fₑc[i] += (δNci * (μʳᵉᶠ - (k + αᶜʰ ⊡ E ⊡ αᶜʰ)*cʳᵉᶠ)) * dΩ
        end

        for i in 1:nbf.μ
            δN∇μi = shape_gradient(cv.μ, qp, i)
            δNμi = shape_value(cv.μ, qp, i)

            for j in 1:nbf.μ
                N∇μj = shape_gradient(cv.μ, qp, j)
                Kₑμμ[i, j] += (δN∇μi ⋅ M ⋅ N∇μj) * dΩ
            end

            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Mₑμc[i, j] += (Ncj * δNμi) * dΩ
            end
        end
    end

    return Kₑ, Mₑ, fₑ
end
