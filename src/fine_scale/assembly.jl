# IDEA: Merge assembly of K and M into one function to save some code?

"""
    TODO
"""
function assemble_K!(setup::RVESetup)
	(; phasesetups, K) = setup
    assembler = start_assemble(K)
    for phase in phasesetups
        assemble_K!(assembler, phase)
    end
    return K
end

function assemble_K!(assembler, setup::PhaseSetup)
    (; dh, cells, cv, Kₑ) = setup
    for cc in CellIterator(dh, cells)
        reinit!(cv.u, cc)
        reinit!(cv.c, cc)
        reinit!(cv.μ, cc)
        fill!(Kₑ, 0)
        assemble_Kₑ!(setup)
        assemble!(assembler, celldofs(cc), Kₑ)
    end
    return assembler
end

"""
    TODO
"""
function assemble_Kₑ!(setup::PhaseSetup)
    (; cv, nbf, material, Kₑ, submatrices) = setup
    (; E, αᶜʰ, k, cʳᵉᶠ, M, μʳᵉᶠ) = material
    (; Kₑuu, Kₑuc, Kₑcu, Kₑcc, Kₑcμ, Kₑμμ) = submatrices
    for qp in 1:getnquadpoints(cv.u)
        dΩ = getdetJdV(cv.u, qp)
        for i in 1:nbf.u
            δNϵi = shape_symmetric_gradient(cv.u, qp, i)
            
            for j in 1:nbf.u
                Nϵj = shape_symmetric_gradient(cv.u, qp, j)
                Kₑuu[i, j]  += (δNϵi ⊡ E ⊡ Nϵj) * dΩ
            end	

            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Kₑuc[i, j] -= (δNϵi ⊡ E ⊡ (αᶜʰ*(Ncj - cʳᵉᶠ))) * dΩ
            end
        end

        for i in 1:nbf.c
            δNci = shape_value(cv.c, qp, i)

            for j in 1:nbf.u
                Nϵj = shape_symmetric_gradient(cv.u, qp, j)
                Kₑcu[i, j] += (δNci * (αᶜʰ ⊡ E ⊡ Nϵj))* dΩ
            end

            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Kₑcc[i, j] -= (δNci * (k*(Ncj-cʳᵉᶠ) + αᶜʰ ⊡ E ⊡(αᶜʰ*(Ncj - cʳᵉᶠ))))* dΩ
            end

            for j in 1:nbf.μ
                Nμj = shape_value(cv.μ, qp, j)
                Kₑcμ[i, j] += (δNci * (Nμj - μʳᵉᶠ)) * dΩ
            end
        end

        for i in 1:nbf.μ
            δN∇μi = shape_gradient(cv.μ, qp, i)

            for j in 1:nbf.μ
                N∇μj = shape_gradient(cv.μ, qp, j)
                Kₑμμ[i, j] += (δN∇μi ⋅ M ⋅ N∇μj) * dΩ
            end
        end
    end
    
    return Kₑ
end

################################################################
# Assembly M
################################################################

"""
    TODO
"""
function assemble_M!(setup::RVESetup)
	(; phasesetups, M) = setup
    assembler = start_assemble(M)
    for phase in phasesetups
        assemble_M!(assembler, phase)
    end
    return M
end

function assemble_M!(assembler, setup::PhaseSetup)
    (; dh, cells, cv, Mₑ) = setup
    for cc in CellIterator(dh, cells)
        reinit!(cv.u, cc)
        reinit!(cv.c, cc)
        reinit!(cv.μ, cc)
        fill!(Mₑ, 0)
        assemble_Mₑ!(setup)
        assemble!(assembler, celldofs(cc), Mₑ)
    end
    return assembler
end

"""
    TODO
"""
function assemble_Mₑ!(setup::PhaseSetup)
    (; cv, nbf, Mₑ, submatrices) = setup
    (; Mₑμc) = submatrices
function assemble_Mₑ!(Mₑ::Matrix, setup::iso_cm_ElementSetup)
    (; cells, cv, nbf) = setup

    Me_μ_c = @view Mₑ[dof_range(cells, :μ), dof_range(cells, :c)]

    for qp in 1:getnquadpoints(cv.μ)
        dΩ = getdetJdV(cv.μ, qp)
        for i in 1:nbf.μ
        dΩ = getdetJdV(cv.μ, qp)
        for i in 1:nbf.μ
            δNμi = shape_value(cv.μ, qp, i)
            for j in 1:nbf.c
            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Mₑμc[i, j] += (Ncj * δNμi) * dΩ
                Me_μ_c[i, j] += (Ncj * δNμi) * dΩ
            end
        end
    end
    
    #@show Me_μ_c
    return Mₑ
end

################################################################
# Preparing system
################################################################
function adjust_f(f::Vector, load::LoadCase)
    (; μ̄) = load
    return vcat(f, [μ̄])
end

function assemble_rve_system(setup::FESetup_base{dim}) where {dim}
    (; dh, ch, Load) = setup

    K = assemble_K(setup)# TODO: Do not assemble everything anew every time
    M = assemble_M(setup)

    
    return K, M
end
