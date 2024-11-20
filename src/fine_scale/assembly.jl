function assemble_K(setup::FESetup_base)
	(; dh, ch, sets, setups) = setup
    K = allocate_matrix(dh, ch)
    assembler = start_assemble(K)
    for (cellset, elementsetup) in zip(sets, setups)
        assemble_K!(assembler, cellset, elementsetup)
    end
    return K
end

function assemble_K!(assembler, cellset::Set{Int}, setup::iso_cm_ElementSetup)
    Kₑ = zeros(sum(setup.nbf), sum(setup.nbf))
    for cc in CellIterator(setup.cells, cellset)
        Ferrite.reinit!(setup.cv.u, cc)
        Ferrite.reinit!(setup.cv.c, cc)
        Ferrite.reinit!(setup.cv.μ, cc)
        fill!(Kₑ, 0)
        Kₑ = assemble_Kₑ!(Kₑ, setup)
        assemble!(assembler, celldofs(cc), Kₑ)
    end
    return assembler
end

function assemble_Kₑ!(Kₑ::Matrix, setup::iso_cm_ElementSetup)
    (; cells, cv, nbf, material) = setup
    (; E, α_ch, k, c_ref, M, μ_ref) = material

    Ke_u_u = @view Kₑ[dof_range(cells, :u), dof_range(cells, :u)]
    Ke_u_c = @view Kₑ[dof_range(cells, :u), dof_range(cells, :c)]

    Ke_c_u = @view Kₑ[dof_range(cells, :c), dof_range(cells, :u)]
    Ke_c_c = @view Kₑ[dof_range(cells, :c), dof_range(cells, :c)]
    Ke_c_μ = @view Kₑ[dof_range(cells, :c), dof_range(cells, :μ)]

    Ke_μ_μ = @view Kₑ[dof_range(cells, :μ), dof_range(cells, :μ)]


    
    for qp in 1:getnquadpoints(cv.u)
        dΩ = getdetJdV(cv.u, qp)
        for i in 1:nbf.u
            δNϵi = shape_symmetric_gradient(cv.u, qp, i)
            
            for j in 1:nbf.u
                Nϵj = shape_symmetric_gradient(cv.u, qp, j)
                Ke_u_u[i, j]  += (δNϵi ⊡ E ⊡ Nϵj) * dΩ
            end	

            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Ke_u_c[i, j] -= (δNϵi ⊡ E ⊡ (α_ch*(Ncj - c_ref))) * dΩ
            end
        end

        for i in 1:nbf.c
            δNci = shape_value(cv.c, qp, i)

            for j in 1:nbf.u
                Nϵj = shape_symmetric_gradient(cv.u, qp, j)
                Ke_c_u[i, j] += (δNci * (α_ch ⊡ E ⊡ Nϵj))* dΩ
            end

            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Ke_c_c[i, j] -= (δNci * (k*(Ncj-c_ref) + α_ch ⊡ E ⊡(α_ch*(Ncj - c_ref))))* dΩ
            end

            for j in 1:nbf.μ
                Nμj = shape_value(cv.μ, qp, j)
                Ke_c_μ[i, j] += (δNci * (Nμj - μ_ref)) * dΩ
            end
        end

        for i in 1:nbf.μ
            δN∇μi = shape_gradient(cv.μ, qp, i)

            for j in 1:nbf.μ
                N∇μj = shape_gradient(cv.μ, qp, j)
                Ke_μ_μ[i, j] += (δN∇μi ⋅ M ⋅ N∇μj) * dΩ
            end
        end
    end
    
    return Kₑ
end

################################################################
# Assembly M
################################################################

function assemble_M(setup::FESetup_base)
	(; dh, ch, sets, setups) = setup
    M = allocate_matrix(dh, ch)
    assembler = start_assemble(M)
    for (cellset, elementsetup) in zip(sets, setups)
        assemble_M!(assembler, cellset, elementsetup)
    end
    return M
end

function assemble_M!(assembler, cellset::Set{Int}, setup::iso_cm_ElementSetup)
    Mₑ = zeros(sum(setup.nbf), sum(setup.nbf))
    for cc in CellIterator(setup.cells, cellset)
        Ferrite.reinit!(setup.cv.u, cc)
        Ferrite.reinit!(setup.cv.c, cc)
        Ferrite.reinit!(setup.cv.μ, cc)
        fill!(Mₑ, 0)
        Mₑ = assemble_Mₑ!(Mₑ, setup)
        assemble!(assembler, celldofs(cc), Mₑ)
    end
    return assembler
end

function assemble_Mₑ!(Mₑ::Matrix, setup::iso_cm_ElementSetup)
    (; cells, cv, nbf) = setup

    Me_μ_c = @view Mₑ[dof_range(cells, :μ), dof_range(cells, :c)]

    for qp in 1:getnquadpoints(cv.μ)
        dΩ = getdetJdV(cv.μ, qp)
        for i in 1:nbf.μ
            δNμi = shape_value(cv.μ, qp, i)
            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
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
