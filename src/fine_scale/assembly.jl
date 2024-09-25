################################################################
# Assembly M
################################################################

function assemble_M(setup::FESetup)
	(; dh, ch, cellsets, elementsetups) = setup
    M = create_sparsity_pattern(dh, ch)
    assembler = start_assemble(M)
    for (cellset, elementsetup) in zip(cellsets, elementsetups)
        assemble_M!(assembler, cellset, elementsetup)
    end
    return M
end

function assemble_M!(assembler, cellset::Set{Int}, setup::BulkElementSetup)
    Mₑ = zeros(setup.nbf, setup.nbf)
    for cc in CellIterator(setup.dh, cellset)
        Ferrite.reinit!(setup.cv, cc)
        fill!(Mₑ, 0)
        Mₑ = assemble_Mₑ!(Mₑ, setup)
        assemble!(assembler, celldofs(cc), Mₑ)
    end
    return assembler
end

function assemble_Mₑ!(Mₑ::Matrix, setup::BulkElementSetup)
    (; cv, nbf, material) = setup
    (; k) = material
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        for i in 1:nbf
            δu = shape_value(cv, qp, i)
            for j in 1:nbf
                u = shape_value(cv, qp, j)
                Mₑ[i, j]  += (u*δu)/k * dΩ
            end	
        end
    end
    return Mₑ
end


################################################################
# Assembly K
################################################################

function assemble_K(setup::FESetup)
	(; dh, ch, sets, setups) = setup
    K = create_sparsity_pattern(dh, ch)
    assembler = start_assemble(K)
    for (cellset, elementsetup) in zip(sets, setups)
        assemble_K!(assembler, cellset, elementsetup)
    end
    return K
end

function assemble_K!(assembler, cellset::Set{Int}, setup::BulkElementSetup)
    Kₑ = zeros(setup.nbf, setup.nbf)
    for cc in CellIterator(setup.cells, cellset)
        Ferrite.reinit!(setup.cv, cc)
        fill!(Kₑ, 0)
        Kₑ = assemble_Kₑ!(Kₑ, setup)
        assemble!(assembler, celldofs(cc), Kₑ)
    end
    return assembler
end

function assemble_Kₑ!(Kₑ::Matrix, setup::BulkElementSetup)
    (; cv, nbf, material) = setup
    (; η) = material
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        for i in 1:nbf
            ∇δu = shape_gradient(cv, qp, i)
            #∇δu = shape_gradient(cv, qp, i)
            for j in 1:nbf
                ∇u = shape_gradient(cv, qp, j)
                Kₑ[i, j]  += ((η ⋅ ∇u) ⋅ ∇δu) * dΩ
            end	
        end
    end
    return Kₑ
end

################################################################
# Assembly C
################################################################

function assemble_C(setup::FESetup{dim}) where {dim}
	(; dh, sets, setups) = setup
    C = spzeros(1, ndofs(dh))
    assemble_C!(C, sets.P, setups.P)
    assemble_C!(C, sets.M, setups.M)
    return C
end

function assemble_C!(C::SparseMatrixCSC, cellset::Set{Int}, setup::BulkElementSetup)
    Cₑ = zeros(1, setup.nbf)
    for cc in CellIterator(setup.cells, cellset)
        Ferrite.reinit!(setup.cv, cc)
        fill!(Cₑ, 0)
        Cₑ = assemble_Cₑ!(Cₑ, setup)
        C[: , celldofs(cc)] .+= Cₑ
    end
    return C
end

function assemble_Cₑ!(Cₑ::Matrix, setup::BulkElementSetup)
    (; cv, nbf) = setup
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        for j in 1:nbf
            u  = shape_value(cv, qp, j)
            Cₑ[1, j]  +=  u * dΩ
        end
    end
    return Cₑ
end

################################################################
# Preparing system
################################################################

function adjust_K(K::SparseMatrixCSC, C::SparseMatrixCSC)
    return hcat( vcat(K,C), copy(hcat(C, spzeros(1,1))') )
end

function adjust_f(f::Vector, load::LoadCase)
    (; μ̄) = load
    return vcat(f, [μ̄])
end

function assemble_rve_system(setup::FESetup{dim}) where {dim}
    (; dh, ch, Load) = setup

    K = assemble_K(setup)# TODO: Do not assemble everything anew every time
    C = assemble_C(setup)
    f = zeros(ndofs(dh))

    K = adjust_K(K, C)
    f = adjust_f(f, Load)
    apply!(K, f, ch)

    return K, f
end

function assemble_dns_system(setup::FESetup{dim}) where {dim}
    (; dh) = setup
    M = assemble_M(setup)
    K = assemble_K(setup)
    f = zeros(ndofs(dh))
    return M, K, f
end