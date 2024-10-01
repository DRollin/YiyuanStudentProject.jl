function doassemble_K_f!(problem::pv_Problem, Δt)
	(; setup, K, f, a, a_old, assembler) = problem
    (; sets, setups) = setup
    

    for (cellset, elementsetup) in zip(sets, setups)
        assemble_K_f!(assembler, cellset, elementsetup, a, a_old, Δt)
    end
    
end

function assemble_K_f!(assembler, cellset::Set{Int}, setup::iso_pv_ElementSetup, a, a_old, Δt)
    Ke = zeros(setup.nbf, setup.nbf)
    fe = zeros(setup.nbf)
    ae = zeros(setup.nbf)
    ae_old = zeros(setup.nbf)
    
    for cc in CellIterator(setup.cells, cellset)
        map!(i->a[i], ae, celldofs(cc))
        map!(i->a_old[i], ae_old, celldofs(cc))

        fill!(Ke, 0)
        fill!(fe, 0)

        reinit!.(setup.cv, cc)
        assemble_Ke_fe!(Ke, fe, setup, cc, ae, ae_old, Δt)
        assemble!(assembler, celldofs(cc), Ke, fe)
    end
    
end

function assemble_Ke_fe!(Ke::Matrix, fe::Vector, setup::iso_pv_ElementSetup, cell, ae, ae_old, Δt)
    (; cv, nbf, material) = setup
    (; E, Κ, α, β) = material
    ae_u = @view ae[dof_range(cell, :μ)]

    for qp in 1:getnquadpoints(cv.u)
        dΩ = getdetJdV(cv.u, qp)
        p = function_value(cv.p, qp, a, dof_range(dh, :p))

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

    K = assemble_K(setup)
    C = assemble_C(setup)
    f = zeros(ndofs(dh))

    K = adjust_K(K, C)
    f = adjust_f(f, Load)
    apply!(K, f, ch)

    return K, f
end