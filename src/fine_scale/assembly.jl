function compute_element_residual!(gₑ, uₑⁿ::AbstractVector, uₑⁿ⁺¹::AbstractVector, Δt::Real, ϵ::Real, dh::DofHandler, cvΦ::CellValues, cvμ::CellValues)
    dofsΦ, dofsμ = dof_range(dh, :Φ), dof_range(dh, :μ)
    	# Using views work with the different parts for the fields
	uΦⁿ   = @view uₑⁿ[dofsΦ]
    uΦⁿ⁺¹ = @view uₑⁿ⁺¹[dofsΦ]
    uμⁿ⁺¹ = @view uₑⁿ⁺¹[dofsμ]
    gₑΦ = @view gₑ[dofsΦ]
    gₑμ = @view gₑ[dofsμ]
    
    for qp in 1:getnquadpoints(cvΦ)
        dΩ    = getdetJdV(cvΦ, qp)
        Φⁿ    = function_value(cvΦ, qp, uΦⁿ)
        Φⁿ⁺¹  = function_value(cvΦ, qp, uΦⁿ⁺¹)
        ∇Φⁿ⁺¹ = function_gradient(cvΦ, qp, uΦⁿ⁺¹)
        μⁿ⁺¹  = function_value(cvμ, qp, uμⁿ⁺¹)
        ∇μⁿ⁺¹ = function_gradient(cvμ, qp, uμⁿ⁺¹)
        for i in 1:getnbasefunctions(cvμ)
            δμ  = shape_value(cvμ, qp, i)
            ∇δμ = shape_gradient(cvμ, qp, i)
            gₑμ[i] += ( δμ*(Φⁿ⁺¹-Φⁿ)/Δt + ϵ*∇δμ⋅∇μⁿ⁺¹ ) * dΩ
        end
        for i in 1:getnbasefunctions(cvΦ)
            δΦ  = shape_value(cvΦ, qp, i)
            ∇δΦ = shape_gradient(cvΦ, qp, i)
            gₑΦ[i] += ( δΦ*μⁿ⁺¹ - ϵ*∇δΦ⋅∇Φⁿ⁺¹ - δΦ*(Φⁿ⁺¹^3 - Φⁿ⁺¹)/ϵ ) * dΩ
        end
    end
    return gₑ
end

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
    Kₑ = zeros(setup.nbf, setup.nbf)
    for cc in CellIterator(setup.dh, cellset)
        Ferrite.reinit!(setup.cv, cc)
        fill!(Kₑ, 0)
        Kₑ = assemble_Kₑ!(Kₑ, setup)
        assemble!(assembler, celldofs(cc), Kₑ)
    end
    return assembler
end

function assemble_Kₑ!(Kₑ::Matrix, setup::iso_cm_ElementSetup)
    (; cells, cv, nbf, material) = setup
    (; E, α_ch, c_ref) = material

    Ke_u_u = @view Kₑ[dof_range(cells, :u), dof_range(cells, :u)]
    Ke_u_c = @view Kₑ[dof_range(cells, :u), dof_range(cells, :c)]

    
    for qp in 1:getnquadpoints(cv.u)
        for i in nbf.u
            δNϵi = shape_symmetric_gradient(cv.u, qp, i)
            
            for j in 1:nbf.u
                Nϵj = shape_symmetric_gradient(cv.u, qp, j)
                Ke_u_u[i, j]  += (δNϵi ⊡ E ⊡ Nϵj) * dΩ
            end	

            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Ke_u_c[i, j] -= (δNϵi ⊡ E ⊡ α_ch*(Ncj - c_ref)) * dΩ
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
    Mₑ = zeros(setup.nbf, setup.nbf)
    for cc in CellIterator(setup.dh, cellset)
        Ferrite.reinit!(setup.cv, cc)
        fill!(Mₑ, 0)
        Mₑ = assemble_Mₑ!(Mₑ, setup)
        assemble!(assembler, celldofs(cc), Mₑ)
    end
    return assembler
end

function assemble_Mₑ!(Mₑ::Matrix, setup::iso_cm_ElementSetup)
    (; cv, nbf, material) = setup
    (; E, α_ch, k, c_ref, M, μ_ref) = material

    Me_μ_μ = @view Mₑ[dof_range(cells, :μ), dof_range(cells, :μ)]
    Me_μ_c = @view Mₑ[dof_range(cells, :μ), dof_range(cells, :c)]

    for qp in 1:getnquadpoints(cv.μ)
        for i in nbf.μ
            δNμi = shape_value(cv.μ, qp, i)
            δN∇μi = shape_gradient(cv.μ, qp, i)

            δNci = shape_value(cv.c, qp, i)

            for j in nbf.μ
                N∇μj = shape_gradient(cv.μ, qp, j)
                Me_μ_μ[i, j] -= (δN∇μi ⋅ M ⋅ N∇μj) * dΩ
            end

            for j in nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Nμj = shape_value(cv.μ, qp, j)
                N∇uj = shape_gradient(cv.u, qp, j)
                #####?????
                δNci * Nμj .= δNci * (μ_ref + k*(Ncj-c_ref)-α_ch ⊡ E ⊡ (N∇uj - α_ch*(Ncj - c_ref)))
                #####?????
                Me_μ_c[i, j] += (Ncj * δNμi) * dΩ
            end
            
            ge = Me_μ_c
        end
    end
    
    
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
    f = zeros(ndofs(dh))
    
    f = adjust_f(f, Load)
    apply!(K, f, ch)
    apply!(M, f, ch)

    return K, M, f
end
