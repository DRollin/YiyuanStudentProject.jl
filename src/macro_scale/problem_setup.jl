function compute_Kₑ!(Kₑ::Matrix, cv::CellValues, problem::RVEProblem{dim}) where {dim}
    nbf  = Ferrite.getnbasefunctions(cv)
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        sens = compute_sensitivities(problem)
        for i in 1:nbf
            ∇δū = shape_gradient(cv, qp, i)
            for j in 1:nbf
                ū  = shape_value(cv, qp, j)
                ∇ū = shape_gradient(cv, qp, j)
                j̄ = get_j̄(sens; ū=ū, ∇ū=∇ū)
                Kₑ[i,j] += (-∇δū⋅j̄) * dΩ
            end
        end
    end
    return Kₑ
end

function compute_Mₑ!(Mₑ::Matrix, cv::CellValues, problem::RVEProblem{dim}) where {dim}
    nbf  = Ferrite.getnbasefunctions(cv)
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        sens = compute_sensitivities(problem)
        for i in 1:nbf
            δū  = shape_value(cv, qp, i)
            ∇δū = shape_gradient(cv, qp, i)
            for j in 1:nbf
                ū  = shape_value(cv, qp, j)
                ∇ū = shape_gradient(cv, qp, j)
                c̄  = get_c̄(sens; ū=ū, ∇ū=∇ū, l=true)
                c̄₂ = get_c̄₂(sens; ū=ū, ∇ū=∇ū, l=true)
                Mₑ[i,j] += (δū*c̄ + ∇δū⋅c̄₂) * dΩ
            end
        end
    end
    return Mₑ
end
#Mₑᴾᴾ, Mₑᴾᴹ, Mₑᴹᴾ, Mₑᴹᴹ

function assemble_K!(K::SparseMatrixCSC, dh::DofHandler, cv::CellValues, problem::RVEProblem{dim}) where {dim}
    assembler = start_assemble(K)
    nbf  = Ferrite.getnbasefunctions(cv)
    Kₑ   = zeros(nbf, nbf)
    
    for cc in CellIterator(dh)
        Ferrite.reinit!(cv, cc)
        fill!(Kₑ, 0)
        compute_Kₑ!(Kₑ, cv, problem)
        assemble!(assembler, celldofs(cc), Kₑ)
    end
    return K
end

function assemble_M!(M::SparseMatrixCSC, dh::DofHandler, cv::CellValues, problem::RVEProblem{dim}) where {dim}
    assembler = start_assemble(M)
    nbf  = Ferrite.getnbasefunctions(cv)
    Mₑ   = zeros(nbf, nbf)
    
    for cc in CellIterator(dh)
        Ferrite.reinit!(cv, cc)
        fill!(Mₑ, 0)
        compute_Mₑ!(Mₑ, cv, problem)
        assemble!(assembler, celldofs(cc), Mₑ)
    end
    return M
end