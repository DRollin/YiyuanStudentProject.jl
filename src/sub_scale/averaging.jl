function compute_homogenized_potential(uₛₒₗ::Vector, dh::DofHandler, setup::BulkElementSetup{dim}, cellset::Set{Int}) where {dim}
    (; cv) = setup
    ū = 0.0
    for cc in CellIterator(dh, cellset)
        Ferrite.reinit!(cv, cc)
        for qp in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, qp)
            u  = function_value(cv, qp, uₛₒₗ[celldofs(cc)])
            ū += u * dΩ
        end
    end
    return ū # Assuming |Ω|=1
end

function compute_homogenized_gradient(uₛₒₗ::Vector, dh::DofHandler, setup::BulkElementSetup{dim}, cellset::Set{Int}) where {dim}
    (; cv) = setup
    ∇ū = zero(Tensor{1,dim})
    for cc in CellIterator(dh, cellset)
        Ferrite.reinit!(cv, cc)
        for qp in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, qp)
            ∇u = function_gradient(cv, qp, uₛₒₗ[celldofs(cc)])
            ∇ū += ∇u * dΩ
        end
    end
    return ∇ū # Assuming |Ω|=1
end