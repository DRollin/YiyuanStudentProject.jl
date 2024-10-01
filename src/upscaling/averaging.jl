function average_bulk_quantities(uₛₒₗ::Vector, cellset::Set{Int}, setup::BulkElementSetup{dim}) where {dim}
    (; cells, cv, material) = setup
    (; η, k, μᵣₑ, cᵣₑ) = material
    j̄  = zero(Tensor{1,dim})
    c̄ = 0.0
    c̄₂ = zero(Tensor{1,dim})
    for cc in CellIterator(cells, cellset)
        Ferrite.reinit!(cv, cc)
        coords = getcoordinates(cc)
        for qp in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, qp)
            u  = function_value(   cv, qp, uₛₒₗ[celldofs(cc)])
            ∇u = function_gradient(cv, qp, uₛₒₗ[celldofs(cc)])
            x  = spatial_coordinate(cv, qp, coords)
            c̄  += (cᵣₑ+(u-μᵣₑ)/k) * dΩ
            c̄₂ += (cᵣₑ+(u-μᵣₑ)/k)*x * dΩ
            j̄  -= η⋅∇u * dΩ
        end
    end
    return c̄, c̄₂, j̄
end