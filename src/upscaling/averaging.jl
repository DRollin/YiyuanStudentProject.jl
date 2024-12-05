function average_bulk_quantities(aₛₒₗ::Vector, setup::PhaseSetup)
    (; dh, cells, cv, material) = setup
    (; E, αᶜʰ, k, cʳᵉᶠ, M, μʳᵉᶠ) = material

    uₛₒₗ = evaluate_at_grid_nodes(dh, aₛₒₗ, :u)
    μₛₒₗ = evaluate_at_grid_nodes(dh, aₛₒₗ, :μ)
    cₛₒₗ = evaluate_at_grid_nodes(dh, aₛₒₗ, :μ)

    j̄  = zero(Tensor{1,dim})
    c̄ = 0.0
    c̄₂ = zero(Tensor{1,dim})
    σ̄  = zero(SymmetricTensor{2,2})

    for cc in CellIterator(dh, cells)
        reinit!(cv.u, cc)
        reinit!(cv.c, cc)
        reinit!(cv.μ, cc)
        coords = getcoordinates(cc)
        for qp in 1:getnquadpoints(cv.μ)
            dΩ = getdetJdV(cv.μ, qp)
            μ = function_value(cv.μ, qp, μₛₒₗ[celldofs(cc)])
            u = function_value(cv.u, qp, uₛₒₗ[celldofs(cc)])
            c = function_value(cv.u, qp, uₛₒₗ[celldofs(cc)])
            ϵ  = function_symmetric_gradient(cv.u, qp, uₛₒₗ[celldofs(cc)])
            ∇μ = function_gradient(cv.μ, qp, μₛₒₗ[celldofs(cc)])
            #x_μ  = spatial_coordinate(cv.μ, qp, coords)
            c̄  += (cʳᵉᶠ+(μ-μʳᵉᶠ)/k) * dΩ
            #c̄₂ += (cʳᵉᶠ+(μ-μʳᵉᶠ)/k)*x_μ * dΩ
            j̄  -= M⋅∇μ * dΩ
            σ̄  += E ⊡ (ϵ - αᶜʰ(c̄ - cʳᵉᶠ)) * dΩ
            
        end
    end
    return σ̄ , j̄
end