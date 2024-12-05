function average_quantities(aₛₒₗ::Vector, setup::PhaseSetup)
    (; dh, cells, cv, material) = setup
    (; E, αᶜʰ, k, cʳᵉᶠ, M, μʳᵉᶠ) = material

    uₛₒₗ = evaluate_at_grid_nodes(dh, aₛₒₗ, :u)
    μₛₒₗ = evaluate_at_grid_nodes(dh, aₛₒₗ, :μ)
    cₛₒₗ = evaluate_at_grid_nodes(dh, aₛₒₗ, :c)

    j̄  = zero(Tensor{1,dim})
    c̄  = 0.0
    c̄₂ = zero(Tensor{1,dim})
    σ̄  = zero(SymmetricTensor{2,2})

    for cc in CellIterator(dh, cells)
        reinit!(cv.u, cc)
        reinit!(cv.c, cc)
        reinit!(cv.μ, cc)
        coords = getcoordinates(cc)
        for qp in 1:getnquadpoints(cv.μ)
            dΩ = getdetJdV(cv.μ, qp)
            uₛₒₗᵉ = @view uₛₒₗ[celldofs(cc)]
            μᵉ, uᵉ, cᵉ = uₛₒₗᵉ[dofrange(dh, :μ)], uₛₒₗᵉ[dofrange(dh, :u)], uₛₒₗᵉ[dofrange(dh, :c)]
            μ = function_value(cv.μ, qp, μᵉ)
            u = function_value(cv.u, qp, uᵉ)
            c = function_value(cv.c, qp, cᵉ)
            ϵ  = function_symmetric_gradient(cv.u, qp, uᵉ)
            ∇μ = function_gradient(cv.μ, qp, μᵉ)
            x  = spatial_coordinate(cv.μ, qp, coords)
            c̄  += (cʳᵉᶠ+(μ-μʳᵉᶠ)/k)   * dΩ
            c̄₂ += (cʳᵉᶠ+(μ-μʳᵉᶠ)/k)*x * dΩ
            j̄  -= M ⋅ ∇μ * dΩ
            σ̄  += E ⊡ (ϵ - αᶜʰ(c̄ - cʳᵉᶠ)) * dΩ
        end
    end
    return σ̄, c̄, c̄₂, j̄
end
