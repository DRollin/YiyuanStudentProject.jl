"""
TODO

"""

function average_quantities(aₛₒₗ, phasedata::GaussPointPhaseData, setup::PhaseSetup{dim}, Δt) where {dim}
    (; c, c₂) = phasedata
    (; dh, cells, cv, material) = setup
    (; E, αᶜʰ, k, cʳᵉᶠ, M, μʳᵉᶠ) = material

    j  = zero(Tensor{1,dim})
    c_ₙ  = c
    c₂_ₙ = c₂
    c_ₙ₊₁  = 0.0
    c₂_ₙ₊₁ = zero(Tensor{1,dim})
    ċ = 0.0
    ċ₂ = zero(Tensor{1,dim})
    σ = zero(SymmetricTensor{2,dim})


    for cc in CellIterator(dh, cells)
        reinit!(cv.u, cc)
        reinit!(cv.c, cc)
        reinit!(cv.μ, cc)
    
        coords = getcoordinates(cc)
        for qp in 1:getnquadpoints(cv.μ)
            dΩ = getdetJdV(cv.μ, qp)
            aₛₒₗᵉ = @view aₛₒₗ[celldofs(cc)]
            μᵉ, uᵉ, cᵉ = aₛₒₗᵉ[dof_range(dh, :μ)], aₛₒₗᵉ[dof_range(dh, :u)], aₛₒₗᵉ[dof_range(dh, :c)]
            c = function_value(cv.c, qp, cᵉ)
            ϵ  = function_symmetric_gradient(cv.u, qp, uᵉ)
            ∇μ = function_gradient(cv.μ, qp, μᵉ)
            x  = spatial_coordinate(cv.c, qp, coords)
            c_ₙ₊₁  += c * dΩ
            c₂_ₙ₊₁ += c * x * dΩ
            j -= M ⋅ ∇μ * dΩ
            σ += (E ⊡ (ϵ - αᶜʰ * (c - cʳᵉᶠ))) * dΩ
        end
    end
    

    
    ċ = (c_ₙ₊₁ - c_ₙ) / Δt
    ċ₂ = (c₂_ₙ₊₁ .- c₂_ₙ) ./ Δt
    c = c_ₙ₊₁
    c₂ = c₂_ₙ₊₁

    return σ, ċ, ċ₂, j
end

"""
TODO

"""

function compute_effective_response(setup::RVESetup{dim}, load::LoadCase{dim}, data::GaussPointData{dim} , Δt) where {dim}
    @info "Upscaling start"

    
    setup.aⁿ .+= data.aⁿ_rve
              
    setup, aⁿ⁺¹ = compute_time_step!(setup, load, Δt)

    data.aⁿ_rve .= aⁿ⁺¹

    (; phasesetups) = setup
    σ̄ , c̄̇, c̄̇₂, j̄ = average_quantities(data.aⁿ_rve, data.phasedata[:P], phasesetups[:P], Δt) .+ average_quantities(data.aⁿ_rve, data.phasedata[:M], phasesetups[:M], Δt) 
   
    @info "End Upscaling"
    return σ̄, c̄̇, c̄̇₂, j̄, data
end