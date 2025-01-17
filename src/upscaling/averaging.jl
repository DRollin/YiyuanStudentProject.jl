"""
TODO

"""
function average_quantities(a::Vector{Float64}, setup::RVESetup{dim}) where {dim}
    return average_quantities(a, setup.phasesetups.P) .+ average_quantities(a, setup.phasesetups.M)
end

function average_quantities(a::Vector{Float64}, setup::PhaseSetup{dim}) where {dim}
    (; dh, cells, cv, material) = setup
    (; E, αᶜʰ, k, cʳᵉᶠ, M, μʳᵉᶠ) = material

    σ̄  = zero(SymmetricTensor{2,dim})
    c̄  = 0.0
    c̄₂ = zero(Tensor{1,dim})
    j̄  = zero(Tensor{1,dim})

    for cc in CellIterator(dh, cells)
        reinit!(cv.u, cc)
        reinit!(cv.c, cc)
        reinit!(cv.μ, cc)
        coords = getcoordinates(cc)
        for qp in 1:getnquadpoints(cv.μ)
            dΩ = getdetJdV(cv.μ, qp)
            aᵉ = @view a[celldofs(cc)]
            μᵉ, uᵉ, cᵉ = aᵉ[dof_range(dh, :μ)], aᵉ[dof_range(dh, :u)], aᵉ[dof_range(dh, :c)]
            
            c  = function_value(cv.c, qp, cᵉ)
            ε  = function_symmetric_gradient(cv.u, qp, uᵉ)
            ∇μ = function_gradient(cv.μ, qp, μᵉ)
            x  = spatial_coordinate(cv.c, qp, coords)
            
            σ̄  += (E ⊡ (ε - αᶜʰ * (c - cʳᵉᶠ))) * dΩ
            c̄  += c * dΩ
            c̄₂ += c * x * dΩ
            j̄  -= M ⋅ ∇μ * dΩ
        end
    end
    return σ̄, c̄, c̄₂, j̄
end

"""
TODO

"""
function compute_effective_response!(gpdata::GaussPointData{dim}, rvesetup::RVESetup{dim}, load::LoadCase{dim}, Δt) where {dim}
    @info "Upscaling start"
    (; aⁿ, aⁿ⁺¹) = rvesetup
    (; c̄ⁿ, c̄₂ⁿ) = gpdata

    aⁿ .= gpdata.aᵣᵥₑⁿ      
    compute_time_step!(rvesetup, load, Δt)
    gpdata.aᵣᵥₑⁿ .= aⁿ⁺¹

    σ̄ⁿ⁺¹, c̄ⁿ⁺¹, c̄₂ⁿ⁺¹, j̄ⁿ⁺¹ = average_quantities(aⁿ⁺¹, rvesetup)
    ċ, ċ₂ = (c̄ⁿ⁺¹ - c̄ⁿ)/Δt, (c̄₂ⁿ⁺¹ - c̄₂ⁿ)/Δt

    gpdata.c̄ⁿ  = c̄ⁿ⁺¹
    gpdata.c̄₂ⁿ = c̄₂ⁿ⁺¹
    @info "End Upscaling"
    return σ̄ⁿ⁺¹, ċ, ċ₂, j̄ⁿ⁺¹
end