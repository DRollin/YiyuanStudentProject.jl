"""
    average_quantities(a::Vector{Float64}, setup::RVESetup{dim}) where {dim}

Computes the average quantities `σ̄``, `c̄`, `c̄₂`, `j̄` for a given solution vector `a` over the phase setups of an RVE. 

# Arguments:
- `a`:  Initialied or updated solution vector of the corresponding RVE problem,
- `setup`: An object `RVESetup` for solving RVE problem

# Implementation Details:
Initialize the average quantities `σ̄``, `c̄`, `c̄₂`, `j̄`.

For all cells update the cell values for all the unknown fields.

For each quadrature point the element volume is computed. Furthermore for each base function in corresponding field,
evaluate the value of `c` concentration, `ε` strain, `∇μ` the gradient of chemical potantial and `x` the spatial coordinate of the current quadrature point.

Then the macro scale consistent macro scale fields are computed as following:


σ̄  = ∫ (E ⊡ (ε - αᶜʰ * (c - cʳᵉᶠ))) * dΩ / Vᵣᵥₑ

c̄ = ∫ c * dΩ / Vᵣᵥₑ

c̄₂ = ∫ (c * x) * dΩ / Vᵣᵥₑ

j̄ = - ∫ (M ⋅ ∇μ) * dΩ / Vᵣᵥₑ



"""
function average_quantities(a::Vector{Float64}, setup::RVESetup{dim}) where {dim}
    return (average_quantities(a, setup.phasesetups.P) .+ average_quantities(a, setup.phasesetups.M)) ./ setup.Vʳᵛᵉ
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
    compute_effective_response!(gpdata::GaussPointData{dim}, rvesetup::RVESetup{dim}, load::LoadCase{dim}, Δt) where {dim}

Computes the effective response of a material at the macroscale using an RVE. This function updates the state variables in the `GaussPointData` and performs upscaling calculations.

# Arguments:
- `gpdata`: An mutable object ``GaussPointData`` contains the state data,
- `rvesetup`: An object ``RVESetup`` helping for RVE solving,
- `load`: An object ``LoadCase`` as implicit boundary condition for RVE problem,
- `Δt`: The time increment for the current time step.


# Implementation Details:
Update the aᵣᵥₑⁿ from `gpdata` to the aⁿ in `rvesetup`.

Solve RVE problem at the current time step. Update the RVE solution in `gpdata`.

Computes upscaled quantities (stress, concentration, flux, etc.) using the `average_quantities` function.

Calculates the rates of change for concentration (`ċ`) and position-weighted concentration (`ċ₂`).
    
Updates `gpdata` with the newly computed averages.

"""
function compute_effective_response!(gpdata::GaussPointData{dim}, rvesetup::RVESetup{dim}, load::LoadCase{dim}, Δt) where {dim}
    @info "Upscaling start"
    (; aⁿ, aⁿ⁺¹) = rvesetup
    (; c̄ⁿ, c̄₂ⁿ) = gpdata

    for i in 1:1
        aⁿ .= gpdata.aᵣᵥₑⁿ      
        compute_time_step!(rvesetup, load, Δt)
        gpdata.aᵣᵥₑⁿ .= aⁿ⁺¹
        @info "RVE step $i"
    end
    

    σ̄ⁿ⁺¹, c̄ⁿ⁺¹, c̄₂ⁿ⁺¹, j̄ⁿ⁺¹ = average_quantities(aⁿ⁺¹, rvesetup)
    ċ, ċ₂ = (c̄ⁿ⁺¹ - c̄ⁿ)/Δt, (c̄₂ⁿ⁺¹ - c̄₂ⁿ)/Δt

    gpdata.c̄ⁿ  = c̄ⁿ⁺¹
    gpdata.c̄₂ⁿ = c̄₂ⁿ⁺¹
    @info "End Upscaling"
    return σ̄ⁿ⁺¹, ċ, ċ₂, j̄ⁿ⁺¹, aⁿ⁺¹
end