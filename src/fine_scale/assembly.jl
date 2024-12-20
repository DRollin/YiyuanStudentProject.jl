"""
        assemble_K_M_f!(setup::RVESetup)

Return the assembled mass matrix `M` and stiffness matrix `K`.

# Arguments:
- `setup`:   A preconfigured object that defined for RVE assembly. Fields from `setup` are used in this function:
    - `phasesetups`: A preconfigured object type that defined for elementweise assembly,
    - `M`:           Initialized global mass matrix to be assembled,
    - `K`:           Initialized global stiffness matrix to be assembled,
    - `f`:           Initialized global right hand side vector to be assembled.

# Implementation Details:
Create a `CSCAssembler` for both K, f, and M.

Do the whole assembly for all the phases:
    for all cells update the cell values for all the unknown fields and fill zeros to the 
    do element assembly,
    return assembled `Kₑ`, `Mₑ` and `fₑ`.

"""
function assemble_K_M_f!(setup::RVESetup)
	(; phasesetups, K, M, f) = setup
    assembler_Kf = start_assemble(K, f)
    assembler_M  = start_assemble(M)
    for phase in phasesetups
        assemble_K_M_f!(assembler_Kf, assembler_M, K, phase)
    end
    return K, M, f
end

function assemble_K_M_f!(assembler_Kf, assembler_M, K, setup::PhaseSetup{dim}) where {dim}
    (; dh, cells, cv, Kₑ, Mₑ, fₑ) = setup
    @info "Assembling system"
    for cc in CellIterator(dh, cells)
        reinit!(cv.u, cc)
        reinit!(cv.c, cc)
        reinit!(cv.μ, cc)
        fill!(Kₑ, 0)
        fill!(Mₑ, 0)
        fill!(fₑ, 0)
        assemble_element!(setup)
        assemble!(assembler_Kf, celldofs(cc), Kₑ, fₑ)
        assemble!(assembler_M, celldofs(cc), Mₑ)
    end
    @info "System assembled"
    return assembler_Kf, assembler_M
end

"""
        assemble_element!(setup::PhaseSetup)

Return the assembled mass matrix for coupling between the chemical potantial and the concentration fields.

# Arguments:
- `setup`:  a preconfigured object type that defined for elementweise assembly. 

# Implementation Details:
For each quadrature point the element volume is computed. Furthermore for each base function in corresponding field,
evaluate the shape function for the test function. The nodal value of a certain unknown field is then computed for each base funtion in that field.

Then the coupling subarrays are computed as:

Kₑuu = ∫(δNϵi ⊡ E ⊡ Nϵj) * dΩ

Kₑuc = ∫-(δNϵi ⊡ E ⊡ (αᶜʰ*(Ncj - cʳᵉᶠ))) * dΩ

Kₑcu = ∫(δNci * (αᶜʰ ⊡ E ⊡ Nϵj))* dΩ

Kₑcc = ∫-(δNci * (k*(Ncj-cʳᵉᶠ) + αᶜʰ ⊡ E ⊡(αᶜʰ*(Ncj - cʳᵉᶠ))))* dΩ

Kₑcμ = ∫(δNci * (Nμj - μʳᵉᶠ)) * dΩ

Kₑμμ = ∫(δN∇μi ⋅ M ⋅ N∇μj) * dΩ

Mₑμc = ∫(Ncj * δNμi) * dΩ

where:
- `Nxj`:     Shape function of field x 
- `δNxi`:    Shape function for x field test function
- `N∇xj`:    Gradient of shape function of field x 
- `δN∇xi`:   Gradient of shape function for test function of field x
- `dΩ`:      Determinant of the Jacobian times the quadrature weight
"""
function assemble_element!(setup::PhaseSetup)
    (; cv, nbf, material, Kₑ, Mₑ, fₑ, subarrays) = setup
    (; E, αᶜʰ, k, cʳᵉᶠ, M, μʳᵉᶠ) = material
    (; Kₑuu, Kₑuc, Kₑcu, Kₑcc, Kₑcμ, Kₑμμ, Mₑμc, fₑu, fₑc) = subarrays

   
    for qp in 1:getnquadpoints(cv.u)
        dΩ = getdetJdV(cv.u, qp)
        μ = function_value(cv.μ, qp, aₑμ)
        for i in 1:nbf.u
            δNϵi = shape_symmetric_gradient(cv.u, qp, i)
            
            for j in 1:nbf.u
                Nϵj = shape_symmetric_gradient(cv.u, qp, j)
                Kₑuu[i, j]  += (δNϵi ⊡ E ⊡ Nϵj) * dΩ
            end	

            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Kₑuc[i, j] -= (δNϵi ⊡ E ⊡ αᶜʰ * Ncj) * dΩ
            end

            fₑu[i] -= (δNϵi ⊡ E ⊡ αᶜʰ * cʳᵉᶠ) * dΩ
        end

        for i in 1:nbf.c
            δNci = shape_value(cv.c, qp, i)

            for j in 1:nbf.u
                Nϵj = shape_symmetric_gradient(cv.u, qp, j)
                Kₑcu[i, j] += (δNci * (Nϵj ⊡ E ⊡ αᶜʰ)) * dΩ
            end

            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Kₑcc[i, j] -= (δNci * (k + αᶜʰ ⊡ E ⊡ αᶜʰ) * Ncj) * dΩ
            end

            for j in 1:nbf.μ
                Nμj = shape_value(cv.μ, qp, j)
                Kₑcμ[i, j] += (δNci * Nμj) * dΩ
            end

            fₑc[i] += (δNci * (μʳᵉᶠ - (k + αᶜʰ ⊡ E ⊡ αᶜʰ)*cʳᵉᶠ)) * dΩ
        end

        for i in 1:nbf.μ
            δN∇μi = shape_gradient(cv.μ, qp, i)
            δNμi = shape_value(cv.μ, qp, i)

            for j in 1:nbf.μ
                N∇μj = shape_gradient(cv.μ, qp, j)
                Kₑμμ[i, j] += (δN∇μi ⋅ M ⋅ N∇μj) * dΩ
            end

            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Mₑμc[i, j] += (Ncj * δNμi) * dΩ
            end
        end
    end

    return Kₑ, Mₑ, fₑ
end
