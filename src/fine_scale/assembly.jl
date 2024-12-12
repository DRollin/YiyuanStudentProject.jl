"""
        assemble_K_M!(setup::RVESetup)

Return the assembled mass matrix `M` and stiffness matrix `K`.

# Arguments:
- `setup`:   A preconfigured object type that defined for RVE assembly. Fields from `setup` are used in this function:
    - `phasesetups`: A preconfigured object type that defined for elementweise assembly,
    - `M`:           Initialized gloable mass matrix to be assembled,
    - `K`:           Initialized gloable stiffness matrix to be assembled.

# Implementation Details:
Create a `CSCAssembler` for both K and M.

Do the whole assembly for all the phases:
    for all cells update the cell values for all the unknown fields and initialize both `Kₑ` and `Mₑ`,
    do element assembly,
    return assembled `Kₑ` and `Mₑ`.

"""
function assemble_K_M!(setup::RVESetup)
	(; phasesetups, K, M) = setup
    assembler_K = start_assemble(K)
    assembler_M = start_assemble(M)
    for phase in phasesetups
        assemble_K_M!(assembler_K, assembler_M, phase)
    end
    return K, M
end

function assemble_K_M!(assembler_K, assembler_M, setup::PhaseSetup)
    (; dh, cells, cv, Kₑ, Mₑ) = setup
    for cc in CellIterator(dh, cells)
        reinit!(cv.u, cc)
        reinit!(cv.c, cc)
        reinit!(cv.μ, cc)
        fill!(Kₑ, 0)
        fill!(Mₑ, 0)
        assemble_element!(setup)
        assemble!(assembler_K, celldofs(cc), Kₑ)
        assemble!(assembler_M, celldofs(cc), Mₑ)
    end
    return assembler_K, assembler_M
end

"""
        assemble_element!(setup::PhaseSetup)

Return the assembled mass matrix for coupling between the chemical potantial and the concentration fields.

# Arguments:
- `setup`:  a preconfigured object type that defined for elementweise assembly. Fields from `setup` are used in this function:

    - `cv`:             Cell values for all unknown fields,
    - `nbf`:            The number of basis functions for all unknown fields,
    - `material`:       The predefined material parameters for both phase particals and matrix:
        - `E`:          fourth order stiffness tensor E,
        - `αᶜʰ`:        isotropic ion intercalation tensor,
        - `k`:          concentration-chemical potantial coefficient,
        - `cʳᵉᶠ`:       reference concentration,
        - `M`:          mobility tensor,
        - `μʳᵉᶠ`        reference chemical potantial.             
    - `Kₑ`:             The initialized element stiffness matrix to be assembled,
    - `Mₑ`:             The initialized element mass matrix to be assembled,
    - `submatrices`:    Submatrices of certain coupled unknown fields, including:
        - `Kₑuu`:   	    Coupling matrix between displacement and displacement fields.
        - `Kₑuc`:   	    Coupling matrix between displacement and concentration fields.
        - `Kₑcu`:   	    Coupling matrix between concentration and displacement fields.
        - `Kₑcc`:   	    Coupling matrix between concentration and concentration fields.
        - `Kₑcμ`:   	    Coupling matrix between concentration and chemical potential fields.
        - `Kₑμμ`:           Coupling matrix between chemical potential and chemical potential fields.
        - `Mₑμc`:   	    Coupling matrix between chemical potential and concentration fields.

# Implementation Details:
For each quadrature point the element volume is computed. Furthermore for each base function in corresponding field,
evaluate the shape function for test function. The nodal value of a certain unknown field is then computed for each base funtion in that field.

Then the coupling submatrices are computed as:

Kₑuu = ∫(δNϵi ⊡ E ⊡ Nϵj) * dΩ

Kₑuc = ∫-(δNϵi ⊡ E ⊡ (αᶜʰ*(Ncj - cʳᵉᶠ))) * dΩ

Kₑcu = ∫(δNci * (αᶜʰ ⊡ E ⊡ Nϵj))* dΩ

Kₑcc = ∫-(δNci * (k*(Ncj-cʳᵉᶠ) + αᶜʰ ⊡ E ⊡(αᶜʰ*(Ncj - cʳᵉᶠ))))* dΩ

Kₑcμ = ∫(δNci * (Nμj - μʳᵉᶠ)) * dΩ

Kₑμμ = ∫(δN∇μi ⋅ M ⋅ N∇μj) * dΩ

Mₑμc = ∫(Ncj * δNμi) * dΩ

where:
- `Nxj`:     Shape function for x field nodal values
- `δNxi`:    Shape function for x field test function
- `N∇xj`:    Gradient of Shape function for x field nodal values
- `δN∇xi`:   Gradient of Shape function for x field test function
- `dΩ`:      Determinant of the Jacobian times the quadrature weight (element volume)

"""
function assemble_element!(setup::PhaseSetup)
    (; cv, nbf, material, Kₑ, Mₑ, submatrices) = setup
    (; E, αᶜʰ, k, cʳᵉᶠ, M, μʳᵉᶠ) = material
    (; Kₑuu, Kₑuc, Kₑcu, Kₑcc, Kₑcμ, Kₑμμ, Mₑμc) = submatrices

    for qp in 1:getnquadpoints(cv.u)
        dΩ = getdetJdV(cv.u, qp)
        for i in 1:nbf.u
            δNϵi = shape_symmetric_gradient(cv.u, qp, i)
            
            for j in 1:nbf.u
                Nϵj = shape_symmetric_gradient(cv.u, qp, j)
                Kₑuu[i, j]  += (δNϵi ⊡ E ⊡ Nϵj) * dΩ
            end	

            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Kₑuc[i, j] -= (δNϵi ⊡ E ⊡ (αᶜʰ*(Ncj - cʳᵉᶠ))) * dΩ
            end
        end

        for i in 1:nbf.c
            δNci = shape_value(cv.c, qp, i)

            for j in 1:nbf.u
                Nϵj = shape_symmetric_gradient(cv.u, qp, j)
                Kₑcu[i, j] += (δNci * (αᶜʰ ⊡ E ⊡ Nϵj))* dΩ
            end

            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Kₑcc[i, j] -= (δNci * (k*(Ncj-cʳᵉᶠ) + αᶜʰ ⊡ E ⊡(αᶜʰ*(Ncj - cʳᵉᶠ))))* dΩ
            end

            for j in 1:nbf.μ
                Nμj = shape_value(cv.μ, qp, j)
                Kₑcμ[i, j] += (δNci * (Nμj - μʳᵉᶠ)) * dΩ
            end
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

    return Kₑ, Mₑ
end

