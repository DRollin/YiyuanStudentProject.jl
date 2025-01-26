"""
        assemble_macro_K!(setup::SolveSetup{dim}, Δt) where {dim}
    
Return the whole ``SolveSetup`` with assembled stiffness matrix ``K``.

# Arguments:
- `setup`:   ``SolveSetup``: A preconfigured object that defined for macro scale problem assembly. Fields from `setup` are used in this function:
    - `assemblysetup`: ``AssemblySetup``: A preconfigured object type that defined for elementweise assembly,
    - `K`:           Initialized global stiffness matrix to be assembled,
    - `aⁿ`:           Initialized global solution vector for calling the element solution vector


# Implementation Details:
Create a `CSCAssembler` for K.

Do the whole assembly:
    for all cells update the cell values for all the unknown fields and fill zeros to the 
    do element assembly, aliening the element solution vector with the global solution vector.
    return assembled `Kₑ`.

"""
function assemble_macro_K!(setup::SolveSetup{dim}, Δt) where {dim}
    (; assemblysetup, K, aⁿ) = setup
    assembler = start_assemble(K)
    assemble_macro_K!(assembler, assemblysetup, aⁿ, Δt)
    return setup
end

function assemble_macro_K!(assembler, setup::AssemblySetup{dim}, aⁿ, Δt) where {dim}
    (; dh, cv, Kₑ, aₑ, gpdata) = setup
    @info "Assembling macro system"
    for cc in CellIterator(dh)
        reinit!(cv.u, cc)
        reinit!(cv.μ, cc)
        fill!(Kₑ, 0)
        aₑ .= aⁿ[celldofs(cc)]
        assemble_macro_element!(setup, Δt, gpdata[cellid(cc)])
        assemble!(assembler, celldofs(cc), Kₑ)
    end
    @info "Macro system assembled"
    return assembler
end

"""
        assemble_macro_element!(setup::AssemblySetup{dim}, Δt, gpdata::Vector{GaussPointData{dim}}) where {dim}

Return the assembled stiffness matrix.

# Arguments:
- `setup`:  a preconfigured object type that defined for elementweise assembly. 

# Implementation Details:
The local unknowns μₑ and uₑ are passed for creating the step dependent macro scale variables for each quadrature point. 
Object ``LoadCase`` is generated using the macro scale variables `μ̄ `, `ζ̄ `and `ε̄ `.
For each quadrature point the element volume is generated. The variationally consistent macro-scale (homogenized) fields σ̄, ċ, ċ₂, j̄ are computed passing the upscaling function ``compute_effective_response!``.
Furthermore for each base function in corresponding field, evaluate the shape function for the test function. 

Then the coupling subarrays are computed as:

Kₑuu = ∫(δNϵi ⊡ σ̄ ) * dΩ

Kₑμμ = ∫(δNμi * ċ - δN∇μi ⋅ (ċ₂ - j̄) ) * dΩ


where:
- `δNxi`:    Shape function for x field test function
- `δN∇xi`:   Gradient of shape function for test function of field x
- `dΩ`:      Determinant of the Jacobian times the quadrature weight
"""
function assemble_macro_element!(setup::AssemblySetup{dim}, Δt, gpdata::Vector{GaussPointData{dim}}) where {dim}
    (; dh, cv, nbf, Kₑ, subarrays, rvesetup, aₑ) = setup
    (; Kₑuu, Kₑμμ) = subarrays

    μₑ = @view(aₑ[dof_range(dh, :μ)])
    uₑ = @view(aₑ[dof_range(dh, :u)])

    for qp in 1:getnquadpoints(cv.u)
        dΩ = getdetJdV(cv.u, qp)

        μ̄ = function_value(cv.μ, qp, μₑ)
        ζ̄ = function_gradient(cv.μ, qp, μₑ)
        ε̄ = function_symmetric_gradient(cv.u, qp, uₑ)

        load = LoadCase{dim}(ε̄, μ̄, ζ̄)

        σ̄, ċ, ċ₂, j̄,  = compute_effective_response!(gpdata[qp], rvesetup, load, Δt)

        for i in 1:nbf.u
            δNϵi = shape_symmetric_gradient(cv.u, qp, i)
            for j in 1:nbf.u
                Kₑuu[i,j] += (δNϵi ⊡ σ̄  ) * dΩ
            end
        end

        for i in 1:nbf.μ
            δN∇μi = shape_gradient(cv.μ, qp, i)
            δNμi = shape_value(cv.μ, qp, i)
            for j in 1:nbf.μ
                Kₑμμ[i,j] += (δNμi * ċ - δN∇μi ⋅ (ċ₂ - j̄) ) * dΩ
            end
        end
    end

    return Kₑ
end

