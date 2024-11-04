function doassemble_K_f!(problem::cm_Problem, Δt)
	(; setup, K, f, a, a_old) = problem
    (; sets, setups) = setup
    
    assembler = start_assemble(K,f)
    for (cellset, elementsetup) in zip(sets, setups)
        assemble_K_f!(assembler, cellset, elementsetup, a, a_old, Δt)
    end
    
end

function assemble_K_f!(assembler, cellset::Set{Int}, setup::iso_cm_ElementSetup, a, a_old, Δt)
    Ke = zeros(sum(setup.nbf), sum(setup.nbf))
    fe = zeros(sum(setup.nbf))
    ae = zeros(sum(setup.nbf))
    ae_old = zeros(sum(setup.nbf))
    
    for cc in CellIterator(setup.cells, cellset)
        reinit!(setup.cv.u, cc)
        reinit!(setup.cv.μ, cc)
        reinit!(setup.cv.c, cc)

        fill!(Ke, 0)
        fill!(fe, 0)

        ae .= a[celldofs(cc)]
        ae_old .= a_old[celldofs(cc)]

        assemble_Ke_fe!(Ke, fe, setup, ae, ae_old, Δt)
        assemble!(assembler, celldofs(cc), Ke, fe)
    end
    
end

function assemble_Ke_fe!(Ke::Matrix, fe::Vector, setup::iso_cm_ElementSetup, ae, ae_old, Δt)
    (; cells, cv, nbf, material) = setup
    (; E, α_ch, k, c_ref, M, μ_ref) = material


    ae_u = @view ae[dof_range(cells, :u)]
    ae_μ = @view ae[dof_range(cells, :μ)]
    ae_c = @view ae[dof_range(cells, :c)]

    ae_u_old = @view ae_old[dof_range(cells, :u)]
    ae_μ_old = @view ae_old[dof_range(cells, :μ)]
    ae_c_old = @view ae_old[dof_range(cells, :c)]

    #@show ae_u

    fe_u = @view fe[dof_range(cells, :u)]
    fe_μ = @view fe[dof_range(cells, :μ)]
    fe_c = @view fe[dof_range(cells, :c)]

    Ke_u_u = @view Ke[dof_range(cells, :u), dof_range(cells, :u)]
    Ke_u_μ = @view Ke[dof_range(cells, :u), dof_range(cells, :μ)]
    Ke_u_c = @view Ke[dof_range(cells, :u), dof_range(cells, :c)]

    Ke_μ_u = @view Ke[dof_range(cells, :μ), dof_range(cells, :u)]
    Ke_μ_μ = @view Ke[dof_range(cells, :μ), dof_range(cells, :μ)]
    Ke_μ_c = @view Ke[dof_range(cells, :μ), dof_range(cells, :c)]
    
    Ke_c_u = @view Ke[dof_range(cells, :c), dof_range(cells, :u)]
    Ke_c_μ = @view Ke[dof_range(cells, :c), dof_range(cells, :μ)]
    Ke_c_c = @view Ke[dof_range(cells, :c), dof_range(cells, :c)]
    

    
    for qp in 1:getnquadpoints(cv.u)

        dΩ = getdetJdV(cv.c, qp)

        #for c
        c = function_value(cv.c, qp, ae_c)
        @show c
        c_old = function_value(cv.c, qp, ae_c_old)
        ċ = (c - c_old) / Δt

        #for μ
        ∇μ = function_gradient(cv.μ, qp, ae_μ)

        #for u
        ϵ = function_symmetric_gradient(cv.u, qp, ae_u)
        
        for i in nbf.u
            δNϵi = shape_symmetric_gradient(cv.u, qp, i)
            
            for j in 1:nbf.u
                Nϵj = shape_symmetric_gradient(cv.u, qp, j)
                Ke_u_u[i, j]  += (δNϵi ⊡ E ⊡ Nϵj) * dΩ
            end	

            for j in 1:nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Ke_u_c[i, j] -= (δNϵi ⊡ E ⊡ α_ch*(Ncj - c_ref)) * dΩ
            end

            fe_u[i] += (δNϵi ⊡ E ⊡ ϵ - δNϵi ⊡ E ⊡ α_ch*(c - c_ref)) * dΩ
        end

        for i in nbf.μ
            δNμi = shape_value(cv.μ, qp, i)
            δN∇μi = shape_gradient(cv.μ, qp, i)

            for j in nbf.μ
                N∇μj = shape_gradient(cv.μ, qp, j)
                Ke_μ_μ[i, j] -= (δN∇μi ⋅ M ⋅ N∇μj) * dΩ
            end

            for j in nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Ke_μ_c[i, j] += (Ncj / Δt * δNμi) * dΩ
            end

            fe_μ[i] += (δNμi * ċ - δN∇μi ⋅ M ⋅ ∇μ) * dΩ
        end

        for i in nbf.c
            δNci = shape_value(cv.c, qp, i)

            for j in nbf.u
                Nϵj = shape_symmetric_gradient(cv.u, qp, j)
                Ke_c_u[i, j] += (α_ch ⊡ E ⊡ Nϵj) * dΩ
            end

            for j in nbf.μ
                Nμj = shape_value(cv.μ, qp, j)
                Ke_c_μ[i, j] += (δNci * Nμj) * dΩ
            end

            for j in nbf.c
                Ncj = shape_value(cv.c, qp, j)
                Ke_c_c[i, j] -= (δNci * (μ_ref + k * (Ncj - c_ref))) * dΩ
            end

        end
    end


end

