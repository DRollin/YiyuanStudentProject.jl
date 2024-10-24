function doassemble_K_f!(problem::pv_Problem, Δt)
	(; setup, K, f, a, a_old) = problem
    (; sets, setups) = setup
    
    assembler = start_assemble(K,f)
    for (cellset, elementsetup) in zip(sets, setups)
        assemble_K_f!(assembler, cellset, elementsetup, a, a_old, Δt)
    end
    
end

function assemble_K_f!(assembler, cellset::Set{Int}, setup::iso_pv_ElementSetup, a, a_old, Δt)
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

function assemble_Ke_fe!(Ke::Matrix, fe::Vector, setup::iso_pv_ElementSetup, ae, ae_old, Δt)
    (; cells, cv, nbf, material) = setup
    (; E, Κ, α, β) = material

    ae_u = @view ae[dof_range(cells, :u)]
    ae_p = @view ae[dof_range(cells, :μ)]
    ae_u_old = @view ae_old[dof_range(cells, :u)]
    ae_p_old = @view ae_old[dof_range(cells, :μ)]

    #@show ae_u

    fe_u = @view fe[dof_range(cells, :u)]
fe_p = @view fe[dof_range(cells, :μ)]
    Ke_u_u = @view Ke[dof_range(cells, :u), dof_range(cells, :u)]
    Ke_u_p = @view Ke[dof_range(cells, :u), dof_range(cells, :μ)]
    Ke_p_u = @view Ke[dof_range(cells, :μ), dof_range(cells, :u)]
    Ke_p_p = @view Ke[dof_range(cells, :μ), dof_range(cells, :μ)]


    for qp in 1:getnquadpoints(cv.u)
        dΩ = getdetJdV(cv.u, qp)
        #@show dΩ
        #for p
        p = function_value(cv.μ, qp, ae_μ)
        #@show p
        p_old = function_value(cv.μ, qp, ae_μ_old)
        ṗ = (p - p_old) / Δt
        ζ = function_gradient(cv.μ, qp, ae_μ)
        #@show ζ
        #for u
        ϵ = function_symmetric_gradient(cv.u, qp, ae_u)
        #@show ϵ
        div_u = tr(ϵ)
        div_u_old = function_divergence(cv.u, qp, ae_u_old)
        div_u̇ = (div_u - div_u_old) / Δt
        σ = E ⊡ ϵ
        
        for i in nbf.u
            δNϵi = shape_symmetric_gradient(cv.u, qp, i)
            div_δNui = shape_divergence(cv.u, qp, i)
            
            for j in 1:nbf.u
                Nϵj = shape_symmetric_gradient(cv.u, qp, j)
                Ke_u_u[i, j]  += (δNϵi ⊡ E ⊡ Nϵj) * dΩ
            end	

            for j in 1:nbf.p
                Npj_u = shape_value(cv.p, qp, j)
                Ke_u_p[i, j] -= (div_δNui * α * Npj_u) * dΩ
            end

            fe_u[i] += (δNϵi ⊡ σ - div_δNui*α*p) * dΩ
        end

        for i in nbf.p
            δNpi = shape_value(cv.p, qp, i)
            δNζi = shape_gradient(cv.p, qp, i)

            for j in nbf.p
                Npj_p = shape_value(cv.p, qp, j)
                Nζj = shape_gradient(cv.p, qp, j)
                Ke_p_p[i, j] += (δNpi * β * Npj_p/Δt + δNζi ⋅ Κ ⋅ Nζj) * dΩ
            end

            for j in nbf.u
                div_Nuj = shape_divergence(cv.u, qp, j)
                Ke_p_u[i, j] += (α * δNpi * div_Nuj) * dΩ
            end

            fe_p[i] += (δNpi * (α * div_u̇ + β * ṗ) + δNζi ⋅ Κ ⋅ ζ) * dΩ
        end
    end


end

