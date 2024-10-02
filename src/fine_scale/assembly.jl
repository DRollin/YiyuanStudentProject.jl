function doassemble_K_f!(problem::pv_Problem, Δt)
	(; setup, K, f, a, a_old) = problem
    (; sets, setups) = setup
    
    assembler = start_assemble(K,f)
    for (cellset, elementsetup) in zip(sets, setups)
        assemble_K_f!(assembler, cellset, elementsetup, a, a_old, Δt)
    end
    
end

function assemble_K_f!(assembler, cellset::Set{Int}, setup::iso_pv_ElementSetup, a, a_old, Δt)
    Ke = zeros(setup.nbf, setup.nbf)
    fe = zeros(setup.nbf)
    ae = zeros(setup.nbf)
    ae_old = zeros(setup.nbf)
    
    for cc in CellIterator(setup.cells, cellset)
        map!(i->a[i], ae, celldofs(cc))
        map!(i->a_old[i], ae_old, celldofs(cc))

        fill!(Ke, 0)
        fill!(fe, 0)

        reinit!.(setup.cv, cc)
        assemble_Ke_fe!(Ke, fe, setup, ae, ae_old, Δt)
        assemble!(assembler, celldofs(cc), Ke, fe)
    end
    
end

function assemble_Ke_fe!(Ke::Matrix, fe::Vector, setup::iso_pv_ElementSetup, ae, ae_old, Δt)
    (; cells, cv, nbf, material) = setup
    (; E, Κ, α, β) = material

    ae_u = @view ae[dof_range(cells, :u)]
    ae_p = @view ae[dof_range(cells, :p)]
    fe_u = @view fe[dof_range(cells, :u)]
    fe_p = @view fe[dof_range(cells, :p)]
    Ke_u_u = @view Ke[[dof_range(cells, :u)], [dof_range(cells, :u)]]
    Ke_u_p = @view Ke[[dof_range(cells, :u)], [dof_range(cells, :p)]]
    Ke_p_u = @view Ke[[dof_range(cells, :p)], [dof_range(cells, :u)]]
    Ke_p_p = @view Ke[[dof_range(cells, :p)], [dof_range(cells, :p)]]


    for qp in 1:getnquadpoints(cv.u)
        dΩ = getdetJdV(cv.u, qp)
        #for p
        p = function_value(cv.p, qp, ae, ae_p)
        p_old = function_value(cv.p, qp, ae_old, ae_p)
        ṗ = (p - p_old) / Δt
        ζ = function_gradient(cv.p, qp, ae_p)
        #for u
        ϵ = function_symmetric_gradient(cv.u, qp, ae, ae_u)
        div_u = tr(ϵ)
        div_u_old = function_divergence(cv.u, qp, ae_old, ae_u)
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

            fe[i] += (δNϵi ⊡ σ - div_δNui*α*p) * dΩ
        end

        for i in nbf.p
            δNpi = shape_value(cv.p, qp, i)
            δNζi = shape_gradient(cv.p, qp, i)

            for j in nbf.p
                Npj_p = shape_value(cv.p, qp, j)
                Nζi = shape_gradient(cv.p, qp, j)
                Ke_p_p[i, j] += (β * δNpi * Npj_p/Δt + δNζi * Κ * Nζi) * dΩ
            end

            for j in nbf.u
                div_Nuj = shape_divergence(cv.u, qp, j)
                Ke_p_u[i, j] += (α * δNpi * div_Nuj) * dΩ
            end
        end
    end


end

