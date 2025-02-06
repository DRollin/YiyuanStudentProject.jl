using YiyuanStudentProject

function example(Δt, t_total)
    #Preparation

    #sub scale
    dim = 3
    P = Material{dim}(; G=4.5e5, K=8.0e5, η=0.5, cʳᵉᶠ=0.0, μʳᵉᶠ=0.0, θʳᵉᶠ=1.0, cᵐ=1.0, α=0.2)
    M = Material{dim}(; G=0.8e5, K=1.7e5, η=1.0, cʳᵉᶠ=0.0, μʳᵉᶠ=0.0, θʳᵉᶠ=1.0, cᵐ=1.0, α=0.0)
    d = 0.2
    ϕ = 0.3
    meshsize = 0.1
    grid = generate_rve_grid(; ϕ=ϕ, d=d, meshsize=meshsize, dx=(1.0, 1.0, 1.0))

    #macro scale
    element_number = (1, 1, 1)
    x̄_l = Vec{3, Float64}((-50, -50, -50))
    x̄_r = Vec{3, Float64}((50, 50, 50))
    grid_macro = generate_grid(Tetrahedron, element_number , x̄_l, x̄_r)

    #solve
    rve = RVE(grid, P, M)
    setup_rve = prepare_setup(rve)
    assemble_K_M_f!(setup_rve)
    bc = MacroBCParams(Vec{3}((0.01, 0.0, 0.0)), Vec{3}((0.0, 0.0, 0.0)), 1.0, 0.0, "left", "right")
    res, res_rve, setup = solve_macro_problem(grid_macro, setup_rve, rve, bc, Δt=Δt, t_total=t_total)

    #file, fig, anim = animate_combined_result(res, res_rve, setup, file_name="Myresult.mp4", scale=1000.0)
end

example(1e-6, 1e-5)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
