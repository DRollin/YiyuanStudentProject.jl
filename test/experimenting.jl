using YiyuanStudentProject
import GLMakie

dim = 3

P = Material{dim}(; G=4.5e5, K=8.0e5, η=0.5, cʳᵉᶠ=0.0, μʳᵉᶠ=0.0, θʳᵉᶠ=298.5, cᵐ=1.0, α=2.0)
M = Material{dim}(; G=0.8e5, K=1.7e5, η=1.0, cʳᵉᶠ=0.0, μʳᵉᶠ=0.0, θʳᵉᶠ=298.5, cᵐ=1.0, α=0.0) #G=0.8e5, K=1.7e5
#load = LoadCase(dim; ε̄₁₁=0.1, μ̄=1.0, ζ̄₁=0.1)

d = 0.2
ϕ = 0.3 
meshsize = 0.1
grid_rve = generate_rve_grid(; ϕ=ϕ, d=d, meshsize=meshsize, dx=(1.0,1.0,1.0))
rve = RVE(grid_rve, P, M)

nx = 2
grid_macro = generate_grid(Tetrahedron, (nx,nx,nx) , Vec{3, Float64}((-nx*50,-nx*50,-nx*50)), Vec{3, Float64}((nx*50,nx*50,nx*50)))

setup_rve = prepare_setup(rve)

assemble_K_M_f!(setup_rve)

#res, setup= solve_time_series(rve, load;  Δt=1e-4, t_total=1e-3)
#animate_result(res, setup; file_name ="Myresult.mp4")

bc = MacroBCParams(Vec{3}((0.01, 0.0, 0.0)), Vec{3}((0.0, 0.0, 0.0)), 1.0, 0.0, "left", "right")

res, res_rve, setup = solve_macro_problem(grid_macro, setup_rve, rve, bc, Δt=1e-6, t_total=1e-5)

file, fig, anim = animate_combined_result(res, res_rve, setup, file_name="Myresult.mp4", scale=1000.0)