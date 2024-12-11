using YiyuanStudentProject
import GLMakie

dim = 3

P = Material{dim}(; G=4.0e5, K=6.67e5, η=1.0, cʳᵉᶠ=1.0, μʳᵉᶠ=1.0, θʳᵉᶠ=1.0, cᵐ=1.0, α=0.0001)
M = Material{dim}(; G=4.0e5, K=6.67e5, η=1.0, cʳᵉᶠ=1.0, μʳᵉᶠ=1.0, θʳᵉᶠ=1.0, cᵐ=1.0, α=0.0001)
load = LoadCase(dim; ε̄₁₁ = 0.25, μ̄ = 1.0, ζ̄₁=1.0)

d = 0.1 
ϕ = 0.1 #packing fraction
meshsize = 0.1
grid = generate_rve_grid(; ϕ=ϕ, d=d, meshsize=meshsize, dx=(1.0,1.0, 1.0))
rve = RVE(grid, P, M)

res, setup = solve_time_series(rve, load)

#plot_grid(grid)

file, fig, anim = animate_result(res, setup, file_name="Myresult.mp4", n=0.5)