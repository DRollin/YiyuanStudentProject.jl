using YiyuanStudentProject
import GLMakie
import LinearSolve, LinearSolvePardiso

dim = 3

P = Material{dim}(; G=4.0e5, K=6.0e5, η=1.0, cʳᵉᶠ=1.0, μʳᵉᶠ=1.0, θʳᵉᶠ=1.0, cᵐ=1.0, α=0.0)
M = Material{dim}(; G=2.0e5, K=3.0e5, η=1.0, cʳᵉᶠ=1.0, μʳᵉᶠ=1.0, θʳᵉᶠ=1.0, cᵐ=1.0, α=0.0)
load = LoadCase(dim; ε̄₁₁=0.01, μ̄=1.0, ζ̄₁=0.1)

d = 0.2
ϕ = 0.3 #packing fraction
meshsize = 0.1
grid = generate_rve_grid(; ϕ=ϕ, d=d, meshsize=meshsize, dx=(1.0,1.0, 1.0))
rve = RVE(grid, P, M)


res, setup = solve_time_series(rve, load; Δt=1e-4, t_total=1e-3)

#plot_grid(grid)

file, fig, anim = animate_result(res, setup, file_name="Myresult.mp4", n=10.0)