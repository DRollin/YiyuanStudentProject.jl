using YiyuanStudentProject

dim = 3

P = Material{dim}(; G=4.0e5, K=6.67e5, η=1.0, cʳᵉᶠ=1.0, μʳᵉᶠ=1.0, θʳᵉᶠ=1.0, cᵐ=1.0, α=1.0)
M = Material{dim}(; G=4.0e5, K=6.67e5, η=1.0, cʳᵉᶠ=1.0, μʳᵉᶠ=1.0, θʳᵉᶠ=1.0, cᵐ=1.0, α=1.0)
load = LoadCase(dim; ζ̄₁=1.0)
#TODO: Check, how to apply μ̄ properly (I think we need some additional constraint for that or a Dirichlet BC on one node)

d = 0.1 
ϕ = 0.1 #packing fraction
meshsize = 0.1
grid = generate_rve_grid(; ϕ=ϕ, d=d, meshsize=meshsize, dx=(1.0,1.0, 1.0))
rve = RVE(grid, P, M)

res, setup = solve_load_case(rve, load)

plot_grid(grid)

file, fig, anim = animate_result(res, setup)