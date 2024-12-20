# # [Example 1](@id example-1)
#
#
# In this example, we solve the transit linear chemo-mechanical problem in FE² framework that is discussed the documentation in a domain with particles ($\Omega^\text{P}$)
# embedded in matrix ($\Omega^\text{M}$). 
#
#md # The full program, without comments, can be found in the next [section](@ref example_1-plain-program).
#
using YiyuanStudentProject
#
# ## Fine scale
# Setting up the dimension, material parameters, `loadcase` originated from macro scale for solving a 3D RVE problem. The material parameters are characterized using simple values
# for the particle `P` and Matrix `M`. Parameters like the average particle diameter `d`, the particle volume fraction `ϕ`, and the `meshsize` are given followed.
dim = 3

P = Material{dim}(; G=4.0e5, K=6.67e5, η=1.0, cʳᵉᶠ=1.0, μʳᵉᶠ=1.0, θʳᵉᶠ=1.0, cᵐ=1.0, α=0.0001)
M = Material{dim}(; G=4.0e5, K=6.67e5, η=1.0, cʳᵉᶠ=1.0, μʳᵉᶠ=1.0, θʳᵉᶠ=1.0, cᵐ=1.0, α=0.0001)
load = LoadCase(dim; ε̄₁₁ = 0.25, μ̄ = 1.0, ζ̄₁=1.0)
d = 0.1 
ϕ = 0.1
meshsize = 0.1
#
# Generate a grid with spherical inclusions for a 3D Representative Volume Element (RVE) in desired meshsize. 
# Material phases like partical and matrix are included as cell sets.
grid = generate_rve_grid(; ϕ=ϕ, d=d, meshsize=meshsize, dx=(1.0,1.0, 1.0))
# Construct an object as setup for solving the RVE problem.
rve = RVE(grid, P, M)
setup = prepare_setup(rve)
# Perform the assembly for constructing the system matrices K and M.
(; aⁿ, aⁿ⁺¹) = setup
assemble_K_M!(setup)
# Solve the time dependent problem applying the Crank-Nicolson Method.
Δt=0.1
setup, aⁿ⁺¹ = compute_time_step!(setup, load, Δt)

#md # ## [Plain program](@id example_1-plain-program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`example_1.jl`](example_1.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```