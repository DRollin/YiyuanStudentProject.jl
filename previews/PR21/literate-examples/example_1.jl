# # [Example 1](@id example-1)
#
# ![](description_problem.png)
# *Figure 1* Schematic illustration of a FE² framwork with two material phases (Particles in blue and matrix in gray) in RVE and a possible boundary condition on macro scale object 
# as well as required parameters for a transient chemo-mechanical problem
#
#
# In this example, we will solve the transient linear chemo-mechanical problem in a FE² framework. The mathmatical model can be found in the documentation. A domain with particles ($\Omega^\text{P}$)
# embedded in matrix ($\Omega^\text{M}$) is prescribed. Parameters for solving the problem are defined for both sub scale (RVE) and macro scale. A macro scale boundary condition is as well selected
# to constrain the object.
#
# A main function will be defined which sets up everything and solves the problem with a given time step size and total simulation time.
#
#md # The full program, without comments, can be found in the next [section](@ref example_1-plain-program).
#
using YiyuanStudentProject
#
# ## Preparation
# ### Sub Scale (RVE)
# Setting up the dimension, material parameters. The material parameters are cataloged with `P` for particles and `M` for matrix and 
# are characterized using simple values for now.
# 
#
# Generate a grid with spherical inclusions for a 3D Representative Volume Element (RVE) in desired meshsize and grid size `dx = (1.0, 1.0, 1.0)`
function sub_scale_pre()
    dim = 3
    P = Material{dim}(; G=4.5e5, K=8.0e5, η=0.5, cʳᵉᶠ=0.0, μʳᵉᶠ=0.0, θʳᵉᶠ=1.0, cᵐ=1.0, α=0.2)
    M = Material{dim}(; G=0.8e5, K=1.7e5, η=1.0, cʳᵉᶠ=0.0, μʳᵉᶠ=0.0, θʳᵉᶠ=1.0, cᵐ=1.0, α=0.0)
    d = 0.2
    ϕ = 0.3
    meshsize = 0.1
    grid = generate_rve_grid(; ϕ=ϕ, d=d, meshsize=meshsize, dx=(1.0, 1.0, 1.0))
    return grid, P, M
end
# ### Macro Scale
# Setting up the grid parameters and generate a macro scale grid. 
#
function macro_scale_pre()
    element_number = (1, 1, 1)
    x̄_l = Vec{3, Float64}((-50, -50, -50))
    x̄_r = Vec{3, Float64}((50, 50, 50))
    grid_macro = generate_grid(Tetrahedron, element_number , x̄_l, x̄_r)
    return grid_macro
end
# ## Solve
# Construct an object `RVESetup` named `setup_rve` with relevent infomation for solving the RVE problem.
#
# Due to the linearity of the problem, gloable stiffness matrix `K` and gloable mass matrix `M` as well as gloable right hand side vector `f` remain the same
# through out the time stepping. Perform the assembly for constructing the system matrices K, M, and right hand side vector f.
#
# Create an appropriate boundary condition using object `MacroBCParams`
#
# Solve the time dependent problem macro scale problem. A corresponding RVE would be solved at each quadrature point in macro scale problem. 
#
# A visualization of results from both fine scale and macro scale problem can use `animate_combined_result`
function solve(grid_rve, P, M, grid_macro, Δt, t_total)
    rve = RVE(grid_rve, P, M)
    setup_rve = prepare_setup(rve)
    assemble_K_M_f!(setup_rve)
    bc = MacroBCParams(Vec{3}((0.01, 0.0, 0.0)), Vec{3}((0.0, 0.0, 0.0)), 1.0, 0.0, "left", "right")
    res, res_rve, setup = solve_macro_problem(grid_macro, setup_rve, rve, bc, Δt=Δt, t_total=t_total)
end
# Construct the main function
function example()
    grid_rve, P, M = sub_scale_pre()
    grid_macro = macro_scale_pre()
    Δt = 1e-6
    t_total = 1e-5
    solve(grid_rve, P, M, grid_macro, Δt, t_total)
end
# Run the simulation
example()
# 
#
# 

#md # ## [Plain program](@id example_1-plain-program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`example_1.jl`](example_1.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```