module YiyuanStudentProject

using Reexport
@reexport using Ferrite
using FerriteGmsh
using SparseArrays, OrderedCollections


using BubbleBath
using Distributions: Uniform
using LaTeXStrings
using WriteVTK
using ForwardDiff

import Makie: Makie, GeometryBasics

#include("C:\\Users\\Yiyua\\.julia\\packages\\Gmsh\\vuFrI\\src\\Gmsh.jl")

# TODO:
# - Improve Plotting
# - Prepare upscaling -> Avergage Flux j̄, and stress σ̄
# - Check boundary conditions for ū, μ̄
#       -> 1. find node id of the node closest to coordinate origin (center of RVE)
#       -> Dirichlet(:u, OrderedSet{Int}([nodeid]), (x,t) -> zero(Vec{dim}))
#          Dirichlet(:μ, OrderedSet{Int}([nodeid]), (x,t) -> μ̄)
# - Change initial state: u=0, c=cref, μ=μref -> apply_analytical!()
# - Add tests and docs

δ(i,j) = i == j ? 1 : 0

include("types.jl")
export Material, RVE, LoadCase, PhaseSetup, FESetup

include("fine_scale/mesh_generation.jl")
export generate_rve_spheres, generate_rve_grid, generate_box_grid, 
    add_sphere_to_gmsh, generate_spheres, get_radius_pdf, get_tags_from_dimtags

include("fine_scale/mesh_characteristics.jl")
export get_volume

include("fine_scale/setup.jl")
export prepare_setup

include("fine_scale/assembly.jl")
export assemble!, assemble_element!

include("fine_scale/solve.jl")
export compute_time_step!, solve_time_series
    
include("fine_scale/plotting.jl")
export plot_grid, plot_rve_grid, plot_mesh_overlay!, plot_potential, select_state, animate_result

#include("fine_scale/averaging.jl")
export compute_homogenized_potential, compute_homogenized_gradient
    # Upscaling
#=
include("upscaling/averaging.jl")
include("upscaling/sensitivities.jl")

export EffectiveResponse,
    Sensitivities, compute_sensitivities, average_bulk_quantities, compute_effective_response

    # Macro scale
    
include("macro_scale/problem_setup.jl")
include("macro_scale/assembly.jl")
include("macro_scale/plotting.jl")
include("macro_scale/solve.jl")


export solve_macro_problem, plot_macro_potential
=#

end
