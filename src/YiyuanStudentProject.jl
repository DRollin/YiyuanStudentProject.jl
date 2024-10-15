module YiyuanStudentProject

using Reexport
@reexport using Ferrite
using FerriteGmsh
using SparseArrays


using BubbleBath
using Distributions: Uniform
using LaTeXStrings
using WriteVTK
using ForwardDiff

import Makie, GeometryBasics

#include("C:\\Users\\Yiyua\\.julia\\packages\\Gmsh\\vuFrI\\src\\Gmsh.jl")

include("types.jl")

include("fine_scale/mesh_generation.jl")
include("fine_scale/mesh_characteristics.jl")

include("fine_scale/problem_setup.jl")
include("fine_scale/assembly.jl")
include("fine_scale/solve.jl")
    
#include("fine_scale/plotting.jl")
#include("fine_scale/averaging.jl")


export pre_Material, iso_pv_Material, 
       RVEProblem, FESetup, LoadCase, iso_pv_ElementSetup,
    generate_rve_spheres, generate_rve_grid, get_volume,
    plot_grid, plot_rve_grid, plot_mesh_overlay!, plot_potential, select_state,
    prepare_setup, solve_load_case, solve_RVE, 
    compute_homogenized_potential, compute_homogenized_gradient, generate_box_grid, add_sphere_to_gmsh, generate_spheres,
    get_radius_pdf, get_tags_from_dimtags


    # Upscaling
#=include("upscaling/averaging.jl")
include("upscaling/sensitivities.jl")

export EffectiveResponse,
    Sensitivities, compute_sensitivities, average_bulk_quantities, compute_effective_response

    # Macro scale
    
include("macro_scale/problem_setup.jl")
include("macro_scale/assembly.jl")
include("macro_scale/plotting.jl")
include("macro_scale/solve.jl")


export solve_macro_problem, plot_macro_potential=#

end
