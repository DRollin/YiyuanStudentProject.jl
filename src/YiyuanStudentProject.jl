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

import LinearSolve, Pardiso

using ProgressLogging
using Logging: global_logger
using TerminalLoggers: TerminalLogger

global_logger(TerminalLogger())



δ(i,j) = i == j ? 1 : 0

include("types.jl")
export Material, RVE, LoadCase, PhaseSetup, FESetup, AssemblySetup, SolveSetup, GaussPointData, MacroBCParams

include("sub_scale/mesh_generation.jl")
export generate_rve_spheres, generate_rve_grid, generate_box_grid, 
        generate_spheres, get_radius_pdf

include("sub_scale/setup.jl")
export prepare_setup

include("sub_scale/assembly.jl")
export assemble_K_M_f!, assemble_element!

include("sub_scale/solve.jl")
export compute_time_step!, solve_time_series
    
include("sub_scale/plotting.jl")
export animate_result

include("upscaling/averaging.jl")
export compute_effective_response, average_quantities

include("macro_scale/setup.jl")
export add_macro_bc!, prepare_macro_setup

include("macro_scale/assembly.jl")
export assemble_macro_K!, assemble_macro_element!

include("macro_scale/solve.jl")
export solve_macro_problem

include("macro_scale/plotting.jl")
export animate_macro_result, animate_combined_result

end
