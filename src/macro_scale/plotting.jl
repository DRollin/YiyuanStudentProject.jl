
"""
    animate_macro_result(res::NamedTuple, setup::RVESetup{dim}, file_name::String="Myresult.mp4", n::Number)

Return an animation showing the evolution of the solution for a 3D RVE simulation.

# Arguments:
- `res`:         A `NamedTuple` containing the simulation results
- `setup`:       The setup object for the RVE simulation, containing: `grid` and `dh` fields
- `file_name`:   The path and name of the output animation file (default: `"Myresult.mp4"`)
- `n`:           Scaling factor for displacement

# Implementation Details:
A plotable mesh is generated using `_prepare_plotable_mesh` with a cut open to show the inner structure.

A Figure is plotted showing the results of for displacement `u`, chemical potantial `μ`, concentration`c` at each time step.
    
"""
function animate_macro_result(res::NamedTuple, setup::SolveSetup{dim}; file_name ="Myresult.mp4", kwargs...) where {dim}
    (; t, a) = res
    fig = Makie.Figure(size=(1200,800))
    tᵒᵇˢ, aᵒᵇˢ = _prepare_plots!(fig, res, setup; kwargs...)

	file = joinpath(file_name)
	anim = Makie.record(fig, file, eachindex(t); framerate=1) do i
        tᵒᵇˢ[] = t[i]
        aᵒᵇˢ[] = a[i]
	end
	return file, fig, anim
end

function _prepare_plots!(pos, res::NamedTuple, setup::SolveSetup{dim};
        scale::Real=1.0, 
        titlestart::String="macroscale solution") where {dim}
    (; grid, dh) = setup
    (; t, a)= res

    μ_all = [evaluate_at_grid_nodes(dh, res.a[i], :μ) for i in eachindex(res.t)]
    μ_min, μ_max = minimum(μ -> minimum(μ), μ_all), maximum(μ -> maximum(μ), μ_all)
         
    #=addcellset!(grid, "sliced open X", x -> x[1] ≥ 0)
    addcellset!(grid, "sliced open Y", x -> x[2] ≥ 0)
    addcellset!(grid, "sliced open Z", x -> x[3] ≤ 0)
    cells = union(getcellset(grid, "sliced open X"), getcellset(grid, "sliced open Y"), getcellset(grid, "sliced open Z"))
    delete!(grid.cellsets, "sliced open X")
    delete!(grid.cellsets, "sliced open Y")
    delete!(grid.cellsets, "sliced open Z")=#
    mesh = _prepare_plotable_mesh(grid)

    tᵒᵇˢ = Makie.Observable(t[1])
    aᵒᵇˢ = Makie.Observable(a[1])

    title = Makie.@lift titlestart*" at t=$( round($(tᵒᵇˢ); sigdigits=4) )"
    Makie.Label(pos[1,1:2], title)

    ax  = Makie.Axis3(pos[2,1], aspect=:equal, title="undeformed grid")
    Makie.mesh!(ax, mesh; color=Makie.RGB(1.0,0.5,0.5), shading=Makie.NoShading)
    Makie.wireframe!(ax, mesh; color=:black)
    
    u = Makie.@lift evaluate_at_grid_nodes(dh, $(aᵒᵇˢ), :u)
    defmesh = Makie.@lift _prepare_plotable_mesh(grid, ( $(u) .* scale ))
    ax = Makie.Axis3(pos[2,2], aspect=:equal, title="deformed grid")
    Makie.mesh!(ax, defmesh; color=Makie.RGB(1.0,0.5,0.5), shading=Makie.NoShading)
    Makie.wireframe!(ax, defmesh; color=:black)

    subpos = pos[3,1]
    μ = Makie.@lift evaluate_at_grid_nodes(dh, $(aᵒᵇˢ), :μ)
    ax = Makie.Axis3(subpos[1,1], aspect=:equal, title="chemical potential")
    colorsettings = (colorrange=(μ_min, μ_max), colormap=:viridis)
    Makie.mesh!(ax, mesh; color=μ, colorsettings..., shading=Makie.NoShading)
    Makie.Colorbar(subpos[1,2]; colorsettings...)

    return tᵒᵇˢ, aᵒᵇˢ
end


function animate_combined_result(res::NamedTuple, setup::SolveSetup{dim}; file_name ="Myresult.mp4", scale=1.0) where {dim}
    (; t, a) = res
    (; rvesetup, gpdata) = setup.assemblysetup

    fig = Makie.Figure(size=(1200,800))
    tᵒᵇˢ, aᵒᵇˢ = _prepare_plots!(fig[1,1], res, setup; scale=scale)

        # Create some dummy RVE data -> TODO: save RVE-responses during the computation and use here as argument of the function
    resᵣᵥₑ = [(cellid=1, qp=i, 
               res=(t=[t[n] for n in eachindex(t)], 
                    a=[deepcopy(rvesetup.aⁿ) for n in eachindex(t)]) 
                     ) for i in 1:2]
    rveᵒᵇˢ = [_prepare_plots!(fig[i+1,1], resᵣᵥₑ[i].res, rvesetup; titlestart="RVE $(i)") for i in eachindex(resᵣᵥₑ)]

	file = joinpath(file_name)
	anim = Makie.record(fig, file, eachindex(t); framerate=1) do i
        tᵒᵇˢ[] = t[i]
        aᵒᵇˢ[] = a[i]
        for j in eachindex(resᵣᵥₑ)
            rveᵒᵇˢ[j][1][] = t[i]
            rveᵒᵇˢ[j][2][] = resᵣᵥₑ[j].res.a[i]
        end
	end
	return file, fig, anim
end