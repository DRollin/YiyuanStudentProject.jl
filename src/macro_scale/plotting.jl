
"""
    animate_macro_result(res::NamedTuple, setup::RVESetup{dim}, file_name::String="Myresult.mp4", n::Number)

Return an animation showing the evolution of the solution for a 3D RVE simulation.

# Arguments:
- `res`:         A `NamedTuple` containing the simulation results
- `setup`:       The setup object for the RVE simulation, containing: `grid` and `dh` fields
- `file_name`:   The path and name of the output animation file (default: `"Myresult.mp4"`)
- `n`:           Scaling factor for displacement

# Implementation Details:
A plotable mesh is generated using `prepare_plotable_mesh` with a cut open to show the inner structure.

A Figure is plotted showing the results of for displacement `u`, chemical potantial `μ`, concentration`c` at each time step.
    
"""
function animate_macro_result(res::NamedTuple, setup::SolveSetup{dim}; file_name ="Myresult.mp4", n=1.0) where {dim}
    (; grid, dh) = setup
    
    μ_all = [evaluate_at_grid_nodes(dh, res.a[i], :μ) for i in eachindex(res.t)]

    μ_min, μ_max = minimum(μ -> minimum(μ), μ_all), maximum(μ -> maximum(μ), μ_all)
     
    #=addcellset!(grid, "sliced open X", x -> x[1] ≥ 0)
    addcellset!(grid, "sliced open Y", x -> x[2] ≥ 0)
    addcellset!(grid, "sliced open Z", x -> x[3] ≤ 0)
    cells = union(getcellset(grid, "sliced open X"), getcellset(grid, "sliced open Y"), getcellset(grid, "sliced open Z"))
    delete!(grid.cellsets, "sliced open X")
    delete!(grid.cellsets, "sliced open Y")
    delete!(grid.cellsets, "sliced open Z")=#

    mesh = prepare_plotable_mesh(grid)

    fig = Makie.Figure(size=(1200,800))


    t = res.t[1]
    tᵒᵇˢ = Makie.Observable(t)
    title = Makie.@lift "solution at t=$( round($(tᵒᵇˢ); sigdigits=4) )"
    Makie.Label(fig[1,1:2], title)

    
    ax  = Makie.Axis3(fig[2,1], aspect=:equal, title="undeformed grid")
    Makie.mesh!(ax, mesh; color=Makie.RGB(1.0,0.5,0.5), shading=Makie.NoShading)
    Makie.wireframe!(ax, mesh; color=:black)
    
    u = evaluate_at_grid_nodes(dh, res.a[1], :u)
    uᵒᵇˢ = Makie.Observable(u)
    defmesh = Makie.@lift prepare_plotable_mesh(grid, ( $(uᵒᵇˢ) .* n ))
    ax = Makie.Axis3(fig[2,2], aspect=:equal, title="deformed grid")
    Makie.mesh!(ax, defmesh; color=Makie.RGB(1.0,0.5,0.5), shading=Makie.NoShading)
    Makie.wireframe!(ax, defmesh; color=:black)

    pos = fig[3,1]
    μ = evaluate_at_grid_nodes(dh, res.a[1], :μ)
    μᵒᵇˢ = Makie.Observable(μ)
    ax = Makie.Axis3(pos[1,1], aspect=:equal, title="chemical potential")
    colorrange = (μ_min, μ_max)
    Makie.mesh!(ax, mesh; color=μᵒᵇˢ, colormap=:viridis, colorrange=colorrange, shading=Makie.NoShading)
    Makie.Colorbar(pos[1,2]; colormap=:viridis, colorrange=colorrange)
      

	file = joinpath(file_name)
	anim = Makie.record(fig, file, eachindex(res.t); framerate=1) do i
        tᵒᵇˢ[] = res.t[i]
        uᵒᵇˢ[] = evaluate_at_grid_nodes(dh, res.a[i], :u)
        μᵒᵇˢ[] = evaluate_at_grid_nodes(dh, res.a[i], :μ)
	end

	return file, fig, anim
end