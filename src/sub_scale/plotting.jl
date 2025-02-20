
function _convert_cells(::Grid, cells::Vector{Tetrahedron})
    return [ GeometryBasics.TriangleFace{Int}(facet...) for cell in cells for facet in Ferrite.facets(cell) ]
end
_convert_cells(grid::Grid, set::OrderedSet{Int}) = _convert_cells(grid, getcells(grid, collect(set)))
_convert_cells(grid::Grid, set::String) = _convert_cells(grid, getcells(grid, set))

function _prepare_plotable_mesh(grid::Grid{dim}, set) where {dim}
    nodes = [ GeometryBasics.Point{dim,Float64}(n.x...) for n in grid.nodes ]
    cells = _convert_cells(grid, set)
    return GeometryBasics.Mesh(nodes, cells)
end
function _prepare_plotable_mesh(grid::Grid{dim}, displ::Vector{<:Vec{dim}}, set) where {dim}
    nodes = [ GeometryBasics.Point{dim,Float64}( (n.x .+ u)...) for (n, u) in zip(grid.nodes, displ) ]
    cells = _convert_cells(grid, set)
    return GeometryBasics.Mesh(nodes, cells)
end
_prepare_plotable_mesh(grid::Grid) = _prepare_plotable_mesh(grid, grid.cells)
_prepare_plotable_mesh(grid::Grid{dim}, u::Vector{<:Vec{dim}}) where {dim} = _prepare_plotable_mesh(grid, u, grid.cells)

"""
    plot_grid(grid::Grid{3})

Ruturn a Makie.Figure to visualize a given 3D finite element `grid` using the Makie plotting library.

"""
function plot_grid(grid::Grid{3})
    mesh = _prepare_plotable_mesh(grid)
    fig = Makie.Figure()
    ax  = Makie.Axis3(fig[1, 1], aspect = :equal, title = "undeformed grid")
    Makie.mesh!(ax, mesh; color=Makie.RGB(0.5,0.5,1.0))
    Makie.wireframe!(ax, mesh; color=:black)
    return fig
end


"""
    animate_result(res::NamedTuple, setup::RVESetup{dim}, file_name::String="Myresult.mp4", kwargs...)

Return an mp4 file with name `file_name`, which by default is "Myresult.mp4". 

This shows the evolution of the solution `res` for a 3D RVE simulation based on input `RVESetup`. 

Further arguments like `scale` for visibel deformation can be used.
   
"""
function animate_result(res::NamedTuple, setup::RVESetup{dim}; file_name ="Myresult.mp4", kwargs...) where {dim}
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

"""
    _prepare_plots!(pos, res::NamedTuple, setup::RVESetup{dim}; scale::Real=1.0, titlestart::String="RVE response") where {dim}

Sets up the necessary plots with a specified layout position `pos` and observables `tᵒᵇˢ` and `aᵒᵇˢ` to visualize the simulation results `res` of a RVE problem `RVESetup`:
deformation `u`, chemical potential `μ`, concentration `c`. 

Optional arguments include a scaling factor `scale` for visualizing the deformation and a title `titlestart`. By default the `scale` is set to `1.0` and `titlestart` to `RVE response`.

"""
function _prepare_plots!(pos, res::NamedTuple, setup::RVESetup{dim};
        scale::Real=1.0, 
        titlestart::String="RVE response") where {dim}

    (; grid, dh) = setup
    (; a, t) = res

    μ_all = [evaluate_at_grid_nodes(dh, a[i], :μ) for i in eachindex(t)]
    c_all = [evaluate_at_grid_nodes(dh, a[i], :c) for i in eachindex(t)]
    μ_min, μ_max = minimum(μ -> minimum(μ), μ_all), maximum(μ -> maximum(μ), μ_all)
    c_min, c_max = minimum(c -> minimum(c), c_all), maximum(c -> maximum(c), c_all)
     
    addcellset!(grid, "sliced open X", x -> x[1] ≥ 0)
    addcellset!(grid, "sliced open Y", x -> x[2] ≥ 0)
    addcellset!(grid, "sliced open Z", x -> x[3] ≤ 0)
    cells = union(getcellset(grid, "sliced open X"), getcellset(grid, "sliced open Y"), getcellset(grid, "sliced open Z"))
    delete!(grid.cellsets, "sliced open X")
    delete!(grid.cellsets, "sliced open Y")
    delete!(grid.cellsets, "sliced open Z")
    cellsᴾ = setdiff(cells, getcellset(grid, "matrix"))
    cellsᴹ = setdiff(cells, getcellset(grid, "particles"))

    mesh  = _prepare_plotable_mesh(grid, cells)
    meshᴾ = _prepare_plotable_mesh(grid, cellsᴾ)
    meshᴹ = _prepare_plotable_mesh(grid, cellsᴹ)
    
    tᵒᵇˢ = Makie.Observable(t[1])
    aᵒᵇˢ = Makie.Observable(a[1])

    title = Makie.@lift titlestart*" at t=$( round($(tᵒᵇˢ); sigdigits=4) )"
    Makie.Label(pos[1,1:2], title)

    ax  = Makie.Axis3(pos[2,1], aspect=:equal, title="RVE undeformed grid")
    Makie.mesh!(ax, meshᴾ; color=Makie.RGB(1.0,0.5,0.5), shading=Makie.NoShading)
    Makie.mesh!(ax, meshᴹ; color=Makie.RGB(0.5,0.5,1.0), shading=Makie.NoShading)
    Makie.wireframe!(ax, mesh; color=:black)
    
    u = Makie.@lift evaluate_at_grid_nodes(dh, $(aᵒᵇˢ), :u)
    defmeshᴾ = Makie.@lift _prepare_plotable_mesh(grid, ( $(u) .* scale ), cellsᴾ)
    defmeshᴹ = Makie.@lift _prepare_plotable_mesh(grid, ( $(u) .* scale ), cellsᴹ)
    ax = Makie.Axis3(pos[2,2], aspect=:equal, title="RVE deformed grid")
    Makie.mesh!(ax, defmeshᴾ; color=Makie.RGB(1.0,0.5,0.5), shading=Makie.NoShading)
    Makie.mesh!(ax, defmeshᴹ; color=Makie.RGB(0.5,0.5,1.0), shading=Makie.NoShading)
    Makie.wireframe!(ax, defmeshᴾ; color=:black)
    Makie.wireframe!(ax, defmeshᴹ; color=:black)

    subpos = pos[3,1]
    μ = Makie.@lift evaluate_at_grid_nodes(dh, $(aᵒᵇˢ), :μ)
    ax = Makie.Axis3(subpos[1,1], aspect=:equal, title="RVE chemical potential")
    colorsettings = (colorrange=(0.0, 1.0), colormap=:viridis)
    Makie.mesh!(ax, mesh; color=μ, colorsettings..., shading=Makie.NoShading)
    Makie.Colorbar(subpos[1,2]; colorsettings...)
   
    subpos = pos[3,2]
    c = Makie.@lift evaluate_at_grid_nodes(dh, $(aᵒᵇˢ), :c)
    ax = Makie.Axis3(subpos[1,1], aspect=:equal, title="RVE concentration")
    colorsettings = (colorrange=(c_min, c_max), colormap=:viridis)
    Makie.mesh!(ax, mesh; color=c, colorsettings..., shading=Makie.NoShading)
    Makie.Colorbar(subpos[1,2]; colorsettings...)

    return tᵒᵇˢ, aᵒᵇˢ
end