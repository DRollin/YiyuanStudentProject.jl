
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

Ruturn a Makie.Figure to visualize a 3D finite element grid using the Makie plotting library.

# Arguments:
- `grid::Grid{3}`: A 3D grid object representing the finite element mesh to be visualized. 


# Implementation Details:
A plotable mesh is generated using `_prepare_plotable_mesh` and input `grid`

Sets up a 3D axis `Axis3` with an equal aspect ratio and a title `undeformed grid`.

Renders the grid as a solid mesh `Makie.mesh!` with a light blue color.

Overlays the grid's edges as a wireframe `Makie.wireframe!` with black edges.

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
    animate_result(res::NamedTuple, setup::RVESetup{dim}, file_name::String="Myresult.mp4", n::Number)

Return an animation showing the evolution of the solution for a 3D RVE simulation.

# Arguments:
- `res`:         A `NamedTuple` containing the simulation results
- `setup`:       The setup object for the RVE simulation, containing: `grid` and `dh` fields
- `file_name`:   The path and name of the output animation file (default: `"Myresult.mp4"`)

# Implementation Details:
Initialize a Makie figure with a size of `(1200, 800)`.

Generate observables for `t` and `a` using `_prepare_plots!`.

Record frames by iterating over the indices of `t` and updates the observables accordingly.

Save the animation as an MP4 file.
    
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
   _prepare_plots!(pos, res::NamedTuple, setup::RVESetup{dim}; 
                    scale::Real=1.0, 
                    titlestart::String="RVE response") where {dim}

    Sets up the necessary plots and observables `tᵒᵇˢ` and `aᵒᵇˢ` for visualizing macroscale simulation results: deformation, chemical potential, concentration.

# Arguments:
- `pos`:    A layout position,
- `res`:    A `NamedTuple` containing:
    - `t`:      Time at each time step,
    - `a`:      Result vector corresponding to the time.
- `setup`: An object `SolveSetup`with the grid and degrees of freedom configuration.

# Implementation Details:
Generate the mesh representation both undeformed and deformed grid.

Determine the boundary values of chemical potential `μ` and concentration `c` across the grid nodes for all time steps.

Add sliced opening to show the internal strucature with particles and matrix with distinguished colors.

Create observables for time `tᵒᵇˢ` and data `aᵒᵇˢ`.
    
Prepare subplots for:
    Undeformed grid with a static wireframe.
    Deformed grid based on the scaled displacement field `u`.
    Chemical potential and concentration using a color-mapped mesh and color bar.

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