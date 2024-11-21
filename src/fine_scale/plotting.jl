
function convert_cells(::Grid, cells::Vector{Tetrahedron})
    return [ GeometryBasics.TriangleFace{Int}(facet...) for cell in cells for facet in Ferrite.facets(cell) ]
end
convert_cells(grid::Grid, set::OrderedSet{Int}) = convert_cells(grid, getcells(grid, collect(set)))
convert_cells(grid::Grid, set::String) = convert_cells(grid, getcells(grid, set))

function prepare_plotable_mesh(grid::Grid{dim}, set) where {dim}
    nodes = [ GeometryBasics.Point{dim,Float64}(n.x...) for n in grid.nodes ]
    cells = convert_cells(grid, set)
    return GeometryBasics.Mesh(nodes, cells)
end
function prepare_plotable_mesh(grid::Grid{dim}, displ::Vector{<:Vec{dim}}, set) where {dim}
    nodes = [ GeometryBasics.Point{dim,Float64}( (n.x .+ u)...) for (n, u) in zip(grid.nodes, displ) ]
    cells = convert_cells(grid, set)
    return GeometryBasics.Mesh(nodes, cells)
end
prepare_plotable_mesh(grid::Grid) = prepare_plotable_mesh(grid, grid.cells)
prepare_plotable_mesh(grid::Grid{dim}, u::Vector{<:Vec{dim}}) where {dim} = prepare_plotable_mesh(grid, u, grid.cells)

"""
    TODO
"""
function plot_grid(grid::Grid{3})
    mesh = prepare_plotable_mesh(grid)
    fig = Makie.Figure()
    ax  = Makie.Axis3(fig[1, 1], aspect = :equal, title = "undeformed grid")
    Makie.mesh!(ax, mesh; color=Makie.RGB(0.5,0.5,1.0))
    Makie.wireframe!(ax, mesh; color=:black)
    return fig
end

# TODO: Idea -> have file name as input argument, so the user can decide where to save the animation
"""
    TODO
"""
function animate_result(res::NamedTuple, setup::RVESetup{dim}) where {dim}
    (; grid, dh) = setup
    # TODO: List of ideas
    # - Maybe introduce a scaling factor for the displacement?
    # - Find a better (larger) size for the Figure
    # - Find initial color ranges which are adapted to the applied loading or the overall max/min?
    # - Improve ticks if necessary to match RVE bounds
    addcellset!(grid, "sliced open X", x -> x[1] ≥ 0)
    addcellset!(grid, "sliced open Y", x -> x[2] ≥ 0)
    addcellset!(grid, "sliced open Z", x -> x[3] ≤ 0)
    cells = union(getcellset(grid, "sliced open X"), getcellset(grid, "sliced open Y"), getcellset(grid, "sliced open Z"))
    delete!(grid.cellsets, "sliced open X")
    delete!(grid.cellsets, "sliced open Y")
    delete!(grid.cellsets, "sliced open Z")

    mesh = prepare_plotable_mesh(grid, cells)
    
    fig = Makie.Figure()

    t = res.t[1]
    tᵒᵇˢ = Makie.Observable(t)
    title = Makie.@lift "solution at t=$( round($(tᵒᵇˢ); sigdigits=4) )"
    Makie.Label(fig[1,1:2], title)

    
    ax  = Makie.Axis3(fig[2,1], aspect=:equal, title="undeformed grid")
    Makie.mesh!(ax, mesh; color=Makie.RGB(0.5,0.5,1.0), shading=Makie.NoShading)
    Makie.wireframe!(ax, mesh; color=:black)
    
    u = evaluate_at_grid_nodes(dh, res.a[1], :u)
    uᵒᵇˢ = Makie.Observable(u)
    deformedmesh = Makie.@lift prepare_plotable_mesh(grid, $(uᵒᵇˢ), cells)
    ax = Makie.Axis3(fig[2,2], aspect=:equal, title="deformed grid")
    Makie.mesh!(ax, deformedmesh; color=Makie.RGB(0.5,0.5,1.0), shading=Makie.NoShading)
    Makie.wireframe!(ax, deformedmesh; color=:black)

    pos = fig[3,1]
    μ = evaluate_at_grid_nodes(dh, res.a[1], :μ)
    μᵒᵇˢ = Makie.Observable(μ)
    ax = Makie.Axis3(pos[1,1], aspect=:equal, title="chemical potential")
    colorrange = Makie.@lift (minimum( $(μᵒᵇˢ) )-1), maximum( $(μᵒᵇˢ) )
    Makie.mesh!(ax, mesh; color=μᵒᵇˢ, colormap=:viridis, colorrange=colorrange, shading=Makie.NoShading)
    Makie.Colorbar(pos[1,2]; colormap=:viridis, colorrange=colorrange)
   
    pos = fig[3,2]
    c = evaluate_at_grid_nodes(dh, res.a[1], :c)
    cᵒᵇˢ = Makie.Observable(c)
    ax = Makie.Axis3(pos[1,1], aspect=:equal, title="concentration")
    colorrange = Makie.@lift (minimum( $(cᵒᵇˢ) )-1), maximum( $(cᵒᵇˢ) )
    Makie.mesh!(ax, mesh; color=cᵒᵇˢ, colormap=:viridis, colorrange=colorrange, shading=Makie.NoShading)
    Makie.Colorbar(pos[1,2]; colormap=:viridis, colorrange=colorrange)    

	file = joinpath("Myresult.mp4")
    #file = joinpath(tempdir(), "Myresult.mp4")
	anim = Makie.record(fig, file, eachindex(res.t); framerate=1) do i
        tᵒᵇˢ[] = res.t[i]
        uᵒᵇˢ[] = evaluate_at_grid_nodes(dh, res.a[i], :u)
        μᵒᵇˢ[] = evaluate_at_grid_nodes(dh, res.a[i], :μ)
        cᵒᵇˢ[] = evaluate_at_grid_nodes(dh, res.a[i], :c)
	end

	return file, fig, anim
end