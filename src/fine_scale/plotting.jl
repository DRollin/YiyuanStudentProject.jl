# Following code is a draft, not tested yet!

function convert_cells(::Grid, cells::Vector{Tetrahedron})
    return [ GeometryBasics.TriangleFace{Int}(facet...) for cell in cells for facet in facets(cell) ]
end
convert_cells(grid::Grid, set::Set{Int}) = convert_cells(grid, getcells(grid, collect(set)))
convert_cells(grid::Grid, set::String) = convert_cells(grid, getcells(grid, set))

function prepare_plotable_mesh(grid::Grid{dim}, set) where {dim}
    nodes = [ GeometryBasics.Point{dim,Float64}(n.x...) for n in grid.nodes ]
    cells = convert_cells(grid, set)
    return GeometryBasics.Mesh(nodes, cells)
end
function prepare_plotable_mesh(grid::Grid{dim}, set, displ::Vector{Vec{dim}}) where {dim}
    nodes = [ GeometryBasics.Point{dim,Float64}( (n.x + u)...) for (n, u) in zip(grid.nodes, displ) ]
    cells = convert_cells(grid, set)
    return GeometryBasics.Mesh(nodes, cells)
end
prepare_plotable_mesh(grid::Grid) = prepare_plotable_mesh(grid, grid.cells)
prepare_plotable_mesh(grid::Grid{dim}, u::Vector{Vec{dim}}) where {dim} = prepare_plotable_mesh(grid, grid.cells, u)

function plot_grid(grid::Grid{3})
    mesh = prepare_plotable_mesh(grid)
    fig = Makie.Figure()
    ax  = Makie.Axis3()
    Makie.mesh!(ax, mesh; color=Makie.RGB(0.5,0.5,1.0))
    Makie.wireframe!(ax, mesh; color=:black)
    return fig
end

function plot_c(grid::Grid{3}, dh::DofHandler, u::Vector{Float64})
    obs = Makie.Observable(u)
    c = @lift evaluate_at_grid_nodes(dh, $obs, :c)
    u = @lift evaluate_at_grid_nodes(dh, $obs, :u)
    mesh = @lift prepare_plotable_mesh(grid, $u)
    fig = Makie.Figure()
    ax  = Makie.Axis3()
    
    Makie.mesh!(ax, mesh; color=c)
    Makie.wireframe!(ax, mesh; color=:black)
    return fig, obs
end