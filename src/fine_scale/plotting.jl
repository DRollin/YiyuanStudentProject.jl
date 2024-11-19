
function convert_cells(::Grid, cells::Vector{Tetrahedron})
    return [ GeometryBasics.TriangleFace{Int}(facet...) for cell in cells for facet in Ferrite.facets(cell) ]
end
convert_cells(grid::Grid, set::Set{Int}) = convert_cells(grid, getcells(grid, collect(set)))
convert_cells(grid::Grid, set::String) = convert_cells(grid, getcells(grid, set))

function prepare_plotable_mesh(grid::Grid{dim}, set) where {dim}
    nodes = [ GeometryBasics.Point{dim,Float64}(n.x...) for n in grid.nodes ]
    cells = convert_cells(grid, set)
    return GeometryBasics.Mesh(nodes, cells)
end
function prepare_plotable_mesh_1(grid::Grid{dim}, set, displ) where {dim}
    nodes = [ GeometryBasics.Point{dim,Float64}( (n.x .+ u)...) for (n, u) in zip(grid.nodes, displ) ]
    cells = convert_cells(grid, set)
    return GeometryBasics.Mesh(nodes, cells)
end
prepare_plotable_mesh(grid::Grid) = prepare_plotable_mesh(grid, grid.cells)
prepare_plotable_mesh_1(grid::Grid{dim}, u) where {dim} = prepare_plotable_mesh_1(grid, grid.cells, u)

function plot_grid(grid::Grid{3})
    mesh = prepare_plotable_mesh(grid)
    fig = Makie.Figure()
    ax  = Makie.Axis3(fig[1, 1], aspect = :equal, title = "undeformed grid")
    Makie.mesh!(ax, mesh; color=Makie.RGB(0.5,0.5,1.0))
    Makie.wireframe!(ax, mesh; color=:black)
    return fig
end

function plot_result(res, dh, field)
    fig1 = Makie.Figure()
    
    t = res[1][1]
    obsₜ = Makie.Observable(t)
    
    if field == "c"
        ax  = Makie.Axis3(fig1[1, 1], aspect = :equal, title = Makie.@lift "Concentration at t=$( round($(obsₜ);sigdigits=4) )")
        result = evaluate_at_grid_nodes(dh, res[1][3], :c)
    elseif field == "μ"
        ax  = Makie.Axis3(fig1[1, 1], aspect = :equal, title = Makie.@lift "chemical potential at t=$( round($(obsₜ);sigdigits=4) )")
        result = evaluate_at_grid_nodes(dh, res[1][3], :μ)
    end

    u = evaluate_at_grid_nodes(dh, res[1][3], :u)
    obs_result = Makie.Observable(result)
    obs_u = Makie.Observable(u)
    mesh = Makie.@lift prepare_plotable_mesh_1(dh.grid, $(obs_u))

    
	colorrange = Makie.@lift (minimum( $(obs_result) )-1), maximum( $(obs_result) )
	Makie.mesh!(ax, mesh; color=obs_result, colormap=:viridis, colorrange=colorrange)
    Makie.wireframe!(ax, mesh; color=:black)
	Makie.Colorbar(fig1[1,2]; colormap=:viridis, colorrange=colorrange)

	file = joinpath(tempdir(), "Myresult.mp4")
	anim = GLMakie.record(fig1, file, eachindex(res); framerate=1) do i
        if field == "c"
		    obsₜ[], obs_result[] = res[i][1], evaluate_at_grid_nodes(dh, res[i][3], :c)
        elseif field == "μ"
            obsₜ[], obs_result[] = res[i][1], evaluate_at_grid_nodes(dh, res[i][3], :μ)
        end
	end

    
	return file, fig1, anim
end