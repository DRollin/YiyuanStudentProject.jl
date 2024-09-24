function convert_nodes(grid::Grid{dim}) where {dim}
    return collect( GeometryBasics.Point{dim,Float64}(n.x...) for n in grid.nodes )
end

function convert_cells(::Grid, cells::Vector{<:Ferrite.AbstractCell})
    bulkcells = filter(c -> !(c isa InterfaceCell), cells)
    return collect( GeometryBasics.TriangleFace{Int}(cell.nodes...) for cell in bulkcells )
end
convert_cells(grid::Grid, set::Set{Int}) = convert_cells(grid, getcells(grid, collect(set)))
convert_cells(grid::Grid, set::String) = convert_cells(grid, getcells(grid, set))

function prepare_plotable_mesh(grid::Grid, set)
    return GeometryBasics.Mesh(convert_nodes(grid), convert_cells(grid, set))
end
prepare_plotable_mesh(grid::Grid) = prepare_plotable_mesh(grid, grid.cells)

function plot_mesh_overlay!(ax::Makie.Axis, grid::Grid; matrix=true, particles=true)
    if matrix
        mesh = prepare_plotable_mesh(grid, "matrix")
        Makie.wireframe!(ax, mesh; color=:gray)
    end
    if particles
        mesh = prepare_plotable_mesh(grid, "particles")
        Makie.wireframe!(ax, mesh; color=:black)
    end
    return ax
end

function prepare_plotable_lines(grid::Grid, cellset::Set{Int})
    lines = Set{Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}}}()
    nodes = grid.nodes
    for cellid in cellset
        n = getcells(grid, cellid).nodes
        push!(lines, (Tuple(nodes[n[1]].x), Tuple(nodes[n[2]].x)))
    end
    return collect(lines)
end
prepare_plotable_lines(grid::Grid, set::String) = prepare_plotable_lines(grid, getcellset(grid, set))

function plot_grid(grid::Grid{2})
    mesh = prepare_plotable_mesh(grid)
    fig, ax, p = Makie.mesh(mesh; color=Makie.RGB(0.5,0.5,1.0))
    Makie.wireframe!(ax, mesh; color=:black)
    return fig, ax, p
end

function plot_rve_grid(grid::Grid{2})
    mesh = prepare_plotable_mesh(grid, "matrix")
    fig, ax, p = Makie.mesh(mesh; color=Makie.RGB(0.5,0.5,1.0))
    Makie.wireframe!(ax, mesh; color=:black)

    mesh = prepare_plotable_mesh(grid, "particles")
    Makie.mesh!(ax, mesh; color=Makie.RGB(1.0,0.5,0.5))
    Makie.wireframe!(ax, mesh; color=:black)
    return fig, ax, p
end

function get_potential(dh::DofHandler, u::Vector, cellset::Set{Int})
    grid = dh.grid
    uᵢ = zeros(length(grid.nodes))
    for cc in CellIterator(dh, cellset)
        cell = getcells(grid, cc.cellid.x)
        dofs = celldofs(cc)
        nodes = [cell.nodes...]
        uᵢ[nodes] .= u[dofs]
    end
    return uᵢ
end

function plot_potential(u::Vector, setup::FESetup)
    (; grid, dh) = setup
    bulkcells = union(getcellset(grid, "matrix"), getcellset(grid, "particles"))
    mesh = prepare_plotable_mesh(grid, bulkcells)
    coloring = get_potential(dh, u, bulkcells)

    fig = Makie.Figure()
    ax = Makie.Axis(fig[1,1]; title="Chemical Potential")
    if minimum(coloring) ≈ maximum(coloring)
        p = Makie.mesh!(ax, mesh; color=coloring, colormap=:heat, colorrange=(0.0,1.0))
    else
        p = Makie.mesh!(ax, mesh; color=coloring, colormap=:heat)
    end
    Makie.Colorbar(fig[1, 2], p)
    return fig, ax, p
end

function plot_dns_potential(sol::Vector, dh::DofHandler, t::Real)
    grid = dh.grid
    obsu = Makie.Observable(sol)
    obst = Makie.Observable(t)

    bulkcells = union(getcellset(grid, "matrix"), getcellset(grid, "particles"))
    u = Makie.@lift get_potential(dh, $obsu, bulkcells)
    mesh = prepare_plotable_mesh(grid)

    title = Makie.@lift "Potentials at t=" * string( round($obst; sigdigits=3) ) * " s"
    
    colorsettings = (colormap=:heat, colorrange=(-0.01, 1.01), lowclip=:blue, highclip=:red)
    axissettings = (height=100, width=600)
    
    fig = Makie.Figure(; size=(800,300))
    Makie.Label(fig[1,1:2], title; valign=:bottom, font=:bold)
    ax = Makie.Axis(fig[2,1]; title=L"\mu", axissettings...)
    
    Makie.mesh!(ax, mesh; color=u, colorsettings...)
    Makie.Colorbar(fig[2,2]; colorsettings...)
    
    return fig, ax, obsu, obst
end
function plot_dns_potential(dh, res::Dict)
    t, u = select_state(res, 0.0)
    return plot_dns_potential(u.u, dh, t)
end

function select_state(res::Dict, t::Real)
    t = argmin(x -> abs(x-t), keys(res))
    return t, res[t]
end