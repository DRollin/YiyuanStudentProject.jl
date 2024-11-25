"""
    TODO
"""
function get_volume(nodes, ::Ferrite.Triangle)
    a = norm(nodes[1].x .- nodes[2].x)
    b = norm(nodes[1].x .- nodes[3].x)
    c = norm(nodes[2].x .- nodes[3].x)
    s = (a + b + c) / 2
    return sqrt( s*(s-a)*(s-b)*(s-c) )
end

function get_volume(grid::Grid, i::Integer)
    cell = getcells(grid, i)
    return get_volume(getnodes(grid, [cell.nodes...]), cell)
end

function get_volume(grid::Grid, cellset::Set{Int})
    vol = 0.0
    for i in cellset
        vol +=  get_volume(grid, i)
    end
    return vol
end
get_volume(grid::Grid, cellset::String) = get_volume(grid, getcellset(grid, cellset))

"""
    TODO
"""
function get_phase_bias(dh::DofHandler{dim}, cv::CellValues, cellset::Set{Int}, V::Real=get_volume(grid, cellset)) where {dim}
    x̂  = zero(Tensor{1,dim})
    x̂² = zero(Tensor{2,dim})
    for cc in CellIterator(dh, cellset)
        Ferrite.reinit!(cv, cc)
        coords = getcoordinates(cc)
        for qp in 1:getnquadpoints(cv)
            dΩ  = getdetJdV(cv, qp)
            x   = spatial_coordinate(cv, qp, coords)
            x̂  += x * dΩ
            x̂² += (x⊗x) * dΩ
        end
    end
    return x̂/V, x̂²/V
end
