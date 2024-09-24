#############################################################################
# Sphere generation
#############################################################################

function generate_spheres(ϕ::Real, pdf::Uniform, domainsize::NTuple{dim,Real}) where {dim}
    if ϕ == 0
        return BubbleBath.Sphere{2}[]
    else
        return bubblebath(pdf, ϕ, domainsize)
    end
end
generate_spheres(ϕ::Real, d::Real, domainsize::NTuple{dim,Real}) where {dim} = generate_spheres(ϕ, get_radius_pdf(d), domainsize)
get_radius_pdf(d::Real) = Uniform(0.9*d, 1.1*d)

function generate_rve_spheres(; ϕ::Real, d::Real, dx::NTuple{dim,<:Real}) where {dim}
    spheres = generate_spheres(ϕ, d, dx)
    return shift_by.(spheres, (-0.5 .* dx,))
end
shift_by(sphere::BubbleBath.Sphere{dim}, shift::NTuple{dim,T}) where {dim,T} = BubbleBath.Sphere(sphere.pos .+ shift, sphere.radius)

#############################################################################
# Meshing
#############################################################################

function get_tags_from_dimtags(v::Vector{Tuple{T,T}}) where {T<:Integer}
    return collect( t[2] for t in v )
end

function add_sphere_to_gmsh(s::BubbleBath.Sphere{2})
    return (2, gmsh.model.occ.addDisk(s.pos..., 0.0, s.radius, s.radius))
end

function generate_box_grid(x₀::NTuple{2,Real}, xₑ::NTuple{2,Real}, ϕ::Real, d::Real, meshsize::Real)
    dx = xₑ .- x₀
    spheres = generate_spheres(ϕ, d, dx)
    spheres = shift_by.(spheres, (x₀,))
    return generate_box_grid(x₀, xₑ, spheres, meshsize)
end
function generate_box_grid(x₀::NTuple{2,Real}, xₑ::NTuple{2,Real}, spheres::Vector{BubbleBath.Sphere{2}}, meshsize::Real)
    dim = 2
    dx = xₑ .- x₀
    @assert all( x₀ .< xₑ ) "Lower bounds ($(x₀)) must be smaller than upper bounds ($(xₑ))!"
    nel = round.((Int,), dx ./ meshsize)
    x = Vec{2, Float64}(x₀)
    y = Vec{2, Float64}(xₑ)
    grid = generate_grid(Triangle, nel , x, y)
    #round.((Int,), dx ./ meshsize)
    function check_if_in_particle(x)
        for s in spheres
            if norm( s.pos .- x ) ≤ s.radius
                return true
            end
        end
        return false
    end
    addcellset!(grid, "particles", check_if_in_particle)
    particles = getcellset(grid, "particles")
    allcells = Set( 1:length(grid.cells))
    addcellset!(grid, "matrix", setdiff(allcells, particles))
    #for parallel setting
    colors = create_coloring(grid)
    return grid, colors
end


#############################################################################
# RVE grid
#############################################################################

generate_rve_grid(; ϕ::Real, d::Real, meshsize::Real, dx::NTuple{2,Real}=(1.0,1.0)) = generate_box_grid(-0.5 .* dx, 0.5 .* dx, ϕ, d, meshsize)
generate_rve_grid(spheres::Vector{BubbleBath.Sphere{2}}; meshsize::Real, dx::NTuple{2,Real}=(1.0,1.0)) = generate_box_grid(-0.5 .* dx, 0.5 .* dx, spheres, meshsize)

#############################################################################
# DNS grid
#############################################################################
#=
function generate_dns_grid(; nᵣᵥₑ::Int, gridᵣᵥₑ::Grid{dim}, lᵣᵥₑ::Real) where {dim}
    basegrid = deepcopy(gridᵣᵥₑ)
    shift_by!(basegrid, Vec{dim}(ones(dim).*(lᵣᵥₑ/2)))
    grid  = deepcopy(basegrid)
    newset = Set{Int}(1:length(grid.cells))
    merge!(grid.cellsets, Dict{String,Set{Int}}("RVE 1" => newset))
    shift = zeros(dim); shift[1] = lᵣᵥₑ; shift = Vec{dim}(shift)
    x̄ = zeros(nᵣᵥₑ)
    x̄[1] = lᵣᵥₑ/2
    for i in 2:nᵣᵥₑ
        x̄[i] = x̄[i-1] + lᵣᵥₑ
        shift_by!(basegrid, shift)
        merge_grids!(grid, basegrid, "right", "left", "RVE $i")
    end
    return grid, x̄
end

function scale_to_origin!(grid::Grid, scale::Real)
    for (i, n) in enumerate(grid.nodes)
        s = scale*n.x
        grid.nodes[i] = shift_by(n, s)
    end
    return grid
end

function shift_by!(grid::Grid{dim}, s::Vec{dim}) where {dim}
    for (i, n) in enumerate(grid.nodes)
        grid.nodes[i] = shift_by(n, s)
    end
    return grid
end
function shift_by(n::Node{dim}, s::Vec{dim}) where {dim}
    return Node(n.x + s)
end

function merge_grids!(grid₁::Grid, grid₂::Grid, faceset₁::String, faceset₂::String, newsetname::String)
    mapping = get_face_node_mapping(grid₁, grid₂, getfaceset(grid₁, faceset₁), getfaceset(grid₂, faceset₂))
    nodeoffset = length(grid₁.nodes)
    celloffset = length(grid₁.cells)
    newset = Set{Int}(celloffset+1:celloffset+length(grid₂.cells))
    merge!(grid₁.cellsets, Dict{String,Set{Int}}(newsetname => newset))
    append!( grid₁.nodes, grid₂.nodes )
    append!( grid₁.cells, adapt_cell.(grid₂.cells; nodeoffset=nodeoffset, mapping=mapping))

    for (name, set) in grid₂.cellsets
        adaptedset = set .+ celloffset
        if name in keys(grid₁.cellsets)
            union!(grid₁.cellsets[name], adaptedset)
        else
            merge!(grid₁.cellsets, Dict{String,Set{Int}}(name => set))
        end
    end

    for (name, set) in grid₂.facesets
        adaptedset = Set{FaceIndex}(adapt_face_index.(set; celloffset=celloffset))
        if name == faceset₁
            grid₁.facesets[name] = adaptedset
        elseif name == faceset₂
            continue
        elseif name in keys(grid₁.facesets)
            union!(grid₁.facesets[name], adaptedset)
        else
            merge!(grid₁.facesets, Dict{String,Set{FaceIndex}}(name => set))
        end
    end
    return grid₁
end

function adapt_face_index(i::FaceIndex; celloffset::Int)
    return FaceIndex((i[1]+celloffset, i[2]))
end

function adapt_cell(cell::C; nodeoffset::Int, mapping::Dict{Int,Int}) where {C<:Ferrite.AbstractCell}
    nodes = [cell.nodes...]
    for (i, n) in enumerate(nodes)
        if n in keys(mapping)
            nodes[i] = mapping[n]
            continue
        end
        nodes[i] += nodeoffset
    end
    return C((nodes...,))
end
function adapt_cell(cell::C; nodeoffset::Int, mapping::Dict{Int,Int}) where {C<:InterfaceCell}
    return C(adapt_cell(cell.here;  nodeoffset=nodeoffset, mapping=mapping), 
             adapt_cell(cell.there; nodeoffset=nodeoffset, mapping=mapping))
end

function get_face_node_mapping(grid₁::Grid, grid₂::Grid, faceset₁::Set{FaceIndex}, faceset₂::Set{FaceIndex})
    nodeids₁ = Set(vcat(collect( begin
            cell = getcells(grid₁, cellid)
            [Ferrite.faces(cell)[face]...,]
        end for (cellid, face) in faceset₁  )...))
    nodeids₂ = Set(vcat(collect( begin
            cell = getcells(grid₂, cellid)
            [Ferrite.faces(cell)[face]...,]
        end for (cellid, face) in faceset₂  )...))
    mapping = Dict{Int,Int}()
    for n₁ in nodeids₁
        foundmatch = false
        for n₂ in nodeids₂
            v₁ = getnodes(grid₁, n₁).x
            v₂ = getnodes(grid₂, n₂).x
            if norm( v₁ - v₂ ) / norm( v₁ ) ≤ 1e-8
                foundmatch = true
                merge!(mapping, Dict(n₂ => n₁))
                delete!(nodeids₂, n₂)
                break
            end
        end
        if ! foundmatch
            throw(ErrorException("No matching node found!"))
        end
    end
    return mapping
end
=#