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
    
    return grid
end