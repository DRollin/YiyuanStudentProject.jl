"""
    TODO
"""
struct Material{dim,T}
    E::Tensor{4,dim,T}
    αᶜʰ::Tensor{2,dim,T}
    k::T
    cʳᵉᶠ::T
    M::Tensor{2,dim,T}
    μʳᵉᶠ::T
end
function Material{dim}(; G::T, K::T, η::T, cʳᵉᶠ::T, μʳᵉᶠ::T, θʳᵉᶠ::T, cᵐ::T, α::T, R::T=8.31446261815324) where {dim, T<:Real}
    E   = SymmetricTensor{4,dim}( (i,j,k,l) -> (K - 2/3*G)*δ(i,j)*δ(k,l) + G*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) )
    αᶜʰ = SymmetricTensor{2,dim}( (i,j) -> δ(i,j)*α )
    M   = SymmetricTensor{2,dim}( (i,j) -> δ(i,j)*η)
    k   = R*θʳᵉᶠ/cᵐ
    return Material{dim,T}(E, αᶜʰ, k, cʳᵉᶠ, M, μʳᵉᶠ)
end

"""
    TODO
"""
struct RVE{dim}
    grid::Grid{dim}
    P::Material{dim}
    M::Material{dim}
end
function RVE(; grid::Grid{dim}, materials::NamedTuple{(:P,:M),Tuple{Material{dim}, Material{dim}}}) where {dim}
    return RVE{dim}(grid, materials.P, materials.M)
end

"""
    TODO
"""
struct LoadCase{dim}
    ε̄::SymmetricTensor{2,dim,Float64}
	μ̄::Float64
    ζ̄::Tensor{1,dim,Float64}
end
function LoadCase{dim}(; ε̄::SymmetricTensor{2,dim,Float64}=zero(SymmetricTensor{2,dim,Float64}), 
                         μ̄::Float64=0.0,
                         ζ̄::Tensor{1,dim,Float64}=zero(Tensor{1,dim,Float64})) where {dim}
    return LoadCase{dim}(ε̄, μ̄, ζ̄)
end
function LoadCase(dim::Int; kwargs...)
	return LoadCase(Val(dim); kwargs... )
end
function LoadCase(::Val{2}; ε̄₁₁::T=0.0, ε̄₁₂::T=0.0, ε̄₂₂::T=0.0,
                            μ̄::T=0.0,
                            ζ̄₁::T=0.0, ζ̄₂::T=0.0) where {T<:Real}
    ε̄ = SymmetricTensor{2,2,Float64}([ε̄₁₁ ε̄₁₂; ε̄₁₂ ε̄₂₂])
    ζ̄ = Tensor{1,2,Float64}([ζ̄₁, ζ̄₂])
    return LoadCase{2}(ε̄, μ̄, ζ̄)
end
function LoadCase(::Val{3}; ε̄₁₁::T=0.0, ε̄₁₂::T=0.0, ε̄₁₃::T=0.0, ε̄₂₃::T=0.0, ε̄₂₂::T=0.0, ε̄₃₃::T=0.0,
                            μ̄::T=0.0,
                            ζ̄₁::T=0.0, ζ̄₂::T=0.0, ζ̄₃::T=0.0) where {T<:Real}
    ε̄ = SymmetricTensor{2,3,Float64}([ε̄₁₁ ε̄₁₂ ε̄₁₃; ε̄₁₂ ε̄₂₂ ε̄₂₃; ε̄₁₃ ε̄₂₃ ε̄₃₃])
    ζ̄ = Tensor{1,3,Float64}([ζ̄₁, ζ̄₂, ζ̄₃])
    return LoadCase{3}(ε̄, μ̄, ζ̄)
end

"""
    TODO
"""
struct PhaseSetup{dim}
    dh::DofHandler
    cells::OrderedSet{Int}
    cv::NamedTuple
    nbf::NamedTuple
    material::Material{dim}
    Kₑ::Matrix{Float64}
    Mₑ::Matrix{Float64}
    submatrices::NamedTuple
end

"""
    TODO
"""
struct RVESetup{dim}
	grid::Grid{dim}
	dh::DofHandler{dim}
	#ch::ConstraintHandler -> Wee need to be able to change the constraints...
    #sets::NamedTuple{(:P,:M),Tuple{Set{Int64},Set{Int64}}}
    phasesetups::NamedTuple{(:P,:M),Tuple{PhaseSetup{dim},PhaseSetup{dim}}}
    K::SparseMatrixCSC{Float64, Int64}
    M::SparseMatrixCSC{Float64, Int64}
    J:: SparseMatrixCSC{Float64, Int64}
    g:: Vector{Float64}
    aⁿ::Vector{Float64}
    aⁿ⁺¹::Vector{Float64}
end
