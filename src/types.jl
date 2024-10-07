#fine_scale
struct iso_pv_Material{dim}
    E::Tensor{4,dim,Float64} #4th order elasticity modulus tensor
    Κ::Tensor{2,dim,Float64} #2nd order constant permeability tensor
    α::Float64 #Biot coefficient
    β::Float64 #intrinsic compliance of the pore fluid
end

function elastic_stiffness(G, K)
    I2 = one(SymmetricTensor{2,3})
    I4sym = symmetric(otimesu(I2,I2))
    I4vol = I2⊗I2/3
    return 2*G*(I4sym-I4vol) + 3*K*I4vol
end

function Material_pre(G::Number, K::Number, κ::Number, α::Number, β::Number, dim)
    return iso_pv_Material{dim}(elastic_stiffness(G, K), Tensor{2,dim,Float64}( (i,j) -> κ ), α, β)
end

struct RVEProblem{dim}
    grid::Grid{dim}
    P::iso_pv_Material{dim}
    M::iso_pv_Material{dim}
end

function RVEProblem(; grid::Grid{dim}, materials::NamedTuple{(:P,:M),Tuple{iso_pv_Material{dim}, iso_pv_Material{dim}}}) where {dim}
    return RVEProblem{dim}(grid, materials.P, materials.M)
end

struct LoadCase{dim}
	μ̄::Float64
    ζ̄::Tensor{1,dim,Float64}
end

function LoadCase_pre(dim, μ̄, ζ̄1, ζ̄2, ζ̄3)
	return LoadCase{dim}( μ̄, Tensor{1,dim,Float64}(( ζ̄1, ζ̄2, ζ̄3)) )
end

struct iso_pv_ElementSetup{dim}
    cells::DofHandler
    cv::NamedTuple
    nbf::NamedTuple
    material::iso_pv_Material{dim}
end


struct FESetup_base{dim}
	grid::Grid{dim}
	dh::DofHandler
	ch::ConstraintHandler
    Load::LoadCase{dim}
    sets::NamedTuple{(:P,:M),Tuple{Set{Int64},Set{Int64}}}
    setups::NamedTuple{(:P,:M),Tuple{iso_pv_ElementSetup{dim},iso_pv_ElementSetup{dim}}}
end


struct pv_Problem{dim}
    setup::FESetup_base{dim}
    K::SparseArrays.SparseMatrixCSC{Float64, Int64}
    f::Vector{Float64}
    a::Vector{Float64}
    a_old::Vector{Float64}
end


#upscaling
struct EffectiveResponse{dim,T}
    c̄::T
    c̄₂::Tensor{1,dim,T}
    j̄::Tensor{1,dim,T}
end

function Base.:-(e₁::EffectiveResponse{dim,T}, e₂::EffectiveResponse{dim,T}) where {dim,T}
    res = collect( getproperty(e₁,s) - getproperty(e₂,s) for s in fieldnames(EffectiveResponse) )
    return ( EffectiveResponse{dim,T}(res...) )
end

struct RVEResponses{dim,T}
    l::EffectiveResponse{dim,T}
    ū::EffectiveResponse{dim,T}
    ∇ū::Union{NTuple{dim,EffectiveResponse{dim,T}}, Missing}
end
function RVEResponses(l::ER, ū::ER, ∇ū::Union{NTuple{dim,ER},Missing}) where{T,dim,ER<:EffectiveResponse{dim,T}}
    return RVEResponses{dim,T}(l, ū, ∇ū)
end

struct SensitivitiesForScalar{dim,T}
    l::T
    ū::T
    ∇ū::Union{Tensor{1,dim,T}, Missing}
end
function SensitivitiesForScalar(res::RVEResponses{dim,T}, field::Symbol) where {dim,T}
    l  = getproperty(res.l, field)
    ū  = getproperty(res.ū, field)
    ∇ū = Tensor{1,dim,T}( [getproperty.(res.∇ū, field)...]  )
    return SensitivitiesForScalar{dim,T}(l, ū, ∇ū)
end

struct SensitivitiesForVector{dim,T}
    l::Tensor{1,dim,T}
    ū::Tensor{1,dim,T}
    ∇ū::Union{Tensor{2,dim,T}, Missing}
end
function SensitivitiesForVector(res::RVEResponses{dim,T}, field::Symbol) where {dim,T}
    l  = getproperty(res.l, field)
    ū  = getproperty(res.ū, field)
    ∇ū = Tensor{2,dim,T}( hcat(getproperty.(res.∇ū, field)...)  )
    return SensitivitiesForVector{dim,T}(l, ū, ∇ū)
end

struct Sensitivities{dim,T}
    c̄̂::SensitivitiesForScalar{dim,T}
    c̄̂₂::SensitivitiesForVector{dim,T}
    j̄̂::SensitivitiesForVector{dim,T}
end
function Sensitivities(res::RVEResponses{dim,T}) where {dim,T}
    c̄̂  = SensitivitiesForScalar(res, :c̄)
    c̄̂₂ = SensitivitiesForVector(res, :c̄₂)
    j̄̂  = SensitivitiesForVector(res, :j̄)
    return Sensitivities{dim,T}(c̄̂, c̄̂₂, j̄̂)
end


#
