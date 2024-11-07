#fine_scale

struct iso_cm_Material{dim}
    E::Tensor{4,dim,Float64} 
    α_ch::Tensor{2,dim,Float64} 
    k::Float64 
    c_ref::Float64 
    M::Tensor{2,dim,Float64}
    μ_ref::Float64 
end

"""
    elastic_stiffness(G, K)

Return the elastic stiffness tensor considering the given `G` the Shear modulus and `K` the Bulk modulus
"""
function elastic_stiffness(G, K)
    I2 = one(SymmetricTensor{2,3})
    I4sym = symmetric(otimesu(I2,I2))
    I4vol = I2⊗I2/3
    return 2*G*(I4sym-I4vol) + 3*K*I4vol
end



"""
    coefficient_k(R, θ_ref, c_m)

Return the concentration-chemical potantial coefficient considering the given `R` the ideal gas constant,
 `θ_ref` the reference temperature and `c_m` the ion concentration 
"""
function coefficient_k(R, θ_ref, c_m)
    return R*θ_ref/c_m
end


 
"""
    Material_pre(G, K, α, R, θ_ref, c_m, c_ref, η, dim)

Return a `iso_cm_Material` material struct with the following material parameters:

    E:      fourth order stiffness tensor E
    α_ch:   isotropic ion intercalation tensor
    k:      concentration-chemical potantial coefficient
    c_ref:  reference concentration
    M:      mobility tensor
    
function `elastic_stiffness` is called with the inputs `G` the Shear modulus and `K` the Bulk modulus to compute the linear elastic fourth stiffness tensor `E`
function `coefficient_k` is called with the inputs `R` the ideal gas constant, `θ_ref` the reference temperature and `c_m` the ion concentration to compute 
the concentration-chemical potantial coefficient `k`
"""
function Material_pre(G, K, α, R, θ_ref, c_m, c_ref, η, μ_ref, dim)
    return iso_cm_Material{dim}(elastic_stiffness(G, K), Tensor{2,dim,Float64}( (i,j) -> δ(i,j)*α ), coefficient_k(R, θ_ref, c_m), c_ref, Tensor{2,dim,Float64}( (i,j) -> δ(i,j)* η), μ_ref)
end


struct RVEProblem{dim}
    grid::Grid{dim}
    P::iso_cm_Material{dim}
    M::iso_cm_Material{dim}
end

function RVEProblem(; grid::Grid{dim}, materials::NamedTuple{(:P,:M),Tuple{iso_cm_Material{dim}, iso_cm_Material{dim}}}) where {dim}
    return RVEProblem{dim}(grid, materials.P, materials.M)
end

struct LoadCase{dim}
	μ̄::Float64
    ζ̄::Tensor{1,dim,Float64}
end

function LoadCase_pre(dim, μ̄, ζ̄1, ζ̄2, ζ̄3)
	return LoadCase{dim}( μ̄, Tensor{1,dim,Float64}(( ζ̄1, ζ̄2, ζ̄3)) )
end

struct iso_cm_ElementSetup{dim}
    cells::DofHandler
    cv::NamedTuple
    nbf::NamedTuple
    material::iso_cm_Material{dim}
end


struct FESetup_base{dim}
	grid::Grid{dim}
	dh::DofHandler
	ch::ConstraintHandler
    Load::LoadCase{dim}
    sets::NamedTuple{(:P,:M),Tuple{Set{Int64},Set{Int64}}}
    setups::NamedTuple{(:P,:M),Tuple{iso_cm_ElementSetup{dim},iso_cm_ElementSetup{dim}}}
    J:: SparseArrays.SparseMatrixCSC{Float64, Int64}
    g:: Vector{Float64}
end


struct cm_Problem{dim}
    setup::FESetup_base{dim}
    K::SparseArrays.SparseMatrixCSC{Float64, Int64}
    M::SparseArrays.SparseMatrixCSC{Float64, Int64}
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
