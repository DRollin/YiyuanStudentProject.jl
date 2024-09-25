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

function Material_pre{dim}(; G::Number, K::Number, κ::Number, α::Number, β::Number) where {dim}
    return iso_pv_Material{dim}(elastic_stiffness(G, K), Tensor{2,dim,Float64}( (i,j) -> κ ), α, β)
end

struct RVEProblem{dim}
    grid::Grid{dim}
    P::BulkMaterial{dim}
    M::BulkMaterial{dim}
end

function RVEProblem(; grid::Grid{dim}, materials::NamedTuple{(:P,:M),Tuple{BulkMaterial{dim},BulkMaterial{dim}}}) where {dim}
    return RVEProblem{dim}(grid, materials.P, materials.M)
end

struct LoadCase{dim}
	μ̄::Float64
    ζ̄::Tensor{1,dim,Float64}
end
function LoadCase{2}(; μ̄=0.0, ζ̄₁=0.0, ζ̄₂=0.0)
	return LoadCase{2}( μ̄, Tensor{1,2,Float64}(( ζ̄₁, ζ̄₂)) )
end

struct iso_pv_ElementSetup{dim}
    cells::DofHandler
    cv::CellValues
    nbf::Int
    material::iso_pv_Material{dim}
    #V::Float64
end


struct FESetup{dim}
	grid::Grid{dim}
	dh::DofHandler
	ch::ConstraintHandler
    Load::LoadCase{dim}
    sets::NamedTuple{(:P,:M),Tuple{Set{Int64},Set{Int64}}}
    setups::NamedTuple{(:P,:M),Tuple{BulkElementSetup{2},BulkElementSetup{2}}}
end


struct ScratchValues{T, TT<:AbstractTensor, Ti, dim}
    Me::Matrix{T}
    Ke::Matrix{T}
    cv::CellValues
    global_dofs::Vector{Int}
    ∇u::Vector{TT}
    coordinates::Vector{Vec{dim,T}}
    assembler_K::Ferrite.AssemblerSparsityPattern{T,Ti}
    assembler_M::Ferrite.AssemblerSparsityPattern{T,Ti}
end