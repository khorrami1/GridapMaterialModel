
struct LinearElastic{D, T} <: AbstractMaterial
    E :: T
    ν :: T
    C_elas :: SymFourthOrderTensorValue{D, T}
    LinearElastic{D, T}(E::T, ν::T) where {D, T} = new{D, T}(E, ν, elastic_tangent(D, E, ν))
end

# a constructor with keyword arguments
# LinearElastic(; D::Int, E::T, ν::T) where {T} = LinearElastic{D, T}(D, T(E), T(ν))

struct LinearElasticState{D, T} <: AbstractMaterialState end

initial_material_state(::LinearElastic{D, T}) where{D, T} = LinearElasticState{D, T}()

get_stress_type(::LinearElasticState{D, T}) where{D, T} = SymTensorValue{D, T}

function material_response(m::LinearElastic{D, T}, ε::SymTensorValue{D, T},
    state::LinearElasticState{D, T}=LinearElasticState{T, D}(), Δt=nothing; cache=nothing, options=nothing) where{T, D}

    σ = m.C_elas ⊙ ε

    return σ, m.C_elas, state
end
