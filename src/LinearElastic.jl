
struct LinearElastic{T, D} <: AbstractMaterial
    E :: T
    ν :: T
    C_elas :: SymFourthOrderTensorValue{D, T}
    LinearElastic(D::Int, E::T, ν::T) where {T} = new{T, D}(E, ν, elastic_tangent(D, E, ν))
end

# a constructor with keyword arguments
LinearElastic(; D::Int, E::T, ν::T) where {T} = LinearElastic(D, E, ν)

struct LinearElasticState{T, D} <: AbstractMaterialState end

initial_material_state(::LinearElastic{T, D}) where{T, D} = LinearElasticState{T, D}()

get_stress_type(::LinearElasticState{T, D}) where{T, D} = SymTensorValue{D, T}

function material_response(m::LinearElastic{T, D}, ε::SymTensorValue{D, T},
    state::LinearElasticState{T, D}=LinearElasticState{T, D}(), Δt=nothing; cache=nothing, options=nothing) where{T, D}

    σ = m.C_elas ⊙ ε

    return σ, m.C_elas, state
end
