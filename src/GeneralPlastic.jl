struct GeneralPlastic{D, T} <: AbstractMaterial 
    C_elas :: SymFourthOrderTensorValue{D, T}
    yieldStress :: Fumction
    yieldSurface :: Function
end

struct GeneralPlasticState{D, T} <: AbstractMaterialState
    λ :: T
    εp :: SymTensorValue{D, T}
end

initial_material_state(::GeneralPlastic{D, T}) where{D, T} = GeneralPlasticState(T(0.0), zero(SymTensorValue{D, T}))

struct GeneralPlasticCache{T<:NLsolve.onceDifferentiable} <: AbstractCache
    nlsolve_cache :: T
end

get_n_scalar_equations(::GeneralPlastic{D, T}) where{D,T} = D*(D+1)//2 + 1

# We don't need ResidualGeneralPlastic

function StateToVector(v::SVector{4, T}, r::GeneralPlasticState{2, T}) where{T}
    # TODO check vector length
    v[1] .= r.λ
    view(v, 2:4) .= r.εp.data 
end

function StateToVector(v::SVector{7, T}, r::GeneralPlasticState{3, T}) where{T}
    # TODO check vector length
    v[1] .= r.λ
    view(v, 2:7) .= r.εp.data 
end

function VectorToState(::Type{GeneralPlasticState{2, T}}, v::SVector{4, T})
    λ = v[1]
    εp = SymTensorValue{2, T, 3}(v[2], v[3], v[4])
    return GeneralPlasticState{2, T}(λ, εp)
end

function VectorToState(::Type{GeneralPlasticState{3, T}}, v::SVector{7, T})
    λ = v[1]
    εp = SymTensorValue{3, T, 6}(v[2], v[3], v[4], v[5], v[6], v[7])
    return GeneralPlasticState{3, T}(λ, εp)
end

mutable struct GenerapPlasticNLOP{D, T} <: NonlinearOperator 
    material :: GeneralPlastic{D, T}
    state :: GeneralPlasticState{D, T}
    ε :: SymTensorValue{D, T}
end

function zero_initial_guess(op::GenerapPlasticNLOP{D, T}) where{T}
    x = allocate_residual(op, T[])
    fill!(x, zero(eltype(x)))
    x
end

function allocate_residual(op::GenerapPlasticNLOP{D, T}, x::AbstractVector{T}) where{T}
    similar(x, T)
end

function material_response(
    material::GeneralPlastic, ε::SymTensorValue{D, T}, state::GeneralPlasticState{D, T}, Δt, cache, extras) where {D, T}

    σ_trial = material.C_elas ⊙ (ε - state.εp)

    Φ = material.yieldSurface(σ_trial) - material.yieldStress(state.λ)

    if Φ <= 0
        return σ_trial, material.C_elas, GeneralPlasticState(state.εp, state.λ)
    else
        op = GenerapPlasticNLOP
        nls = NLSolver(show_trace=false, method=:newton)
        x0 = zero_initial_guess(op)
        solve!(x0, nls, op)

        newState = VectorToState(state, x0)
        ∂f∂σ = grad_wrt_entries(op.material.yieldSurface, σ)
    end

end

function residual!(r::AbstractVector{T}, op::GenerapPlasticNLOP{D,T}, x::AbstractVector{T}) where{T}
    state = VectorToState(GeneralPlasticState{D,T}, x)
    σ = op.material.C_elas ⊙ (op.ε - state.εp)
    ∂f∂σ = grad_wrt_entries(op.material.yieldSurface, σ)
    Rεp = state.εp - op.state.εp - (state.λ - op.state.λ)*∂f∂σ
    r[1] .= op.material.yieldSurface(σ) - op.material.yieldStress(state.λ)
    r[2:end] .= Rεp.data
end
