struct GeneralPlastic{D, T} <: AbstractMaterial 
    C_elas :: SymFourthOrderTensorValue{D, T}
    yieldStress :: Function
    yieldSurface :: Function
end

struct GeneralPlasticState{D, T} <: AbstractMaterialState
    λ :: T
    εp :: SymTensorValue{D, T}
end

initial_material_state(::GeneralPlastic{D, T}) where{D, T} = GeneralPlasticState(T(0.0), zero(SymTensorValue{D, T}))

# struct GeneralPlasticCache{T<:NLsolve.onceDifferentiable} <: AbstractCache
#     nlsolve_cache :: T
# end

# get_n_scalar_equations(::GeneralPlastic{D, T}) where{D,T} = D*(D+1)//2 + 1

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

function VectorToState(::Type{GeneralPlasticState{2, T}}, v::SVector{4, T}) where{T}
    λ = v[1]
    εp = SymTensorValue{2, T, 3}(v[2], v[3], v[4])
    return GeneralPlasticState{2, T}(λ, εp)
end

function VectorToState(::Type{GeneralPlasticState{3, T}}, v::SVector{7, T}) where{T}
    λ = v[1]
    εp = SymTensorValue{3, T, 6}(v[2], v[3], v[4], v[5], v[6], v[7])
    return GeneralPlasticState{3, T}(λ, εp)
end

mutable struct GenerapPlasticNLOP{D, T} <: NonlinearOperator 
    material :: GeneralPlastic{D, T}
    state :: GeneralPlasticState{D, T}
    ε :: SymTensorValue{D, T}
end

function zero_initial_guess(op::GenerapPlasticNLOP{D, T}) where{D, T}
    x = allocate_residual(op, T[])
    fill!(x, zero(eltype(x)))
    x
end

function allocate_residual(op::GenerapPlasticNLOP{D, T}, x::AbstractVector{T}) where{D, T}
    similar(x, T)
end

function material_response(op::GenerapPlasticNLOP{D,T}, Δt, cache, extras; get_Cep=true) where {D, T}

    σ_trial = op.material.C_elas ⊙ (ε - op.state.εp)

    Φ = op.material.yieldSurface(σ_trial) - op.material.yieldStress(op.state.λ)

    if Φ <= 0
        return σ_trial, op.material.C_elas, GeneralPlasticState(op.state.εp, op.state.λ)
    else
        nls = NLSolver(show_trace=false, method=:newton)
        x0 = zero_initial_guess(op)
        StateToVector(x0, op.state)
        solve!(x0, nls, op)

        newState = VectorToState(state, x0)
        σ = op.material.C_elas ⊙ (op.ε - newState.εp)
        if (get_Cep)
            ∂f∂σ = grad_wrt_entries(op.material.yieldSurface, σ)
            H = grad_wrt_entries(op.material.yieldStress, newState.λ)
            C_f_σ = op.material.C_elas ⊙ ∂f∂σ
            f_σ_C = ∂f∂σ ⊙ op.material.C_elas
            C_ep = op.material.C_elas - (C_f_σ ⊗ f_σ_C)/(H + ∂f∂σ ⊙ C_f_σ)
            return σ, C_ep, newState
        else
            return σ, newState
        end
    end

end

function residual!(r::AbstractVector{T}, op::GenerapPlasticNLOP{D,T}, x::AbstractVector{T}) where{D, T}
    state = VectorToState(GeneralPlasticState{D,T}, x)
    σ = op.material.C_elas ⊙ (op.ε - state.εp)
    ∂f∂σ = grad_wrt_entries(op.material.yieldSurface, σ)
    Rεp = state.εp - op.state.εp - (state.λ - op.state.λ)*∂f∂σ
    r[1] .= op.material.yieldSurface(σ) - op.material.yieldStress(state.λ)
    r[2:end] .= Rεp.data
end
