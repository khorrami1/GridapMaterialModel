# using Test
# using LinearAlgebra
# using StaticArrays
# using ForwardDiff
# using Gridap
# using Gridap.Algebra

# ---------------------------------------------------------------------
# We use fully qualified Gridap.Algebra methods ("Option B") because
# this is the safest and avoids namespace/interface extension issues.
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# Material definition
# ---------------------------------------------------------------------

struct GeneralPlastic{D,T,YS,YF} <: AbstractMaterial
    C_elas::SymFourthOrderTensorValue{D,T}
    yieldStress::YS
    yieldSurface::YF
end

GeneralPlastic(
    C_elas::SymFourthOrderTensorValue{D,T},
    yieldStress,
    yieldSurface
) where {D,T} =
    GeneralPlastic{D,T,typeof(yieldStress),typeof(yieldSurface)}(
        C_elas, yieldStress, yieldSurface
    )

struct GeneralPlasticState{D,T} <: AbstractMaterialState
    λ::T
    εp::SymTensorValue{D,T}
end

initial_material_state(::GeneralPlastic{D,T}) where {D,T} =
    GeneralPlasticState{D,T}(zero(T), zero(SymTensorValue{D,T}))

# ---------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------

n_statevars(::Val{2}) = 4
n_statevars(::Val{3}) = 7

# Copy state -> vector (mutable vector, not SVector)
function StateToVector!(v::AbstractVector, s::GeneralPlasticState{2})
    @assert length(v) == 4
    v[1] = s.λ
    v[2] = s.εp.data[1]
    v[3] = s.εp.data[2]
    v[4] = s.εp.data[3]
    return v
end

function StateToVector!(v::AbstractVector, s::GeneralPlasticState{3})
    @assert length(v) == 7
    v[1] = s.λ
    @inbounds for i in 1:6
        v[i+1] = s.εp.data[i]
    end
    return v
end

# Reconstruct state from vector.
# IMPORTANT for ForwardDiff:
# use S = eltype(v), NOT the storage type of the material/operator.
function VectorToState(::Val{2}, v::AbstractVector)
    @assert length(v) == 4
    S = eltype(v)
    λ = v[1]
    εp = SymTensorValue{2,S,3}(v[2], v[3], v[4])
    return GeneralPlasticState{2,S}(λ, εp)
end

function VectorToState(::Val{3}, v::AbstractVector)
    @assert length(v) == 7
    S = eltype(v)
    λ = v[1]
    εp = SymTensorValue{3,S,6}(v[2], v[3], v[4], v[5], v[6], v[7])
    return GeneralPlasticState{3,S}(λ, εp)
end

# ---------------------------------------------------------------------
# Nonlinear operator for the local constitutive solve
# ---------------------------------------------------------------------

mutable struct GeneralPlasticNLOP{D,T,M,S,E} <: NonlinearOperator
    material::M
    state::S
    ε::E
end

GeneralPlasticNLOP(
    material::GeneralPlastic{D,T},
    state::GeneralPlasticState{D,T},
    ε::SymTensorValue{D,T}
) where {D,T} =
    GeneralPlasticNLOP{D,T,typeof(material),typeof(state),typeof(ε)}(
        material, state, ε
    )

# ---------------------------------------------------------------------
# ForwardDiff-safe derivative helpers
# ---------------------------------------------------------------------

# Gradient of a scalar-valued yield surface f(σ) wrt tensor entries.
# We rebuild σ from a vector of entries so ForwardDiff can differentiate it.

function yield_surface_gradient(f, σ::SymTensorValue{2})
    s = SVector(σ.data)
    g = ForwardDiff.gradient(z -> f(SymTensorValue{2,eltype(z),3}(z...)), s)
    return SymTensorValue{2,eltype(g),3}(g...)
end

function yield_surface_gradient(f, σ::SymTensorValue{3})
    s = SVector(σ.data)
    g = ForwardDiff.gradient(z -> f(SymTensorValue{3,eltype(z),6}(z...)), s)
    return SymTensorValue{3,eltype(g),6}(g...)
end

# Derivative of yieldStress(λ) wrt λ
hardening_derivative(f, λ) = ForwardDiff.derivative(f, λ)

# ---------------------------------------------------------------------
# Gridap nonlinear operator interface (fully qualified methods)
# ---------------------------------------------------------------------

function Gridap.Algebra.allocate_residual(
    op::GeneralPlasticNLOP{D,T},
    x::AbstractVector
) where {D,T}
    return zeros(eltype(x), n_statevars(Val(D)))
end

function Gridap.Algebra.zero_initial_guess(op::GeneralPlasticNLOP{D,T}) where {D,T}
    x = zeros(T, n_statevars(Val(D)))
    StateToVector!(x, op.state)
    return x
end

function Gridap.Algebra.allocate_jacobian(
    op::GeneralPlasticNLOP{D,T},
    x::AbstractVector
) where {D,T}
    n = n_statevars(Val(D))
    return zeros(eltype(x), n, n)
end

function Gridap.Algebra.residual!(
    r::AbstractVector{S},
    op::GeneralPlasticNLOP{D,T},
    x::AbstractVector{S}
) where {D,T,S}

    state = VectorToState(Val(D), x)

    σ = op.material.C_elas ⊙ (op.ε - state.εp)

    ∂f∂σ = yield_surface_gradient(op.material.yieldSurface, σ)

    Δλ = state.λ - op.state.λ
    Rεp = state.εp - op.state.εp - Δλ * ∂f∂σ

    r[1] = op.material.yieldSurface(σ) - op.material.yieldStress(state.λ)
    r[2:end] .= Rεp.data

    return r
end

function Gridap.Algebra.jacobian!(
    J::AbstractMatrix{S},
    op::GeneralPlasticNLOP{D,T},
    x::AbstractVector{S}
) where {D,T,S}
    f!(r, xx) = Gridap.Algebra.residual!(r, op, xx)
    y = similar(x)
    ForwardDiff.jacobian!(J, f!, y, x)
    return J
end

# ---------------------------------------------------------------------
# Material response
# ---------------------------------------------------------------------

function material_response(
    op::GeneralPlasticNLOP{D,T},
    Δt,
    cache,
    extras;
    get_Cep=true
) where {D,T}

    σ_trial = op.material.C_elas ⊙ (op.ε - op.state.εp)
    Φ = op.material.yieldSurface(σ_trial) - op.material.yieldStress(op.state.λ)

    # Elastic step
    if Φ <= zero(T)
        return σ_trial,
               op.material.C_elas,
               GeneralPlasticState{D,T}(op.state.λ, op.state.εp)
    end

    # Plastic step: solve local nonlinear system
    nls = NLSolver(show_trace=false, method=:newton)

    x0 = Gridap.Algebra.zero_initial_guess(op)
    cache = solve!(x0, nls, op, cache)

    newState = VectorToState(Val(D), x0)
    σ = op.material.C_elas ⊙ (op.ε - newState.εp)

    if get_Cep
        ∂f∂σ = yield_surface_gradient(op.material.yieldSurface, σ)
        H = hardening_derivative(op.material.yieldStress, newState.λ)

        C_f_σ = op.material.C_elas ⊙ ∂f∂σ
        f_σ_C = ∂f∂σ ⊙ op.material.C_elas

        C_ep = op.material.C_elas - (C_f_σ ⊗ f_σ_C) / (H + ∂f∂σ ⊙ C_f_σ)

        return σ, C_ep, newState
    else
        return σ, newState
    end
end