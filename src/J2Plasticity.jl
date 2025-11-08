struct J2Plasticity{D, T} <: AbstractMaterial
    G::T  # Shear modulus
    K::T  # Bulk modulus
    σ0::T # Initial yield limit
    H::T  # Hardening modulus
    C_elas::SymFourthOrderTensorValue{D, T}  # Elastic stiffness tensor
    function J2Plasticity{D, T}(E::T, ν::T, σ0::T, H::T) where{D, T}
        G = E/(2.0*(1.0+ν))
        K = E/(3.0*(1.0-2.0*ν))
        C_elas = elastic_tangent(D, E, ν)
        return new{D, T}(G, K, σ0, H, C_elas)
    end
end

# keyword argument constructor
J2Plasticity{D, T}(; E, ν, σ0, H) where{D, T} = J2Plasticity{D, T}(E, ν, σ0, H) 

struct J2PlasticityState{D, T} <: AbstractMaterialState
    ϵp :: SymTensorValue{D, T}
    k :: T
end

function initial_material_state(::J2Plasticity{D, T}) where{D, T}
    J2PlasticityState(zero(SymTensorValue{D, T}), T(0))
end

function material_response(
    material::J2Plasticity, ϵ::SymTensorValue{D, T}, state::J2PlasticityState{D, T}, Δt, cache, extras) where {D, T}
    ## unpack some material parameters
    G = material.G
    H = material.H

    ## We use (•)ᵗ to denote *trial*-values
    σᵗ = material.C_elas ⊙ (ϵ - state.ϵp) # trial-stress
    sᵗ = dev(σᵗ)         # deviatoric part of trial-stress
    J₂ = 0.5 * sᵗ ⊙ sᵗ  # second invariant of sᵗ
    σᵗₑ = sqrt(3.0*J₂)   # effective trial-stress (von Mises stress)
    σʸ = material.σ0 + H * state.k # Previous yield limit

    φᵗ  = σᵗₑ - σʸ # Trial-value of the yield surface

    if φᵗ < 0.0 # elastic loading
        return σᵗ, material.C_elas, state
    else # plastic loading
        h = H + 3G
        μ =  φᵗ / h   # plastic multiplier

        c1 = 1 - 3G * μ / σᵗₑ
        s = c1 * sᵗ           # updated deviatoric stress
        σ = s + vol(σᵗ)       # updated stress

        ## Compute algorithmic tangent stiffness ``D = \frac{\Delta \sigma }{\Delta \epsilon}``
        κ = H * (state.k + μ) # drag stress
        σₑ = material.σ0 + κ  # updated yield surface

        # δ(i,j) = i == j ? 1.0 : 0.0
        # Isymdev(i,j,k,l)  = 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) - 1.0/3.0*δ(i,j)*δ(k,l)
        # Q(i,j,k,l) = Isymdev(i,j,k,l) - 3.0 / (2.0*σₑ^2) * s[i,j]*s[k,l]
        # b = (3G*μ/σₑ) / (1.0 + 3G*μ/σₑ)

        # Dtemp(i,j,k,l) = -2G*b * Q(i,j,k,l) - 9G^2 / (h*σₑ^2) * s[i,j]*s[k,l]
        # C_ep = material.C_elas + SymmetricTensor{4, 3}(Dtemp)
        n = 1.5 * s / sqrt(s ⊙ s)
        Celas_n = material.C_elas ⊙ n
        C_ep = material.C_elas - (Celas_n⊗(n⊙material.C_elas))/(H/(3*G) + n ⊙ Celas_n)

        ## Return new state
        Δϵᵖ = 3/2 * μ / σₑ * s # plastic strain
        ϵp = state.ϵp + Δϵᵖ    # plastic strain
        k = state.k + μ        # hardening variable
        return σ, C_ep, J2PlasticityState(ϵp, k)
    end
end