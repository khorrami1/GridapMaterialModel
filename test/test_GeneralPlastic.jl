# ---------------------------------------------------------------------
# Test set
# ---------------------------------------------------------------------

@testset "test_GeneralPlastic" begin

    E = 69.0e3
    ν = 0.3
    Dim = 3
    T = Float64

    # Assumes you already have this in your codebase
    C_elas = elastic_tangent(Dim, E, ν)

    # Generic hardening law (ForwardDiff-compatible)
    yieldStress(ϵp) = 376.9 * (0.0059 + ϵp)^0.152

    # -------------------------------------------------------------
    # Hill48 yield function
    # -------------------------------------------------------------
    struct Yield_Hill48{T}
        F::T
        G::T
        H::T
        M::T
        N::T
        L::T
    end

    function (f::Yield_Hill48)(S1, S2, S3, S4, S5, S6)
        return sqrt(
            (f.F*S1*S1 + f.G*S4*S4 + f.H*S6*S6) /
            (f.F*f.G + f.F*f.H + f.G*f.H) +
            2*S5*S5/f.L +
            2*S3*S3/f.M +
            2*S2*S2/f.N
        )
    end

    # Tensor overload (important: material passes σ as a tensor)
    function (f::Yield_Hill48)(S::SymTensorValue{3})
        return f(S.data...)
    end

    R0  = 0.84
    R45 = 0.64
    R90 = 1.51

    F = R0 / (R90 * (R0 + 1))
    G = 1 / (1 + R0)
    H = R0 / (1 + R0)

    N = (R0 + R90) * (1 + 2*R45) / (2 * R90 * (1 + R0))
    L = N
    M = N

    yield_Hill48 = Yield_Hill48(F, G, H, M, N, L)

    yieldFunction(S::SymTensorValue{Dim}) = yield_Hill48(S)

    material = GeneralPlastic(C_elas, yieldStress, yieldFunction)

    # -------------------------------------------------------------
    # Test 1: initial material state
    # -------------------------------------------------------------
    state0 = initial_material_state(material)

    @test state0.λ == 0.0
    @test state0.εp == zero(SymTensorValue{Dim,Float64})

    # -------------------------------------------------------------
    # Test 2: purely elastic step
    # -------------------------------------------------------------
    ε_el = SymTensorValue{3,Float64,6}(1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0)
    op_el = GeneralPlasticNLOP(material, state0, ε_el)

    σ_el, Cep_el, state_el = material_response(op_el, 1.0, nothing, nothing; get_Cep=true)

    σ_trial_el = C_elas ⊙ ε_el
    Φ_el = yieldFunction(σ_trial_el) - yieldStress(state0.λ)

    @test Φ_el <= 0.0
    @test σ_el ≈ σ_trial_el atol=1e-10 rtol=1e-10
    @test state_el.λ ≈ 0.0 atol=1e-12
    @test state_el.εp ≈ zero(SymTensorValue{Dim,Float64}) atol=1e-12
    @test Cep_el ≈ C_elas atol=1e-10 rtol=1e-10

    # -------------------------------------------------------------
    # Test 3: plastic step
    # -------------------------------------------------------------
    # Use a clearly large strain so the trial stress exceeds the yield stress
    ε_pl = SymTensorValue{3,Float64,6}(2.0e-2, 0.0, 0.0, 0.0, 0.0, 0.0)
    op_pl = GeneralPlasticNLOP(material, state0, ε_pl)

    σ_trial_pl = C_elas ⊙ ε_pl
    Φ_pl = yieldFunction(σ_trial_pl) - yieldStress(state0.λ)

    @test Φ_pl > 0.0

    σ_pl, Cep_pl, state_pl = material_response(op_pl, 1.0, nothing, nothing; get_Cep=true)

    # Plastic multiplier must increase
    @test state_pl.λ > 0.0

    # Plastic strain must be nonzero
    @test norm(state_pl.εp.data) > 0.0

    # Final stress should lie on the yield surface
    f_pl = yieldFunction(σ_pl)
    σy_pl = yieldStress(state_pl.λ)
    @test isapprox(f_pl, σy_pl; atol=1e-8, rtol=1e-8)

    # Return-mapped stress should not exceed trial stress norm
    @test norm(σ_pl.data) <= norm(σ_trial_pl.data) + 1e-8

    # Algorithmic tangent should have same type as elastic tangent
    @test typeof(Cep_pl) == typeof(C_elas)

    # -------------------------------------------------------------
    # Test 4: state <-> vector conversion
    # -------------------------------------------------------------
    x = zeros(Float64, 7)
    StateToVector!(x, state_pl)
    state_chk = VectorToState(Val(3), x)

    @test state_chk.λ ≈ state_pl.λ atol=1e-12
    @test state_chk.εp ≈ state_pl.εp atol=1e-12

    # -------------------------------------------------------------
    # Test 5: residual at converged solution should be near zero
    # -------------------------------------------------------------
    r = Gridap.Algebra.allocate_residual(op_pl, x)
    Gridap.Algebra.residual!(r, op_pl, x)
    @test norm(r) < 1e-7

    # -------------------------------------------------------------
    # Test 6: ForwardDiff Jacobian can be formed
    # -------------------------------------------------------------
    J = Gridap.Algebra.allocate_jacobian(op_pl, x)
    Gridap.Algebra.jacobian!(J, op_pl, x)
    @test size(J) == (7,7)
    @test all(isfinite, J)

end