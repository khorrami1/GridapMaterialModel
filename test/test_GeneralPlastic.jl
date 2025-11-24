
@testset "test_GeneralPlastic" begin
    
    E = 69.0e3
    ν = 0.3
    Dim = 3
    C_elas = elastic_tangent(Dim, E, ν)
    yieldStress(ϵ) = 376.9*(0.0059 + ϵ)^0.152

    struct Yield_Hill48{T}
        F::T
        G::T
        H::T
        M::T
        N::T
        L::T
    end 

    
    function (f::Yield_Hill48)(S1, S2, S3, S4, S5, S6)
        return sqrt( (f.F*S1*S1 + f.G*S4*S4 + f.H*S6*S6) / 
            (f.F*f.G + f.F*f.H + f.G*f.H) + 2*S5*S5/f.L + 2*S3*S3/f.M + 2*S2*S2/f.N )
    end

    R0 =  0.84
    R45 = 0.64
    R90 = 1.51

    F = R0/(R90*(R0+1))
    G = 1/(1+R0)
    H = R0/(1+R0)

    N = (R0+R90)*(1+2*R45)/(2*R90*(1+R0))
    L = N # must be checked!
    M = N # must be checked!

    yield_Hill48 = Yield_Hill48(F, G, H, M, N, L)
    # yieldFunction1(S::SymTensorValue{Dim,T}) = yield_Hill48(S)
    yieldFunction(S1, S2, S3, S4, S5, S6) = yield_Hill48(S1, S2, S3, S4, S5, S6)

    material = GeneralPlastic{Dim, Float64}(C_elas, yieldStress, yieldFunction)

    

end