
function elastic_tangent(D::Int, E::T, ν::T) where {T}
    λ = E * ν / ((1 + ν) * (1 - 2ν))
    μ = E / (2 * (1 + ν))

    # Helper to get index in symmetric tensor basis
    function sym_index(i,j)
        i ≤ j ? (i,j) : (j,i)
    end

    # Generate all symmetric index pairs
    idxs = [(i,j) for i in 1:D for j in i:D]
    L = length(idxs)^2
    data = zeros(T, L)

    # Fill data in Voigt-like order
    for (a, (i,j)) in enumerate(idxs)
        for (b, (k,l)) in enumerate(idxs)
            val = λ * (i==j && k==l) + μ * ((i==k && j==l) + (i==l && j==k))
            data[(a-1)*length(idxs) + b] = val
        end
    end

    return SymFourthOrderTensorValue{D,T}(data...)
end

function dev(S)
    S - vol(S)
end

function vol(S)
    1.0/3.0 * tr(S) * one(S)
end

@inline function grad_wrt_entries(f::Function, S::SymTensorValue{2,T}) where {T}
    return T(FD.derivative(S1->f(S1, S[2], S[3]), S[1]),
             FD.derivative(S2->f(S[1], S2, S[3]), S[2]),
             FD.derivative(S3->f(S[1], S[2], S3), S[3]))
end

@inline function grad_wrt_entries(f::Function, S::SymTensorValue{3,T}) where {T}
    return T(FD.derivative(S1->f(S1, S[2], S[3], S[4], S[5], S[6]), S[1]),
             FD.derivative(S2->f(S[1], S2, S[3], S[4], S[5], S[6]), S[2]),
             FD.derivative(S3->f(S[1], S[2], S3, S[4], S[5], S[6]), S[3]),
             FD.derivative(S4->f(S[1], S[2], S[3], S4, S[5], S[6]), S[4]),
             FD.derivative(S5->f(S[1], S[2], S[3], S[4], S5, S[6]), S[5]),
             FD.derivative(S6->f(S[1], S[2], S[3], S[4], S[5], S6), S[6]))
end