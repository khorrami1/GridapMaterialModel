
kroneckarDelta(i, j) = i==j ? 1 : 0

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