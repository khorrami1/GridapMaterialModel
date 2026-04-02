using Gridap
using Gridap.Algebra
using ForwardDiff

struct NLOpTest <: NonlinearOperator 
    cfg::ForwardDiff.JacobianConfig
    history::Float64
end

function Gridap.Algebra.residual!(r::AbstractVector, op::NLOpTest, x::AbstractVector)
    r[1] = (x[1] - 2) * (x[1] - 1)
    r[2] = x[2] - op.history
    return r
end

function Gridap.Algebra.allocate_residual(op::NLOpTest, x::AbstractVector)
    return zeros(eltype(x), 2)
end

function Gridap.Algebra.zero_initial_guess(op::NLOpTest)
    return zeros(Float64, 2)
end

function Gridap.Algebra.allocate_jacobian(op::NLOpTest, x::AbstractVector)
    return zeros(eltype(x), 2, 2)
end

# function Gridap.Algebra.jacobian!(A::AbstractMatrix, op::NLOpTest, x::AbstractVector)
#     A[1,1] = (x[1] - 1) + (x[1] - 2)   # d/dx of (x-2)(x-1)
#     A[1,2] = 0.0
#     A[2,1] = 0.0
#     A[2,2] = 1.0
#     return A
# end


function Gridap.Algebra.jacobian!(A::AbstractMatrix, op::NLOpTest, x::AbstractVector)
    # f!(r, xx) = Gridap.Algebra.residual!(r, op, xx)
    ForwardDiff.jacobian!(A, f!, x, similar(x), op.cfg)
    return A
end

f!(r, xx) = Gridap.Algebra.residual!(r, op, xx)

cfg = ForwardDiff.JacobianConfig(f!, ones(2), zeros(2),
    ForwardDiff.Chunk{1}(), nothing)

op = NLOpTest(cfg, 3.0)
nls = NLSolver(show_trace=true, method=:newton)

x0 = Gridap.Algebra.zero_initial_guess(op)
cache = nothing
@time cache = solve!(x0, nls, op, cache)

println("solution = ", x0)

res = similar(x0)

@time op = NLOpTest(cfg, 3.2)
@time cache = solve!(x0, nls, op, cache)
println("solution = ", x0)

