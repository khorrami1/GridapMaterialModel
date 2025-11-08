using GridapMaterialModel
using Test
using Gridap.TensorValues
using NLsolve
using Rotations
import ForwardDiff
using StaticArrays

include("test_wrappers.jl")
include("test_LinearElastic.jl")
include("test_J2Plasticity.jl")