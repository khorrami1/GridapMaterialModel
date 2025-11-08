module GridapMaterialModel

using Gridap.TensorValues
using NLsolve
using Rotations
import ForwardDiff
using StaticArrays

abstract type AbstractMaterial end

abstract type AbstractMaterialState end

abstract type AbstractResiduals end

abstract type StrainMeasure end

abstract type AbstractCache end

function material_response end

function initial_material_state end

function get_cache(::AbstractMaterial)
    nothing
end

function update_cache! end

# including files
include("wrappers.jl")
include("LinearElastic.jl")
include("Plastic.jl")

# exporting items
export initial_material_state, get_cache, material_response
export elastic_tangent
export LinearElastic, LinearElasticState

end
