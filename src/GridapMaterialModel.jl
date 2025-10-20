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

# exporting items

end
